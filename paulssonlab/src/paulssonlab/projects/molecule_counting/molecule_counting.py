import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.integrate import simps, trapz
from scipy.optimize import curve_fit
import seaborn as sns
from sklearn.linear_model import LinearRegression


###These functions are all involved in the initial processing of the time series data.
def decay_rates(data, t_end):
    return data.apply(
        lambda y: (y.values[0] - y.values[t_end - 1])
        / np.sum(y[:t_end].values - y.values[t_end - 1]),
        axis=1,
    ).values


def filter_traces(data, ub, lb, dt=1, t_max=800):
    data["dy_max"] = data.apply(
        lambda y: (y.values[: t_max - dt] / y.values[dt:t_max]).max(), axis=1
    )
    data["dy_min"] = data.apply(
        lambda y: (y.values[: t_max - dt] / y.values[dt:t_max]).min(), axis=1
    )
    return data[(data.dy_max < ub) & (data.dy_min > lb)]


def process_df(full_data, fp_dir, bad_fovs=[], trunc=0, dy_max=10, dy_min=0):
    traces = []
    rp_list = []
    rp_df = pd.DataFrame()

    for i in full_data[fp_dir].keys():
        if i not in bad_fovs:
            print(i)

            regionprops = full_data[fp_dir][i]["regionprops"]
            measurements = full_data[fp_dir][i]["traces"]
            area = full_data[fp_dir][i]["regionprops"]["area"]
            y = measurements["mean"][1:, 1:]

            traces.append(y)
            rp_list.append(regionprops)
            print(measurements["mean"].shape)

    traces = np.concatenate(traces)
    rp_df = pd.concat(rp_list, sort=False)
    rp_df.reset_index(inplace=True)
    print(traces.shape)

    data = pd.DataFrame(traces)
    data = pd.concat([data, rp_df], axis=1, sort=False)

    t_end = traces.shape[1] - trunc
    data["k"] = decay_rates(data, t_end)
    data = filter_traces(data, dy_max, dy_min, 10, 950)

    print(data.shape)
    plt.figure()
    data.k.hist(bins=30)
    plt.title("Histogram of Photobleaching Rates")

    plt.figure()
    data[0].hist(bins=30)
    plt.title("Histogram of Initial Intensities")
    return data


###These function are all involved in statistical inference and plotting
def sum_exp(t, *p):
    if isinstance(p[0], (tuple, list, np.ndarray)):
        p = p[0]

    n = int((len(p) - 1) / 2)
    y = p[0]
    C = p[1 : n + 1]
    k = p[n + 1 :]

    for i in range(n):
        y += C[i] * np.exp(-k[i] * t)
    return y


def filter_df(df, prop_dict):
    processed_df = df.copy()
    for prop in prop_dict:
        processed_df = processed_df[
            (processed_df[prop] > prop_dict[prop][0])
            & (processed_df[prop] < prop_dict[prop][1])
        ]
    processed_df = processed_df
    return processed_df


def get_stats(df, n):
    y = df.filter(regex="^[0-9]+$", axis=1)
    #     y = y.rolling(window=1).mean()
    mu = y.mean()

    tp = np.arange(0, mu.size)
    lb = (n + 1) * [-1] + n * [0]
    ub = (n + 1) * [1] + n * [0.1]
    params, _ = curve_fit(
        sum_exp,
        tp,
        mu / mu[0],
        p0=0.01 * np.ones(2 * n + 1),
        bounds=(lb, ub),
        maxfev=5000,
    )

    bg = mu[0] * params[0]
    C = mu[0] * params[1 : n + 1]
    k = params[n + 1 :]
    k, C = np.split(np.array(sorted(zip(k, C), reverse=True)), 2, axis=1)

    plt.figure(figsize=(12, 10))
    plt.semilogy(tp / 10, mu / mu[0], linewidth=10, alpha=0.7, label="Empirical Mean")
    plt.semilogy(
        tp / 10,
        sum_exp(tp, params),
        linewidth=10,
        alpha=0.7,
        label=f"{n} Exponential Fit",
    )
    plt.xlabel("Time", fontsize=28)
    plt.ylabel("Fraction Photobleached", fontsize=28)
    plt.legend(fontsize=24)
    #     plt.ylim(bottom=1000)
    #     plt.xlim(left=0,right=100)
    print(f"background intensity is {bg:.0f}")
    mu -= bg
    y -= bg

    print(y.shape)
    p = (mu / mu[0]).values
    I0 = y.iloc[:, :1].mean(axis=1).values
    ybar = np.outer(p, I0.T).T
    sigma2 = (y - ybar) ** 2
    return bg, C, k, y, mu, sigma2


def nstep_params(k):
    n = k.size
    K = np.zeros((n, n))
    for i in range(1, n + 1):
        K[i - 1, :i] = np.array(
            [
                np.prod(k[: i - 1]) / np.prod(np.delete(k[:i], j) - k[j])
                for j in range(i)
            ]
        )
    return K


def nu_fit_seq(tp, C, k, fbar):
    K = nstep_params(k)
    E = np.exp(-np.outer(k, tp))
    P = np.dot(K, E)

    C_prime = np.dot(C.T, np.linalg.inv(K))
    C_bar = np.repeat(C_prime.T, tp.size, axis=1) - np.dot(C_prime, P)

    X = (C_bar * P).T
    y = fbar
    reg = LinearRegression(fit_intercept=True).fit(X, y)
    print(f"Regression score is {reg.score(X,y):.3f}")
    nu = reg.coef_

    #     nu = [10,4]
    f_bg = reg.intercept_
    f_pred = reg.predict(X) - f_bg
    fbar[1:] = fbar[1:] - f_bg
    V = np.multiply(nu, X)
    return f_pred, nu, V, P


def nu_fit(tp, C, k, fbar):
    P = np.exp(np.outer(-1 * k, tp))
    X = np.multiply(C, P * (1 - P)).T
    y = fbar
    print(X.shape, y.shape)
    reg = LinearRegression(fit_intercept=True).fit(X, y)
    nu = reg.coef_
    f_bg = reg.intercept_
    f_pred = reg.predict(X) - f_bg
    fbar[1:] = fbar[1:] - f_bg
    V = np.multiply(nu, X)
    return f_pred, nu, V, P


def pb_inference(C, k, mu, sigma2, q=1):
    cq = -1 / (1 / 2 * q**2 - 1 / 3 * q**3)
    tp = np.arange(0, mu.size)
    n = len(C)

    pbar = mu / mu[0]
    f = pd.DataFrame(sigma2.values)
    fbar = f.mean(axis=0)
    #     fbar =  fbar - .8*fbar.iloc[-50:].min()
    #     fbar[1:] =  fbar[1:] - (fbar.iloc[1]*pbar**2)[1:]
    #     fbar[1:] =  fbar[1:] - (fbar[1]*pbar**2)[1:]

    f_predicted, nu, var_list, p_list = nu_fit_seq(tp, C, k, fbar)

    for i in range(n):
        print(
            f"nu_{i} = {nu[i]:.2f}, C_{i} = {C[i][0]:.0f}, k_{i} = {k[i][0]:.4f} inferred molecule count = {C[i][0]/nu[i]:.0f}"
        )
    return nu, f_predicted, var_list, p_list, fbar, tp


def fluct_plot(df, n=1, q=1):
    cmap = cm.get_cmap("coolwarm")
    bg, C, k, y, mu, sigma2 = get_stats(df, n)

    nu, f_predicted, var_list, p_list, fbar, tp = pb_inference(C, k, mu, sigma2)
    pbar = mu / mu[0]
    plt.figure(figsize=(12, 10))
    plt.scatter(1 - pbar, fbar, color="k", label="Single-cell data")

    print(simps(fbar, 1 - pbar), simps(f_predicted, 1 - pbar))
    plt.plot(
        1 - pbar,
        f_predicted,
        linewidth=6,
        alpha=0.7,
        color="g",
        label="Theoretical curve",
    )

    plt.fill_between(
        1 - pbar, f_predicted, alpha=0.3, color="c", label="Integrated area"
    )

    plt.xlabel(r"Fraction Photobleached", fontsize=28)
    plt.ylabel(r"$\hat{\sigma}^2$", fontsize=28)
    plt.ylim(bottom=0, top=2500)
    plt.xlim(left=0, right=1)
    plt.legend(fontsize=20)

    plt.figure(figsize=(12, 8))
    plt.plot(tp, f_predicted, linewidth=6, color="g")
    plt.plot(tp, fbar, linewidth=6, alpha=0.7, color="b")
    for i in range(n):
        p = p_list[i]
        plt.plot(tp, var_list, linewidth=6, alpha=0.5, label="Theoretical curve")
    #         plt.fill_between(1-p, var,alpha=.3, label='Integrated area')
    return f_predicted, fbar, pbar, mu, sigma2, y


def kmaps(df):
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(24, 12))

    cmap = cm.get_cmap("coolwarm")

    fit_df = df[["centroid_x", "centroid_y", "k", 0]].copy()
    fit_df["c_k"] = (fit_df.k - fit_df.k.min()) / (fit_df.k.max() - fit_df.k.min())
    fit_df["c_I0"] = (fit_df[0] - fit_df[0].min()) / (fit_df[0].max() - fit_df[0].min())

    sns.regplot(
        data=fit_df,
        x="centroid_y",
        y="centroid_x",
        fit_reg=False,
        ax=axs[0],
        scatter_kws={"color": cmap(fit_df["c_k"]), "alpha": 0.5},
    )

    sns.regplot(
        data=fit_df,
        x="centroid_y",
        y="centroid_x",
        fit_reg=False,
        ax=axs[1],
        scatter_kws={"color": cmap(fit_df["c_I0"]), "alpha": 0.5},
    )

    axs[0].set_xlim([0, 2000])
    axs[0].set_ylim([0, 2000])
    axs[1].set_xlim([0, 2000])
    axs[1].set_ylim([0, 2000])
