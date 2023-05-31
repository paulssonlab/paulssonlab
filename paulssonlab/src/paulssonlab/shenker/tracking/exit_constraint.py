import cv2
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import skimage
from gurobipy import *
from PIL import Image
from skimage import morphology

# from line_profiler import LineProfiler


# @profile
def main2(mainhuh):
    # @profile
    def get_coordinates(img1, mask_num):
        output2 = np.where(img1 == mask_num)
        mask_1 = [mask_num for j in range(len(output2[0]))]
        listOfCoordinates = list(zip(output2[0], output2[1], mask_1))

        return listOfCoordinates

    # @profile
    def get_array(n_mask_img, imgs):
        arr_array = [[]] * n_mask_img
        for i in range(n_mask_img):
            arr1 = get_coordinates(
                imgs, i + 1
            )  # list(zip(np.where(imgs==i+1)[0], np.where(imgs==i+1)[1], [i+1 for j in range (len(np.where(imgs==i+1)[0]))]))

            arr_array[i] = arr1
        return arr_array

    def get_tree(n_mask_img1, arr_array, n_mask_img1a, arra_array):
        tree = []

        for i in range(0, n_mask_img1):
            for k in range(len(arr_array[i])):
                for l in range(0, n_mask_img1a):
                    if (
                        arra_array[l][0][0] == arr_array[i][k][0]
                        and arra_array[l][0][1] == arr_array[i][k][1]
                    ):
                        tree += [[l + 1, i + 1]]
        res = []
        for k in tree:
            if k not in res:
                res.append(k)
        return res

    def get_center(threshold_number, mask_number):
        mask_number = mask_number - 1
        array1 = arr[mn][threshold_number][mask_number]
        max_y = max(array1, key=lambda item: item[0])
        min_y = min(array1, key=lambda item: item[0])
        max_x = max(array1, key=lambda item: item[1])
        min_x = min(array1, key=lambda item: item[1])

        arr_center = [(max_y[0] + min_y[0]) / 2, (max_x[1] + min_x[1]) / 2]
        arr_center_x = arr_center[1]

        arr_lengths = [max_y[0], min_y[0], max_x[1], min_x[1]]
        arr_hist = img_ini[
            int((arr_lengths[0] + arr_lengths[1]) / 5),
            arr_lengths[3] : arr_lengths[2] + 1,
        ]
        arr_hist2 = img_ini[
            int((arr_lengths[0] + arr_lengths[1]) / 1.25),
            arr_lengths[3] : arr_lengths[2] + 1,
        ]
        number_of_peaks = scipy.signal.find_peaks(arr_hist[:, 0])
        number_of_peaks2 = scipy.signal.find_peaks(arr_hist2[:, 0])

        # pp = len(number_of_peaks[0])

        if len(number_of_peaks[0]) == len(number_of_peaks2[0]):
            pp = (len(number_of_peaks[0]) + len(number_of_peaks2[0])) / 2
        else:
            pp = len(number_of_peaks[0]) + 1  # len(number_of_peaks2[0]))

        len_x = max_x[1] - min_x[1]
        right_edge = max_x[1]
        return arr_center_x, pp, len_x, right_edge

    def tree2(x):
        tree2_1 = get_tree(n_mask_img[mn][0], arr[mn][0], n_mask_img[mn][3], arr[mn][3])
        tree2_2 = get_tree(n_mask_img[mn][1], arr[mn][1], n_mask_img[mn][3], arr[mn][3])
        tree2_3 = get_tree(n_mask_img[mn][2], arr[mn][2], n_mask_img[mn][3], arr[mn][3])

        return tree2_1, tree2_2, tree2_3

    arr_m = [0] * mainhuh
    df_list = [pd.DataFrame()] * len(arr_m)
    s = [[7] * 1000] * len(arr_m)
    s2 = [[]] * len(arr_m)
    a_exit = [[7] * 1000] * len(arr_m)
    a_div2 = [[7] * 1000] * len(arr_m)
    a_div3 = [[7] * 1000] * len(arr_m)
    imgs_mat = [[[]]] * len(arr_m)
    th_mat = [[[]]] * len(arr_m)
    n_constr = [0] * len(arr_m)
    n_constr_a = [0] * len(arr_m)
    n_constr_b = [0] * len(arr_m)
    n_constr_c = [0] * len(arr_m)
    n_mask_img = [[]] * len(arr_m)
    model_tree = Model()
    arr = [[[]]] * len(arr_m)
    len_df = [0] * len(arr_m)
    for mn in range(len(arr_m)):
        img_ini = cv2.imread("%d.jpg" % (mn + 0))  # +29+2)) #12
        # cv2.imshow('img_ini', img_ini)
        # cv2.waitKey(0)
        # cv2.destroyAllWindows()
        th = [[]] * 4
        imgs = [[]] * 4
        reta, th[0] = cv2.threshold(
            img_ini, 150, 255, cv2.THRESH_BINARY
        )  # 160#187#188#180#####180
        ret, th[1] = cv2.threshold(
            img_ini, 204, 255, cv2.THRESH_BINARY
        )  # 180#190#195#185###190
        reta, th[2] = cv2.threshold(
            img_ini, 200, 255, cv2.THRESH_BINARY
        )  # 200#205#205#######200
        reta, th[3] = cv2.threshold(
            img_ini, 207, 255, cv2.THRESH_BINARY
        )  # 220#215#212#######208
        n_mask_img[mn] = [0] * (len(imgs))
        arr[mn] = [[]] * (len(imgs))
        costs_get_arr2 = []
        center_arr = []
        lenx_arr = []
        peaks = []
        r_e = []
        df = pd.DataFrame()
        for i in range(len(th)):
            img = skimage.morphology.label(
                th[i],
                neighbors=None,
                background=None,
                return_num=False,
                connectivity=None,
            )
            imgs[i] = img
            n_mask_img[mn][i] = np.amax(imgs[i])
            arr2 = get_array(n_mask_img[mn][i], imgs[i])
            arr[mn][i] = arr2
            for j in range(1, (n_mask_img[mn][i]) + 1):
                cc = get_center(i, j)
                costs_get_arr2 += [[i, j]]
                center_arr += [cc[0]]
                lenx_arr += [cc[2]]
                peaks += [cc[1]]
                r_e += [cc[3]]

        df["threshold_number"] = [x[0] for x in costs_get_arr2]
        df["mask"] = [x[1] for x in costs_get_arr2]
        df_list[mn] = df
        df_list[mn]["center"] = center_arr
        df_list[mn]["length_x"] = lenx_arr
        df_list[mn]["t"] = mn
        df_list[mn]["peaks"] = peaks
        df_list[mn]["right_edge"] = r_e

        matrix2 = []
        a_div2[mn - 1] = [[0] * len(df_list[mn])] * len(df_list[mn - 1])
        cost_move = [[0] * len(df_list[mn])] * len(df_list[mn - 1])

        for j in range(len(df_list[mn - 1])):
            row = [0] * len(df_list[mn])
            for i in range(len(df_list[mn])):
                cost_move[j][i] = (
                    0.033 * (abs(center_arr[i] - df_list[mn - 1]["center"].values[j]))
                    + 0.050 * abs(peaks[i] + df_list[mn - 1]["peaks"].values[j] - 2)
                    + +(0.001 * df_list[mn - 1]["length_x"].values[j])
                    + 0.036 * abs(lenx_arr[i] - df_list[mn - 1]["length_x"].values[j])
                )

                if (
                    lenx_arr[i] < 10 or df_list[mn - 1]["length_x"].values[j] < 10
                ):  # or center_arr[i] < df_list[mn-1]['center'].values[j]-10):# > 100):
                    # print(cost_move[j][i],center_arr[i] )
                    cost_move[j][i] = cost_move[j][i] + 10000
                    # print(cost_move[j][i])
                else:
                    cost_move[j][i] = cost_move[j][i]
                xm = model_tree.addVar(
                    0,
                    1,
                    obj=cost_move[j][i],
                    vtype=GRB.BINARY,
                    name="x_%d_%d_%d" % (mn - 1, j, i),
                )

                row[i] = xm
            matrix2 = row
            a_div2[mn - 1][j] = matrix2

        matrix3 = []

        if mn > 0:
            a_div3[mn - 1] = [[[0] * len(df_list[mn])] * len(df_list[mn])] * len(
                df_list[mn - 1]
            )

            cost_divide = (
                0.340
                * abs(
                    df_list[mn]["center"].values[np.newaxis, np.newaxis, :]
                    - df_list[mn - 1]["center"].values[:, np.newaxis, np.newaxis]
                )
                + +0.004
                * (
                    df_list[mn]["peaks"].values[np.newaxis, :, np.newaxis]
                    + df_list[mn]["peaks"].values[np.newaxis, np.newaxis, :]
                    + df_list[mn - 1]["peaks"].values[:, np.newaxis, np.newaxis]
                    - 3
                )
                + +0.0024
                * abs(
                    df_list[mn]["center"].values[np.newaxis, np.newaxis, :]
                    - df_list[mn]["center"].values[np.newaxis, :, np.newaxis]
                )
                + +(-0.0012 * df_list[mn - 1]["length_x"].values[j])
                + +0.0004
                * (
                    df_list[mn]["length_x"].values[np.newaxis, :, np.newaxis]
                    - df_list[mn - 1]["length_x"].values[:, np.newaxis, np.newaxis]
                )
            )  # - 1.2*df_list[mn-1]['length_x'].values[:,np.newaxis,np.newaxis]#+\

            for k in range(len(df_list[mn - 1])):
                row = [0] * len(df_list[mn])

                c_i = [0] * len(df_list[mn])

                for i in range(len(df_list[mn])):
                    ro2 = [0] * len(df_list[mn])
                    c_i[i] = center_arr[i]

                    if c_i[i] != max(center_arr):
                        index_i = [
                            i2
                            for i2 in range(len(df_list[mn]))
                            if center_arr[i2] > c_i[i]
                        ]
                        for j in index_i:  #
                            if (
                                df_list[mn]["length_x"].values[i] < 10
                                or df_list[mn - 1]["length_x"].values[k] < 10
                                or df_list[mn]["length_x"].values[j] < 10
                            ):  # or abs(df_list[mn]['center'].values[i]   - (df_list[mn-1]['center'].values[k])) >100  or (df_list[mn]['length_x'].values[i]   > (df_list[mn-1]['length_x'].values[k])) ):
                                cost_divide[k, i, j] = cost_divide[k, i, j] + 10000000
                            else:
                                cost_divide[k, i, j] = cost_divide[k, i, j]

                            xm = model_tree.addVar(
                                0,
                                1,
                                obj=cost_divide[k, i, j],
                                vtype=GRB.BINARY,
                                name="x_%d_%d_%d_%d" % (mn - 1, k, i, j),
                            )  # name="x_{0}_{1}_{2}_{3}".format(mn-1,k,i,j))

                            ro2[j] = xm
                    row[i] = ro2
                matrix3 = row
                a_div3[mn - 1][k] = matrix3

            a_exit[mn - 1] = [0] * len(df_list[mn - 1])
            for k in range(0, len(df_list[mn - 1])):
                x_exit = 0
                if df_list[mn - 1]["length_x"].values[k] < 10:
                    x_exit = 100000
                a1 = model_tree.addVar(
                    0,
                    1,
                    obj=(1 + x_exit),
                    vtype=GRB.BINARY,
                    name="a_exit%d,%d" % (mn - 1, k),
                )  ### 20 for it to be not chosen
                a_exit[mn - 1][k] = a1

        ####################################################################################################
        lmn_1 = len(df_list[mn - 1])
        lmn = len(df_list[mn])

        if mn > 1:
            s[mn - 1] = [7] * 100
            for i in range(0, lmn_1):
                aac2ab = quicksum(
                    a_div2[mn - 2][k][i]
                    for k in range(0, len(df_list[mn - 2]))
                    if type(a_div2[mn - 2][k][i]) != int
                )
                aac3ab = quicksum(
                    a_div3[mn - 2][k][i][j]
                    for j in range(0, lmn_1)
                    for k in range(0, len(df_list[mn - 2]))
                    if type(a_div3[mn - 2][k][i][j]) != int
                )
                aac5ab = quicksum(
                    a_div3[mn - 2][k][j][i]
                    for j in range(0, lmn_1)
                    for k in range(0, len(df_list[mn - 2]))
                    if type(a_div3[mn - 2][k][j][i]) != int
                )
                s1 = (
                    quicksum(
                        a_div3[mn - 1][i][j][k]
                        for j in range(0, lmn)
                        for k in range(0, lmn)
                        if type(a_div3[mn - 1][i][j][k]) != int
                    )
                    + quicksum(
                        a_div2[mn - 1][i][k]
                        for k in range(0, lmn)
                        if type(a_div2[mn - 1][i][k]) != int
                    )
                    + a_exit[mn - 1][i]
                )
                con2a = model_tree.addConstr(aac2ab + aac3ab + aac5ab - s1 == 0)
                s[mn - 1][i] = s1

        if mn == 1:
            s[mn - 1] = [7] * 10000000
            for i in range(0, len(df_list[mn - 1])):
                s1 = (
                    quicksum(
                        a_div3[mn - 1][i][j][k]
                        for j in range(0, lmn)
                        for k in range(0, lmn)
                        if type(a_div3[mn - 1][i][j][k]) != int
                    )
                    + quicksum(
                        a_div2[mn - 1][i][k]
                        for k in range(0, lmn)
                        if type(a_div2[mn - 1][i][k]) != int
                    )
                    + a_exit[mn - 1][i]
                )

                s[mn - 1][i] = s1

        if mn > 0:
            n_constr[mn - 1] = n_mask_img[mn - 1][len(th_mat[mn - 1]) - 1]
            n_constr_c[mn - 1] = (
                n_mask_img[mn - 1][len(th_mat[mn - 1]) - 2]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 3]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 4]
                - 1
            )
            n_constr_a[mn - 1] = (
                n_mask_img[mn - 1][len(th_mat[mn - 1]) - 3]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 4]
                - 1
            )
            n_constr_b[mn - 1] = n_mask_img[mn - 1][len(th_mat[mn - 1]) - 4] - 1

        if mn > 0:
            mn = mn - 1
            # tr = get_tree(0)
            tr = tree2(0)
            tree2_1 = tr[0]  # tree2(0,3)
            tree2_2 = tr[1]  # tree2(1,3)
            tree2_3 = tr[2]  # tree2(2,3)
            for i in range(1, n_constr[mn] + 1):
                con1 = model_tree.addConstr(
                    s[mn][
                        (tree2_1[int(np.where([x[0] == i for x in tree2_1])[0])][1]) - 1
                    ]
                    + s[mn][
                        (n_constr_b[mn])
                        + (tree2_2[int(np.where([x[0] == i for x in tree2_2])[0])][1])
                    ]
                    + s[mn][
                        n_constr_a[mn]
                        + (tree2_3[int(np.where([x[0] == i for x in tree2_3])[0])][1])
                    ]
                    + s[mn][n_constr_c[mn] + i]
                    == 1
                )

            mn = mn + 1

        ######EXIT  CONSTRAINT #####

        if mn > 0:
            c_a = [0] * len(df_list[mn - 1])

            index_i = []
            for i in range(len(df_list[mn - 1])):
                sum_see2 = LinExpr()
                c_a[i] = df_list[mn - 1]["right_edge"].values[i]
                index_i = [
                    i2a
                    for i2a in range(len(df_list[mn - 1]))
                    if (df_list[mn - 1]["center"].values[i2a] > c_a[i])
                ]  # and i2>k  ]
                sum_see2 = sum(
                    quicksum(
                        a_div3[mn - 1][ia][j][k]
                        for j in range(0, lmn)
                        for k in range(0, lmn)
                        if type(a_div3[mn - 1][ia][j][k]) != int
                    )
                    + quicksum(
                        a_div2[mn - 1][ia][k]
                        for k in range(0, lmn)
                        if type(a_div2[mn - 1][ia][k]) != int
                    )
                    for ia in index_i
                )
                if type(sum_see2) != int:
                    k = len(df_list[mn - 1])
                    con_exit = model_tree.addConstr(
                        ((k * a_exit[mn - 1][i]) + sum_see2) <= k
                    )

        imgs_mat[mn] = imgs
        th_mat[mn] = th
        len_df[mn] = len(df_list[mn])

    model_tree.update()
    model_tree.Params.Threads = 1
    # model_tree.Params.SiftMethod=1
    model_tree.optimize()
    model_tree.write("model.lp")

    if model_tree.status == GRB.INFEASIBLE:
        vars1 = model_tree.getVars()
        # ubpen = [1.0]*model_tree.numVars
        print(vars1)
        model_tree.feasRelax(1, False, vars1, None, None, model_tree.getConstrs(), None)
        model_tree.optimize()

    asg_arr = []
    for v in model_tree.getVars():
        if v.x == 1:  # (v.varName[0]=='x' and v.x==1):
            asg_arr += [v.varName]

    print(">>", asg_arr)

    asg_arr2a = []
    for i in range(len(asg_arr)):
        if type(asg_arr[i][2]) == int:
            if int(asg_arr[i][2]) == (len(arr_m) - 2):
                asg_arr2a += [asg_arr[i]]

    asg_arr3a = [[]] * len(asg_arr2a)
    for i in range(len(asg_arr2a)):
        if len(asg_arr2a[i].split("_")) == 5:
            asg_arr3a[i] = asg_arr2a[i].split("_")[-2:]
        if len(asg_arr2a[i].split("_")) == 4:
            asg_arr3a[i] = [asg_arr2a[i].split("_")[-1]]

    asg_arr4a = []
    for i in range(0, len(asg_arr3a)):
        for j in range(len(asg_arr3a[i])):
            asg_arr4a += [asg_arr3a[i][j]]

    s_arr = []

    s2 = [[]] * (len(arr_m) - 1)
    for mno in range(0, len(arr_m) - 1):
        s2[mno] = [0] * (len(df_list[mno]))
        for i in range(len(df_list[mno])):
            if type(s[mno][i]) != int:
                s2[mno][i] = LinExpr.getValue(s[mno][i])
            s_arr += [s2[mno][i]]

    make = [0] * len(df_list[len(arr_m) - 1])
    for i in range(0, len(df_list[len(arr_m) - 1])):
        for j in range(0, len(asg_arr4a)):
            if i == j:
                make[i] = 1
            else:
                make[i] = 0

    df_list[len(arr_m) - 1]["s_i"] = make
    df_list_concat = pd.concat([df_list[mn] for mn in range(0, len(arr_m))])
    s_arr2 = s_arr + make
    df_list_concat["s_i"] = s_arr2
    df_list_concat = df_list_concat.reset_index(drop=True)

    def get_img_we_want2(i_df, df_img):
        imgsa = imgs_mat[i_df]
        tha = th_mat[i_df]
        threshold_number = df_img["threshold_number"][0]
        mask_number = df_img["mask"][0]
        img3 = np.where(
            imgsa[threshold_number] == mask_number,
            tha[threshold_number],
            tha[threshold_number] * 0,
        )

        if len(df_img) > 1:
            for i in range(1, len(df_img)):
                threshold_number = df_img["threshold_number"][i]
                mask_number = df_img["mask"][i]
                img3 = img3 + np.where(
                    imgsa[threshold_number] == mask_number,
                    tha[threshold_number],
                    tha[threshold_number] * 0,
                )

        return Image.fromarray(img3)

    for jk in range(0, len(arr_m)):
        df_list[jk] = df_list_concat[df_list_concat["t"] == jk]
        df_new = df_list[jk][df_list[jk]["s_i"] == 1]
        df_new = df_new.reset_index(drop=True)
        if jk < len(arr_m) - 1:
            img_mask_video = get_img_we_want2(jk, df_new)
            img_mask_video.save("%d.png" % jk)

    df_list_concat = df_list_concat[df_list_concat["s_i"] == 1]
    df_list_concat = df_list_concat.reset_index(drop=False)

    index2 = df_list_concat["index"].values

    t1 = df_list_concat["t"].values
    for j in range(len(arr_m)):
        sm = sum(len_df[m] for m in range(0, j))
        for i in range(len(t1)):
            if t1[i] == j:
                index2[i] -= sm
    df_list_concat["index"] = index2

    if len(df_list_concat) != len(asg_arr):
        asg_arr += [0]

    l_asg = len(asg_arr)
    for i in range(l_asg):
        asg_arr[i] = str(asg_arr[i])
        asg_arr[i] = asg_arr[i].split("_")

    df_list_concat["assign"] = asg_arr
    df_list_concat = df_list_concat.reset_index(drop=False)

    Mod = [""] * l_asg
    for i in range(l_asg):
        if len(asg_arr[i]) == 4:
            Mod[i] = "move"
        elif len(asg_arr[i]) == 5:
            Mod[i] = "divide"
        elif len(asg_arr[i]) == 1:
            Mod[i] = "last"
        elif len(asg_arr[i]) == 2:
            Mod[i] = "exit"  ##

    df_list_concat["MoD"] = Mod

    import glob

    images = []
    filenames = [img for img in glob.glob("*.png")]
    filenames.sort()
    int_files = [int(filenames[i][:-4]) for i in range(len(filenames))]
    int_files.sort()
    final_files = ["{}.png".format((str(i))) for i in int_files]

    images = []
    for img in final_files:
        n = cv2.imread(img)
        images.append(n)

    img0 = cv2.imread("0.png")
    height, width, layers = img0.shape

    video = cv2.VideoWriter(
        "video_mask_1a.mp4", -1, 5, (width, height)
    )  # second last is fps should be a factor of total time points

    for i in range(0, len(arr_m) - 2):
        video.write(images[i])

    cv2.destroyAllWindows()
    video.release()

    concatenated2 = np.concatenate([x for x in images[0 : len(arr_m) - 1]], axis=0)

    cen_whi = [[]] * l_asg
    for i in range(0, l_asg):
        if Mod[i] == "move":
            cen_whi[i] = asg_arr[i][-1:]
        elif Mod[i] == "divide":
            cen_whi[i] = asg_arr[i][-2:]
        elif Mod[i] == "exit" or Mod[i] == "last":
            cen_whi[i] = 0
    df_list_concat["centers_which"] = cen_whi

    cen_val = [[]] * l_asg  ###len(df_list_concat)

    for i in range(0, l_asg):  ###len(df_list_concat)):
        cen_val[i] = []
        for k in range(0, l_asg):  ### len(df_list_concat)):
            if type(cen_whi[i]) != int and cen_whi[i] != 0:
                for j in range(0, len(cen_whi[i])):
                    if (
                        df_list_concat["t"].values[k]
                        == df_list_concat["t"].values[i] + 1
                    ) and (df_list_concat["index"].values[k] == int(cen_whi[i][j])):
                        cen_val[i] += [str((df_list_concat["center"].values[k]))]
    df_list_concat["center_value"] = cen_val

    df_list_concat.to_csv("df_label_list_new_yay.csv")
    for i in range(0, len(arr_m) - 1):
        df_see = df_list_concat[df_list_concat["t"] == i]
        df_see = df_see.reset_index(drop=True)
        for k in range(len(df_see)):
            for j in range(len(df_see["center_value"][k])):
                x1 = int(df_see["center"].values[k])
                x2 = int(float(df_see["center_value"].values[k][j]))
                concatenated2 = cv2.line(
                    concatenated2,
                    (x1, int(38.5) + height * (i)),
                    (x2, int(38.5 + (height * (i + 1)))),
                    (255, 0, 0),
                    5,
                )

    concatenated2 = Image.fromarray(concatenated2)
    concatenated2.save("mask_img.jpg")


img_out = main2(69)  # 18
