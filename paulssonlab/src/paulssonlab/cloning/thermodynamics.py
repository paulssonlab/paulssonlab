import Bio.SeqUtils.MeltingTemp as MeltingTemp


def tm(seq, method="q5", **kwargs):
    if method == "q5":
        return MeltingTemp.Tm_NN(
            seq, dnac1=500, nn_table=MeltingTemp.DNA_NN3, saltcorr=7, Na=150
        )
    elif method == "phusion":
        pass
    elif method == "":
        pass
    else:
        raise ValueError(f"unknown method: {method}")


def ta_from_tms(*tms, method="q5"):
    if method == "q5":
        # in the following,
        # l is min(len(primer1), len(primer2))
        # s is min(tm(primer1), tm(primer2))
        # NEB Tm calc JAVASCRIPT CODE:
        # l > 7 && (o = s + 1),
        # o > 72 && (o = 72);
        return min(min(tms) + 1, 72)
    elif method == "phusion":
        # NEB Tm calc JAVASCRIPT CODE:
        # (o = .93 * s + 7.5) > 72 && (o = 72);
        return min(0.93 * min(tms) + 7.5, 72)
    elif method == "":
        pass
    else:
        raise ValueError(f"unknown method: {method}")


def ta(*seqs, method="q5", **kwargs):
    tms = [tm(seq, method=method, **kwargs) for seq in seqs]
    return ta_from_tms(*tms, method=method)
