from gurobipy import *


# @profile
def main2(mainhuh, matrix_cost_move, matrix_cost_divide, imgd):

    import skimage

    # from PIL import Image
    import cv2
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy
    import pandas as pd
    from skimage import morphology

    from PIL import Image, ImageDraw, ImageFilter

    import skimage
    import random
    import math

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

            arr1 = get_coordinates(imgs, i + 1)

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

        if len(number_of_peaks[0]) == len(number_of_peaks2[0]):
            pp = (len(number_of_peaks[0]) + len(number_of_peaks2[0])) / 2
        else:
            pp = len(number_of_peaks[0]) + 1

        len_x = max_x[1] - min_x[1]
        right_edge = max_x[1]
        left_edge = min_x[1]
        return int(arr_center_x), pp, len_x, right_edge, left_edge

    def tree2(x):

        tree2_1 = get_tree(
            n_mask_img[mn][0], arr[mn][0], n_mask_img[mn][12], arr[mn][12]
        )
        tree2_2 = get_tree(
            n_mask_img[mn][1], arr[mn][1], n_mask_img[mn][12], arr[mn][12]
        )
        tree2_3 = get_tree(
            n_mask_img[mn][2], arr[mn][2], n_mask_img[mn][12], arr[mn][12]
        )
        tree2_4 = get_tree(
            n_mask_img[mn][3], arr[mn][3], n_mask_img[mn][12], arr[mn][12]
        )
        tree2_5 = get_tree(
            n_mask_img[mn][4], arr[mn][4], n_mask_img[mn][12], arr[mn][12]
        )
        tree2_6 = get_tree(
            n_mask_img[mn][5], arr[mn][5], n_mask_img[mn][12], arr[mn][12]
        )
        tree2_7 = get_tree(
            n_mask_img[mn][6], arr[mn][6], n_mask_img[mn][12], arr[mn][12]
        )
        tree2_8 = get_tree(
            n_mask_img[mn][7], arr[mn][7], n_mask_img[mn][12], arr[mn][12]
        )
        tree2_9 = get_tree(
            n_mask_img[mn][8], arr[mn][8], n_mask_img[mn][12], arr[mn][12]
        )
        tree2_10 = get_tree(
            n_mask_img[mn][9], arr[mn][9], n_mask_img[mn][12], arr[mn][12]
        )
        tree2_11 = get_tree(
            n_mask_img[mn][10], arr[mn][10], n_mask_img[mn][12], arr[mn][12]
        )
        tree2_12 = get_tree(
            n_mask_img[mn][11], arr[mn][11], n_mask_img[mn][12], arr[mn][12]
        )

        return (
            tree2_1,
            tree2_2,
            tree2_3,
            tree2_4,
            tree2_5,
            tree2_6,
            tree2_7,
            tree2_8,
            tree2_9,
            tree2_10,
            tree2_11,
            tree2_12,
        )

    arr_m = [0] * mainhuh
    df_list = [pd.DataFrame()] * len(arr_m)
    s = [[7] * 10000] * len(arr_m)
    s2 = [[]] * len(arr_m)
    a_exit = [[7] * 10000] * len(arr_m)
    a_div2 = [[7] * 10000] * len(arr_m)
    a_div3 = [[7] * 10000] * len(arr_m)
    imgs_mat = [[[]]] * len(arr_m)
    th_mat = [[[]]] * len(arr_m)
    n_constr = [0] * len(arr_m)
    n_constr_a = [0] * len(arr_m)
    n_constr_b = [0] * len(arr_m)
    n_constr_c = [0] * len(arr_m)
    n_constr_d = [0] * len(arr_m)
    n_constr_a1 = [0] * len(arr_m)
    n_constr_b1 = [0] * len(arr_m)
    n_constr_c1 = [0] * len(arr_m)
    n_constr_d1 = [0] * len(arr_m)
    n_constr_a2 = [0] * len(arr_m)
    n_constr_b2 = [0] * len(arr_m)
    n_constr_c2 = [0] * len(arr_m)
    n_constr_d2 = [0] * len(arr_m)
    n_mask_img = [[]] * len(arr_m)
    model_tree = Model()
    arr = [[[]]] * len(arr_m)
    len_df = [0] * len(arr_m)
    for mn in range(len(arr_m)):

        img_ini = cv2.imread("%d.jpg" % (mn + 0))
        th = [[]] * 13
        imgs = [[]] * 13

        reta, th[0] = cv2.threshold(
            img_ini, 140, 255, cv2.THRESH_BINARY
        )  # 160#187#188#180#####180#175
        ret, th[1] = cv2.threshold(
            img_ini, 145, 255, cv2.THRESH_BINARY
        )  # 180#190#195#185###190#180
        reta, th[2] = cv2.threshold(
            img_ini, 150, 255, cv2.THRESH_BINARY
        )  # 200#205#205#######200#190
        reta, th[3] = cv2.threshold(
            img_ini, 155, 255, cv2.THRESH_BINARY
        )  # 220#215#212#######208#207
        reta, th[4] = cv2.threshold(
            img_ini, 160, 255, cv2.THRESH_BINARY
        )  # 220#215#212#######208#207##_##2
        reta, th[5] = cv2.threshold(
            img_ini, 165, 255, cv2.THRESH_BINARY
        )  # 160#187#188#180#####180#175
        ret, th[6] = cv2.threshold(
            img_ini, 170, 255, cv2.THRESH_BINARY
        )  # 180#190#195#185###190#180
        reta, th[7] = cv2.threshold(
            img_ini, 175, 255, cv2.THRESH_BINARY
        )  # 200#205#205#######200#190
        reta, th[8] = cv2.threshold(
            img_ini, 180, 255, cv2.THRESH_BINARY
        )  # 220#215#212#######208#207
        reta, th[9] = cv2.threshold(
            img_ini, 185, 255, cv2.THRESH_BINARY
        )  # 160#187#188#180#####180#175
        ret, th[10] = cv2.threshold(
            img_ini, 190, 255, cv2.THRESH_BINARY
        )  # 180#190#195#185###190#180
        reta, th[11] = cv2.threshold(
            img_ini, 195, 255, cv2.THRESH_BINARY
        )  # 200#205#205#######200#190
        reta, th[12] = cv2.threshold(
            img_ini, 200, 255, cv2.THRESH_BINARY
        )  # 220#215#212#######208#207

        n_mask_img[mn] = [0] * (len(imgs))
        arr[mn] = [[]] * (len(imgs))
        costs_get_arr2 = []
        center_arr = []
        lenx_arr = []
        peaks = []
        r_e = []
        l_e = []
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
                l_e += [cc[4]]

        df["threshold_number"] = [x[0] for x in costs_get_arr2]
        df["mask"] = [x[1] for x in costs_get_arr2]
        df_list[mn] = df
        df_list[mn]["center"] = center_arr
        df_list[mn]["length_x"] = lenx_arr
        df_list[mn]["t"] = mn
        df_list[mn]["peaks"] = peaks
        df_list[mn]["right_edge"] = r_e
        df_list[mn]["left_edge"] = l_e

        matrix2 = []
        a_div2[mn - 1] = [[0] * len(df_list[mn])] * len(df_list[mn - 1])
        cost_move = [[0] * len(df_list[mn])] * len(df_list[mn - 1])

        for j in range(len(df_list[mn - 1])):
            row = [0] * len(df_list[mn])
            for i in range(len(df_list[mn])):

                cost_move[j][i] = (
                    matrix_cost_move[0]
                    * (abs(center_arr[i] - df_list[mn - 1]["center"].values[j]))
                    + matrix_cost_move[1]
                    * abs(peaks[i] + df_list[mn - 1]["peaks"].values[j] - 2)
                    + +(matrix_cost_move[2] * df_list[mn - 1]["length_x"].values[j])
                )

                if lenx_arr[i] < 35 or df_list[mn - 1]["length_x"].values[j] < 35:
                    cost_move[j][i] = cost_move[j][i] + 1000000
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
                matrix_cost_divide[0]
                * abs(
                    df_list[mn]["center"].values[np.newaxis, np.newaxis, :]
                    - df_list[mn - 1]["center"].values[:, np.newaxis, np.newaxis]
                )
                + matrix_cost_divide[1]
                * abs(
                    df_list[mn]["center"].values[np.newaxis, np.newaxis, :]
                    - 2 * df_list[mn - 1]["center"].values[:, np.newaxis, np.newaxis]
                    + df_list[mn]["center"].values[np.newaxis, :, np.newaxis]
                )
                + +matrix_cost_divide[2]
                * (
                    df_list[mn]["peaks"].values[np.newaxis, :, np.newaxis]
                    + df_list[mn]["peaks"].values[np.newaxis, np.newaxis, :]
                    + df_list[mn - 1]["peaks"].values[:, np.newaxis, np.newaxis]
                    - 3
                )
                + +matrix_cost_divide[3]
                * abs(
                    df_list[mn]["center"].values[np.newaxis, np.newaxis, :]
                    - df_list[mn]["center"].values[np.newaxis, :, np.newaxis]
                )
                + +(
                    matrix_cost_divide[4]
                    * df_list[mn - 1]["length_x"].values[:, np.newaxis, np.newaxis]
                )
                + +matrix_cost_divide[5]
                * (
                    2 * df_list[mn]["length_x"].values[np.newaxis, :, np.newaxis]
                    - df_list[mn - 1]["length_x"].values[:, np.newaxis, np.newaxis]
                )
            )  # - 1.2*df_list[mn-1]['length_x'].values[:,np.newaxis,np.newaxis]#+\

            for k in range(len(df_list[mn - 1])):
                row = [0] * len(df_list[mn])
                c_i = [0] * len(df_list[mn])
                for i in range(len(df_list[mn])):
                    ro2 = [0] * len(df_list[mn])
                    c_i[i] = r_e[i]
                    if c_i[i] != max(center_arr):
                        index_i = [
                            i2
                            for i2 in range(len(df_list[mn]))
                            if center_arr[i2] > c_i[i]
                            and center_arr[i2] < (c_i[i] + (1.2 * lenx_arr[i]))
                        ]
                        for j in index_i:

                            if (
                                df_list[mn]["length_x"].values[i] < 35
                                or df_list[mn - 1]["length_x"].values[k] < 35
                                or df_list[mn]["length_x"].values[j] < 35
                            ):
                                cost_divide[k, i, j] = cost_divide[k, i, j] + 1000000
                            else:
                                cost_divide[k, i, j] = cost_divide[k, i, j]
                            xm = model_tree.addVar(
                                0,
                                1,
                                obj=cost_divide[k, i, j],
                                vtype=GRB.BINARY,
                                name="x_%d_%d_%d_%d" % (mn - 1, k, i, j),
                            )
                            ro2[j] = xm
                    row[i] = ro2
                matrix3 = row
                a_div3[mn - 1][k] = matrix3

            a_exit[mn - 1] = [0] * len(df_list[mn - 1])

            for k in range(0, len(df_list[mn - 1])):
                x_exit = 0
                if df_list[mn - 1]["length_x"].values[k] < 35:
                    x_exit = 10000000
                if (df_list[mn - 1]["center"].values[k] > 1410) and (
                    len(df_list[mn]) < len(df_list[mn - 1])
                ):
                    x_exit = x_exit  # -100000
                a1 = model_tree.addVar(
                    0,
                    1,
                    obj=(80.58 + x_exit),
                    vtype=GRB.BINARY,
                    name="aexit_%d_%d" % (mn - 1, k),
                )
                a_exit[mn - 1][k] = a1

        lmn_1 = len(df_list[mn - 1])
        lmn = len(df_list[mn])

        if mn > 1:
            s[mn - 1] = [7] * 10000
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
            s[mn - 1] = [7] * 10000
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

        ##_##3
        if mn > 0 and mn < len(arr_m):
            n_constr[mn - 1] = n_mask_img[mn - 1][len(th_mat[mn - 1]) - 1]

            n_constr_d[mn - 1] = (
                n_mask_img[mn - 1][len(th_mat[mn - 1]) - 2]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 3]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 4]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 5]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 6]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 7]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 8]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 9]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 10]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 11]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 12]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 13]
                - 1
            )
            n_constr_c[mn - 1] = (
                n_mask_img[mn - 1][len(th_mat[mn - 1]) - 3]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 4]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 5]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 6]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 7]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 8]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 9]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 10]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 11]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 12]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 13]
                - 1
            )
            n_constr_a[mn - 1] = (
                n_mask_img[mn - 1][len(th_mat[mn - 1]) - 4]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 5]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 6]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 7]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 8]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 9]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 10]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 11]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 12]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 13]
                - 1
            )
            n_constr_b[mn - 1] = (
                n_mask_img[mn - 1][len(th_mat[mn - 1]) - 5]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 6]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 7]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 8]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 9]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 10]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 11]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 12]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 13]
                - 1
            )
            n_constr_d1[mn - 1] = (
                n_mask_img[mn - 1][len(th_mat[mn - 1]) - 6]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 7]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 8]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 9]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 10]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 11]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 12]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 13]
                - 1
            )
            n_constr_c1[mn - 1] = (
                n_mask_img[mn - 1][len(th_mat[mn - 1]) - 7]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 8]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 9]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 10]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 11]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 12]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 13]
                - 1
            )
            n_constr_a1[mn - 1] = (
                n_mask_img[mn - 1][len(th_mat[mn - 1]) - 8]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 9]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 10]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 11]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 12]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 13]
                - 1
            )
            n_constr_b1[mn - 1] = (
                n_mask_img[mn - 1][len(th_mat[mn - 1]) - 9]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 10]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 11]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 12]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 13]
                - 1
            )
            n_constr_d2[mn - 1] = (
                n_mask_img[mn - 1][len(th_mat[mn - 1]) - 10]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 11]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 12]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 13]
                - 1
            )
            n_constr_c2[mn - 1] = (
                n_mask_img[mn - 1][len(th_mat[mn - 1]) - 11]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 12]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 13]
                - 1
            )
            n_constr_a2[mn - 1] = (
                n_mask_img[mn - 1][len(th_mat[mn - 1]) - 12]
                + n_mask_img[mn - 1][len(th_mat[mn - 1]) - 13]
                - 1
            )
            n_constr_b2[mn - 1] = n_mask_img[mn - 1][len(th_mat[mn - 1]) - 13] - 1

        if mn > 0 and mn < len(arr_m):
            mn = mn - 1

            tr = tree2(0)
            tree2_1 = tr[0]
            tree2_2 = tr[1]
            tree2_3 = tr[2]
            tree2_4 = tr[3]
            tree2_5 = tr[4]
            tree2_6 = tr[5]
            tree2_7 = tr[6]
            tree2_8 = tr[7]
            tree2_9 = tr[8]
            tree2_10 = tr[9]
            tree2_11 = tr[10]
            tree2_12 = tr[11]

            for i in range(1, n_constr[mn] + 1):
                con1 = model_tree.addConstr(
                    s[mn][
                        (tree2_1[int(np.where([x[0] == i for x in tree2_1])[0])][1]) - 1
                    ]
                    + s[mn][
                        (n_constr_b2[mn])
                        + (tree2_2[int(np.where([x[0] == i for x in tree2_2])[0])][1])
                    ]
                    + s[mn][
                        n_constr_a2[mn]
                        + (tree2_3[int(np.where([x[0] == i for x in tree2_3])[0])][1])
                    ]
                    + s[mn][
                        n_constr_c2[mn]
                        + (tree2_4[int(np.where([x[0] == i for x in tree2_4])[0])][1])
                    ]
                    + s[mn][
                        (n_constr_d2[mn])
                        + (tree2_5[int(np.where([x[0] == i for x in tree2_5])[0])][1])
                    ]
                    + s[mn][
                        n_constr_b1[mn]
                        + (tree2_6[int(np.where([x[0] == i for x in tree2_6])[0])][1])
                    ]
                    + s[mn][
                        n_constr_a1[mn]
                        + (tree2_7[int(np.where([x[0] == i for x in tree2_7])[0])][1])
                    ]
                    + s[mn][
                        (n_constr_c1[mn])
                        + (tree2_8[int(np.where([x[0] == i for x in tree2_8])[0])][1])
                    ]
                    + s[mn][
                        n_constr_d1[mn]
                        + (tree2_9[int(np.where([x[0] == i for x in tree2_9])[0])][1])
                    ]
                    + s[mn][
                        n_constr_b[mn]
                        + (tree2_10[int(np.where([x[0] == i for x in tree2_10])[0])][1])
                    ]
                    + s[mn][
                        (n_constr_a[mn])
                        + (tree2_11[int(np.where([x[0] == i for x in tree2_11])[0])][1])
                    ]
                    + s[mn][
                        n_constr_c[mn]
                        + (tree2_12[int(np.where([x[0] == i for x in tree2_12])[0])][1])
                    ]
                    + s[mn][n_constr_d[mn] + i]
                    == 1
                )

            mn = mn + 1

        ######EXIT  CONSTRAINT #####

        if mn > 0 and mn < len(arr_m):
            c_a = [0] * len(df_list[mn - 1])
            sum_see2 = [LinExpr()] * len(df_list[mn - 1])
            index_ia = [None] * len(df_list[mn - 1])
            for i in range(len(df_list[mn - 1])):
                c_a[i] = df_list[mn - 1]["right_edge"].values[i]
                if c_a[i] < max(df_list[mn - 1]["left_edge"]):
                    index_ia[i] = [
                        i2a
                        for i2a in range(len(df_list[mn - 1]))
                        if (df_list[mn - 1]["left_edge"].values[i2a] + 20 > (c_a[i]))
                    ]
                    sum_see2[i] = quicksum(
                        (s[mn - 1][kl] - a_exit[mn - 1][kl]) for kl in index_ia[i]
                    )
                    k = len(df_list[mn - 1])
                    con_exit = model_tree.addConstr(
                        ((k * a_exit[mn - 1][i]) + sum_see2[i]) <= k
                    )

        imgs_mat[mn] = imgs
        th_mat[mn] = th
        len_df[mn] = len(df_list[mn])

    model_tree.update()
    model_tree.Params.Threads = 1
    model_tree.optimize()
    model_tree.write("model.lp")

    asg_arr = []
    for v in model_tree.getVars():
        if v.x == 1:
            asg_arr += [v.varName]

    asg = pd.DataFrame()
    asg["0"] = asg_arr
    asg["int"] = [0] * len(asg)
    asg["t"] = [0] * len(asg)
    spl = [0] * len(asg)
    s_asg_arr = [[0] * 10] * len(asg)
    for ik in range(len(asg)):
        s_asg_arr[ik] = str(asg_arr[ik])
        spl = s_asg_arr[ik].split("_")
        asg["int"][ik] = spl[2]
        asg["t"][ik] = spl[1]

    asg_arrmm = [pd.DataFrame] * len(arr_m)

    for i in range(len(arr_m)):
        asg_arrmm[i] = asg[asg["t"] == i]
        asg_arrmm[i] = asg_arrmm[i].sort_values(by="int")

    asg_concat = pd.concat([asg_arrmm[mn] for mn in range(0, len(arr_m))])
    asg_arr = np.asarray(asg_concat["0"])
    asg_concat.to_csv("asg.csv")
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
    print(len(df_list_concat), len(asg_arr))

    if len(df_list_concat) > len(asg_arr):
        asg_arr += [0]

    if len(df_list_concat) < len(asg_arr):
        n = len(asg_arr) - len(df_list_concat)
        asg_arr = list(asg_arr)
        del asg_arr[-n:]

    l_asg = len(asg_arr)
    for i in range(l_asg):
        asg_arr2 = str(asg_arr[i])
        asg_arr[i] = asg_arr2.split("_")

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
        elif len(asg_arr[i]) == 3:
            Mod[i] = "exit"

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

    concatenated2 = np.concatenate([x for x in images[0 : len(arr_m) - 1]], axis=0)

    cen_whi = [[]] * l_asg
    for i in range(0, l_asg):
        if Mod[i] == "move":
            cen_whi[i] = asg_arr[i][-1:]
        elif Mod[i] == "divide":
            cen_whi[i] = asg_arr[i][-2:]
        elif Mod[i] == "exit" or Mod[i] == "last":
            cen_whi[i] = -1
    df_list_concat["centers_which"] = cen_whi

    cen_val = [[]] * l_asg

    for i in range(0, l_asg):
        cen_val[i] = []
        for k in range(0, l_asg):
            if type(cen_whi[i]) != int and cen_whi[i] != 0:
                for j in range(0, len(cen_whi[i])):
                    if (
                        df_list_concat["t"].values[k]
                        == df_list_concat["t"].values[i] + 1
                    ) and (df_list_concat["index"].values[k] == int(cen_whi[i][j])):
                        cen_val[i] += [str((df_list_concat["center"].values[k]))]
    df_list_concat["center_value"] = cen_val

    df_list_concat.to_csv("df_label_list_new_yay5.csv")
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
    concatenated2.save("mask_img_225%d.jpg" % imgd)


def scorefunc(imgd):
    import pandas as pd
    import numpy as np

    df = pd.read_csv("df_label_list_new_yay5.csv")
    df2 = pd.DataFrame()

    df2["index"] = df["index"]
    df2["MoD"] = df["MoD"]
    df2["center"] = df["center"]
    df2["centerwhich"] = df["center_value"]
    df2["t"] = df["t"]

    max_t = max(df2["t"]) + 1

    arr_center = [[]] * max_t
    arr_center_which = [[]] * max_t
    arr_Mod = [[]] * max_t
    for i in range(0, max_t):
        arr_center[i] = []
        arr_Mod[i] = []
        arr_center_which[i] = []
        for j in range(len(df2)):
            if df2["t"][j] == i:

                if arr_Mod[i] != "exit":
                    df2["centerwhich"][j] = df2["centerwhich"][j][2:]
                    df2["centerwhich"][j] = df2["centerwhich"][j][:-2]

                arr_center[i] += [df2["center"][j]]
                arr_center_which[i] += [df2["centerwhich"][j]]
                arr_Mod[i] += [df2["MoD"][j]]

    df3 = pd.DataFrame()

    df3["t"] = [0] * max_t
    df3["centers"] = [[]] * max_t
    df3["centerwhichs"] = [[]] * max_t
    df3["MoD"] = [[]] * max_t
    for i in range(max_t):
        df3["t"] = i

        Y = arr_center[i]
        X = arr_Mod[i]

        for j in range(len(Y)):
            yj = Y[j]
            xj = X[j]

            Y[j] = int(yj)

        X = np.array(X)
        Y = np.array(Y)

        inds = Y.argsort()
        sortedX = X
        sortedY = Y

        df3["centers"][i] = sortedY
        df3["MoD"][i] = sortedX
        df3["centerwhichs"][i] = arr_center_which[i]

    df2.to_csv("check.csv")
    df3 = df3.drop(columns=["t"])
    df3.to_csv("output_segmentation_tracking5.csv")

    import pandas as pd
    import numpy as np

    df_ground = pd.read_csv("ground_truth_b.csv")
    df_test = pd.read_csv("output_segmentation_tracking5.csv")

    df_test = pd.DataFrame(df_test)
    df_ground = pd.DataFrame(df_ground)
    df_ground = df_ground[:-1]

    for i in range(len(df_test)):
        df_test["centerwhichs"][i] = df_test["centerwhichs"][i][1:][:-1].split(", ")
        df_ground["to"][i] = df_ground["to"][i][1:][:-1].split(" ")
        df_test["centers"][i] = df_test["centers"][i][1:][:-1].split(" ")
        df_ground["center"][i] = df_ground["center"][i][1:][:-1].split(", ")

        new = [0] * len(df_test["centerwhichs"][i])
        for j in range(len(df_test["centerwhichs"][i])):
            aa = df_test["centerwhichs"][i][j][1:][:-1]
            new[j] = aa
        df_test["centerwhichs"][i] = new

    for i in range(len(df_test)):
        df_test["MoD"][i] = df_test["MoD"][i][1:][:-1].split(" ")
        df_ground["assign"][i] = df_ground["assign"][i][1:][:-1].split(", ")

    df_test["coms"] = [[]] * len(df_test)
    comb_arr = [[]] * len(df_test)

    for i in range(len(df_test) - 1):
        comb_arr[i] = [None] * (len(df_test["centers"][i]))
        countm = 0
        for j in range(len(df_test["centers"][i])):
            if df_test["MoD"][i][j][2] == "i":
                comb_arr[i][j] = np.array(
                    [
                        float(df_test["centers"][i][j]),
                        float(df_test["centerwhichs"][i][j + countm]),
                        float(df_test["centerwhichs"][i][j + countm + 1]),
                    ]
                )
                # print('see',i,j, comb_arr[i][j] , type(comb_arr[i][j]))
                countm += 1

            elif df_test["MoD"][i][j][2] == "o":

                a = df_test["centers"][i][j]
                b = df_test["centerwhichs"][i][j + countm]
                # if (len(a) >1 and len(b) >  1 ):
                a = float(a)
                b = float(b)
                # print(a,b,'-', a+b)
                comb_arr[i][j] = np.array([a, b])

                print("see", i, j, comb_arr[i][j], type(comb_arr[i][j]))

        comb_arr[i] = np.asarray(comb_arr[i])
        print("see2", i, comb_arr[i], type(comb_arr[i]))
        df_test["coms"][i] = comb_arr[i]

    df_ground["coms"] = [[]] * len(df_ground)
    comb_arr_g = [[]] * len(df_ground)

    for i in range(len(df_ground)):
        comb_arr_g[i] = [None] * (len(df_ground["center"][i]))
        countm2 = 0
        for j in range(len(df_ground["center"][i])):
            # print((df_ground['assign'][i][j]) )
            if df_ground["assign"][i][j][2] == "i":
                comb_arr_g[i][j] = np.array(
                    [
                        df_ground["center"][i][j],
                        df_ground["to"][i][j + countm2],
                        df_ground["to"][i][j + countm2 + 1],
                    ]
                )
                # print('see',i,j, comb_arr_g[i][j])
                countm2 += 1

            elif df_ground["assign"][i][j][2] == "o":

                an = df_ground["center"][i][j]
                bn = df_ground["to"][i][j + countm2]
                i  # f (len(an) >1 and len(bn) >  1 ):
                an = float(an)
                bn = float(bn)
                # print(a,b,'-', a+b)
                comb_arr_g[i][j] = np.array([an, bn])

                # print('see2',i,j, comb_arr_g[i][j])

        df_ground["coms"][i] = comb_arr_g[i]

    for i in range(len(df_test) - 1):

        Y = df_test["centers"][i]
        X = df_test["coms"][i]
        z = df_test["MoD"][i]
        X = np.asarray(X)
        Y = np.asarray(Y)
        z = np.asarray(z)
        inds = Y.argsort()
        # print(i,inds)
        sortedX = X[inds]
        sortedY = Y[inds]
        sortedz = z[inds]
        # print(sortedX,sortedY)
        df_test["centers"][i] = sortedY
        df_test["coms"][i] = sortedX
        df_test["MoD"][i] = sortedz

    score = [0] * (len(df_test))

    for i in range(len(df_test) - 1):
        score2 = 0
        score[i] = 0
        if len(df_test["centers"][i]) != len(df_ground["center"][i]):
            score2 = score2 + 10 * abs(
                len(df_test["centers"][i]) - len(df_ground["center"][i])
            )

        elif len(df_test["centers"][i]) == len(df_ground["center"][i]):
            for j in range(len(df_test["centers"][i])):
                # print(df_test['MoD'][i][j],'>', df_ground['assign'][i][j],'>')
                if (
                    df_test["MoD"][i][j][2] == "o"
                    and df_ground["assign"][i][j][2] == "o"
                ):
                    # print(df_test['coms'][i][j], df_test['coms'][i][j][1],'>>', df_ground['coms'][i][j][1])
                    # print(df_ground['coms'][i][j], df_test['coms'][i][j][0],'>>', df_ground['coms'][i][j][0])
                    # score2 = score2 + (.1) *  abs(float(df_test['centers'][i][j]) - float(df_ground['center'][i][j]))
                    score2 = score2 + abs(
                        float(df_test["coms"][i][j][1])
                        - float(df_ground["coms"][i][j][1])
                    )

                elif (
                    df_test["MoD"][i][j][2] == "i"
                    and df_ground["assign"][i][j][2] == "i"
                ):
                    # print(df_test['coms'][i][j], df_test['coms'][i][j][1],'>>', df_ground['coms'][i][j][1])
                    # print(df_ground['coms'][i][j], df_test['coms'][i][j][0],'>>', df_ground['coms'][i][j][0])
                    # score2 = score2 + (.1) *  abs(float(df_test['centers'][i][j]) - float(df_ground['center'][i][j]))
                    score2 = score2 + abs(
                        float(df_test["coms"][i][j][1])
                        + float(df_test["coms"][i][j][2])
                        - float(df_ground["coms"][i][j][1])
                        - float(df_ground["coms"][i][j][2])
                    )
                    # print(abs(float(df_test['coms'][i][j][1]) + float(df_test['coms'][i][j][2]) - float(df_ground['coms'][i][j][1]) - float(df_ground['coms'][i][j][2])),'>>>>>>>')
                else:
                    score2 = score2 + 100

        score[i] = score2
    score[len(df_test) - 1] = sum(score)

    score2 = score
    score = pd.DataFrame(score)
    score["score"] = score2
    score.to_csv("score64%d.csv" % imgd)

    # print(score)

    df_ground.to_csv("df_ground.csv")
    df_test.to_csv("df_test.csv")

    return score


# matrix_cost_move2 = [0, 0, 0]
# matrix_cost_divide2 = [0,0,0,0,0,0]


matrix_cost_move2 = [0, 0.050, 0.01]
matrix_cost_divide2 = [0, 0, 0.004, 0.0024, -10.1112, 300]


img_out = main2(7, matrix_cost_move2, matrix_cost_divide2, 0)
sc_i = scorefunc(0)
sc_i = sc_i[0][len(sc_i) - 1]

arr_score = [-1] * 100
arr_score[0] = sc_i
for i in range(1, 100):

    matrix_cost_move2[0] = matrix_cost_move2[0] + i * 1  # 8
    # matrix_cost_move2[1] = matrix_cost_move2[1] + i*0.0050
    # matrix_cost_move2[2] = matrix_cost_move2[2]+  i*0.001

    # matrix_cost_divide2[0] = matrix_cost_divide2[0]+ i*5
    # matrix_cost_divide2[1] = matrix_cost_divide2[1]+ i*0
    # matrix_cost_divide2[2] = matrix_cost_divide2[2]+ i*0.00004
    # matrix_cost_divide2[3] = matrix_cost_divide2[3]+ i*0.00024
    # matrix_cost_divide2[4] = matrix_cost_divide2[4]+ i*1
    # matrix_cost_divide2[5] = matrix_cost_divide2[5]+ i*30

    img_out = main2(7, matrix_cost_move2, matrix_cost_divide2, i)
    sc = scorefunc(i)
    sc = sc[0][len(sc) - 1]
    arr_score[i] = sc

    """if (arr_score[i] < arr_score[i-1]):
        print(arr_score[i], arr_score[i-1])
        continue
    else:
        print(arr_score[i], arr_score[i-1])
        break"""

    """if (arr_score[i] < arr_score[i-1]):
        print(arr_score[i], arr_score[i-1])
        continue
    else:
        print(arr_score[i], arr_score[i-1])
        matrix_cost_divide2[0] = matrix_cost_divide2[0]+ i*1#5
        
        if (arr_score[i] < arr_score[i-1]):
            print(arr_score[i], arr_score[i-1])
            continue
        else:
            print(arr_score[i], arr_score[i-1])

            break"""

    if arr_score[i] < arr_score[i - 1]:
        print(arr_score[i], arr_score[i - 1])
        continue
    else:
        print(arr_score[i], arr_score[i - 1])
        matrix_cost_divide2[0] = matrix_cost_divide2[0] + i * 1  # 5

        if arr_score[i] < arr_score[i - 1]:
            print(arr_score[i], arr_score[i - 1])
            continue
        else:
            matrix_cost_divide2[5] = matrix_cost_divide2[5] + i * 1
            print(arr_score[i], arr_score[i - 1])

            if arr_score[i] < arr_score[i - 1]:
                print(arr_score[i], arr_score[i - 1])
                continue
            else:
                print(arr_score[i], arr_score[i - 1])
                break

# matrix_cost_move2 = [84.226, 0.050, 0.01]
# matrix_cost_divide2 = [50.950,0,0.004,0.0024,-10.1112, 300.104]
