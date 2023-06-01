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
            arr1 = get_coordinates(imgs, i + 1)
            arr_array[i] = arr1
        return arr_array

    # @profile
    def get_coord_min_max(threshold_number, mask_number):
        mask_number = mask_number - 1
        array1 = arr[mn][threshold_number][mask_number]
        max_y = max(array1, key=lambda item: item[0])
        min_y = min(array1, key=lambda item: item[0])
        max_x = max(array1, key=lambda item: item[1])
        min_x = min(array1, key=lambda item: item[1])

        arr_lengths = [max_y[0], min_y[0], max_x[1], min_x[1]]
        arr_hist = img_ini[
            (arr_lengths[0] + arr_lengths[1]) // 2, arr_lengths[3] : arr_lengths[2] + 1
        ]
        plt.plot(arr_hist)
        return arr_lengths, arr_hist, plt.plot(arr_hist)

    # @profile
    def get_center(threshold_number, mask_number):
        mask_number = mask_number - 1
        array1 = arr[mn][threshold_number][mask_number]
        max_y = max(array1, key=lambda item: item[0])
        min_y = min(array1, key=lambda item: item[0])
        max_x = max(array1, key=lambda item: item[1])
        min_x = min(array1, key=lambda item: item[1])

        arr_center = [(max_y[0] + min_y[0]) / 2, (max_x[1] + min_x[1]) / 2]
        arr_center_x = arr_center[1]
        return arr_center_x

    # @profile
    def get_lenx(threshold_number, mask_number):
        mask_number = mask_number - 1
        array1 = arr[mn][threshold_number][mask_number]
        max_x = max(array1, key=lambda item: item[1])
        min_x = min(array1, key=lambda item: item[1])
        len_x = max_x[1] - min_x[1]
        return len_x

    # @profile
    def get_img(threshold_number, mask_number):
        img3 = np.where(imgs[threshold_number] == mask_number, th[threshold_number], 7)
        return Image.fromarray(img3).show()

    # @profile
    def get_cost(threshold_number, mask_number):
        mask_number = mask_number
        arr_get, arr_hist2, img2 = get_coord_min_max(threshold_number, mask_number)
        number_of_peaks = scipy.signal.find_peaks(arr_hist2[:, 0])
        length = arr_get[2] - arr_get[3]
        if length < 20:
            length = length + 100000
        else:
            length = length / 100
        cost = [len(number_of_peaks[0]), length]
        costs = (cost[0] * 50 * 2) + (cost[1] * 2)
        return cost, costs

    def get_center2(threshold_number, mask_number):
        mask_number = mask_number - 1
        array1 = arr[mn][threshold_number][mask_number]

        # max_x = max(array1, key = lambda item : item[1])
        # min_x = min(array1, key = lambda item : item[1])
        range_x = [array1[0][1], array1[-1][1]]
        # range_x = [min_x[1], max_x[1]]
        # print(range_x, range_x2)
        return range_x

    def tree3a(node_a, node_b):
        tree3 = []

        for j in range(n_mask_img[mn][node_a]):
            ag = get_center2(node_a, j)
            for i in range(n_mask_img[mn][node_b]):
                if get_center(node_b, i) > ag[0] and get_center(node_b, i) <= ag[1]:
                    tree3 += [[i + 1, j + 1]]

        # print('>',tree3)
        return tree3

    arr_m = [0] * mainhuh
    df_list = [pd.DataFrame()] * len(arr_m)
    s = [[-1] * 1000] * len(arr_m)
    s2 = [[]] * len(arr_m)
    a_exit = [[-1] * 1000] * len(arr_m)
    a_div2 = [[-1] * 1000] * len(arr_m)
    a_div3 = [[-1] * 1000] * len(arr_m)
    imgs_mat = [[[]]] * len(arr_m)
    th_mat = [[[]]] * len(arr_m)
    n_constr = [0] * len(arr_m)
    # n_constr2=[0]*len(arr_m)
    n_constr_a = [0] * len(arr_m)
    n_constr_b = [0] * len(arr_m)
    n_constr_c = [0] * len(arr_m)
    n_mask_img = [[]] * len(arr_m)
    model_tree = Model()
    arr = [[[]]] * len(arr_m)
    len_df = [0] * len(arr_m)
    for mn in range(len(arr_m)):
        img_ini = cv2.imread("%d.jpg" % (mn + 13))
        th = [[]] * 4
        imgs = [[]] * 4
        reta, th[0] = cv2.threshold(img_ini, 187, 255, cv2.THRESH_BINARY)  # 160
        ret, th[1] = cv2.threshold(img_ini, 190, 255, cv2.THRESH_BINARY)  # 180
        reta, th[2] = cv2.threshold(img_ini, 205, 255, cv2.THRESH_BINARY)  # 200
        reta, th[3] = cv2.threshold(img_ini, 215, 255, cv2.THRESH_BINARY)  # 220

        n_mask_img[mn] = [0] * (len(imgs))
        arr[mn] = [[]] * (len(imgs))
        costs_get_arr2 = []
        center_arr = []
        lenx_arr = []
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
                costs_get_arr2 += [[0, i, j]]
                center_arr += [get_center(i, j)]
                lenx_arr += [get_lenx(i, j)]

        df["threshold_number"] = [x[1] for x in costs_get_arr2]
        df["mask"] = [x[2] for x in costs_get_arr2]
        df["s_i"] = [0 for x in costs_get_arr2]
        df_list[mn] = df
        df_list[mn]["center"] = center_arr
        df_list[mn]["length_x"] = lenx_arr
        df_list[mn]["t"] = mn

        # tree2_1 = tree2(0,3)
        # print(tree2_1)
        ############################## ASSIGNMENT VARIABLES ####################################################################

        matrix2 = []
        a_div2[mn - 1] = [[0] * len(df_list[mn])] * len(df_list[mn - 1])
        cost_move = [[0] * len(df_list[mn])] * len(df_list[mn - 1])
        for j in range(len(df_list[mn - 1])):
            row = [0] * len(df_list[mn])
            for i in range(len(df_list[mn])):
                cost_move[j][i] = 60 * abs(
                    lenx_arr[i] - df_list[mn - 1]["length_x"][j]
                ) + 250 * abs(center_arr[i] - df_list[mn - 1]["center"][j])
                xm = model_tree.addVar(
                    0,
                    1,
                    obj=cost_move[j][i],
                    vtype=GRB.BINARY,
                    name="x_{0}_{1}_{2}".format(mn - 1, j, i),
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
                20
                * (
                    df_list[mn]["length_x"].values[np.newaxis, :, np.newaxis]
                    + df_list[mn]["length_x"].values[np.newaxis, np.newaxis, :]
                    - df_list[mn - 1]["length_x"].values[:, np.newaxis, np.newaxis]
                )
                + 300
                * abs(
                    df_list[mn]["center"].values[np.newaxis, :, np.newaxis]
                    - (df_list[mn - 1]["center"].values[:, np.newaxis, np.newaxis])
                )
                + 100
                * np.abs(
                    df_list[mn]["center"].values[np.newaxis, np.newaxis, :]
                    - df_list[mn]["center"].values[np.newaxis, :, np.newaxis]
                )
                + 50
                * abs(
                    df_list[mn]["center"].values[np.newaxis, np.newaxis, :]
                    + (
                        df_list[mn]["center"].values[np.newaxis, :, np.newaxis]
                        - 2
                        * (df_list[mn - 1]["center"].values[:, np.newaxis, np.newaxis])
                    )
                )
            )

            for k in range(len(df_list[mn - 1])):
                row = [0] * len(df_list[mn])
                for i in range(len(df_list[mn])):
                    ro2 = [0] * len(df_list[mn])
                    for j in range(len(df_list[mn])):
                        xm = model_tree.addVar(
                            0,
                            1,
                            obj=cost_divide[k, i, j],
                            vtype=GRB.BINARY,
                            name="x_{0}_{1}_{2}_{3}".format(mn - 1, k, i, j),
                        )
                        ro2[j] = xm
                    row[i] = ro2
                matrix3 = row
                a_div3[mn - 1][k] = matrix3

            # cost_exit = []
            a_exit[mn - 1] = [0] * len(df_list[mn - 1])
            for k in range(0, len(df_list[mn - 1])):
                a1 = model_tree.addVar(
                    0, 1, obj=20000, vtype=GRB.BINARY, name="a_exit%d,%d" % (mn - 1, k)
                )
                a_exit[mn - 1][k] = a1

        if mn > 1:
            aac2a = [LinExpr(0)] * len(df_list[mn - 1])
            aac3a = [LinExpr(0)] * len(df_list[mn - 1])
            aac4a = [LinExpr(0)] * len(df_list[mn - 1])
            s[mn - 1] = [-1] * 1000
            for i in range(0, len(df_list[mn - 1])):
                aac2a[i] = sum(
                    a_div2[mn - 2][k][i] for k in range(0, len(df_list[mn - 2]))
                )
                aac3a[i] = sum(
                    a_div3[mn - 2][k][i][j] + a_div3[mn - 2][k][j][i]
                    for j in range(0, len(df_list[mn - 1]))
                    for k in range(0, len(df_list[mn - 2]))
                )
                aac4a[i] = (
                    sum(
                        np.sum(a_div3[mn - 1][i][j]) for j in range(0, len(df_list[mn]))
                    )
                    + np.sum(a_div2[mn - 1][i])
                    + a_exit[mn - 1][i]
                )
                s1 = aac4a[i]
                con2a = model_tree.addConstr(aac2a[i] + aac3a[i] - aac4a[i] == 0)
                s[mn - 1][i] = s1

        if mn == 1:
            s[mn - 1] = [-1] * 1000
            for i in range(0, len(df_list[mn - 1])):
                s1 = (
                    sum(
                        np.sum(a_div3[mn - 1][i][j]) for j in range(0, len(df_list[mn]))
                    )
                    + np.sum(a_div2[mn - 1][i])
                    + a_exit[mn - 1][i]
                )
                s[mn - 1][i] = s1

            ####################################################################################################

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
            mn = mn - 1  # to manage tree

            tree2_1 = tree3a(0, 3)
            # print( type(tree2_1))
            tree2_2 = tree3a(1, 3)
            tree2_3 = tree3a(2, 3)
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

        ####################################################################################################

        imgs_mat[mn] = imgs
        th_mat[mn] = th
        len_df[mn] = len(df_list[mn])

        df_list[mn - 1]["moves"] = a_div2[mn - 1]
        df_list[mn - 1]["divides"] = a_div3[mn - 1]
        df_list[mn - 1]["exits"] = a_exit[mn - 1]

    model_tree.update()

    # Set the TuneResults parameter to 1
    # model_tree.Params.tuneResults = 1

    # Tune the model
    # model_tree.tune()

    # if model_tree.tuneResultCount > 0:

    # Load the best tuned parameters into the model
    #   model_tree.getTuneResult(0)

    # Write tuned parameters to a file
    #  model_tree.write('tune.prm')

    # Solve the model using the tuned parameters
    # model_tree.optimize()

    # model_tree.Params.Threads = 1
    model_tree.optimize()
    # model_tree.write("file.lp")

    asg_arr2 = []
    for v in model_tree.getVars():
        if v.x == 1:
            asg_arr2 += [v.varName]

    asg_arr2a = []
    for i in range(len(asg_arr2)):
        if int(asg_arr2[i][2]) == (len(arr_m) - 2):
            asg_arr2a += [asg_arr2[i]]

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

    index2 = df_list_concat["index"]

    t1 = df_list_concat["t"]
    for j in range(len(arr_m)):
        sm = sum(len_df[m] for m in range(0, j))
        for i in range(len(t1)):
            if t1[i] == j:
                index2.values[i] -= sm
    df_list_concat["index"] = index2

    asg_arr = []

    asg_arr = asg_arr2

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
        if len(asg_arr[i]) == 5:
            Mod[i] = "divide"
        if len(asg_arr[i]) == 1:
            Mod[i] = "last"
        if len(asg_arr[i]) == 2:
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

    concatenated2 = np.concatenate([x for x in images[0 : len(arr_m)]], axis=0)

    cen_whi = [[]] * l_asg
    for i in range(0, l_asg):
        if Mod[i] == "move":
            cen_whi[i] = asg_arr[i][-1:]
        if Mod[i] == "divide":
            cen_whi[i] = asg_arr[i][-2:]
        if Mod[i] == "exit" or Mod[i] == "last":
            cen_whi[i] = 0
    df_list_concat["centers_which"] = cen_whi

    cen_val = [[]] * len(df_list_concat)

    for i in range(0, len(df_list_concat)):
        cen_val[i] = []
        for k in range(0, len(df_list_concat)):
            if type(cen_whi[i]) != int and cen_whi[i] != 0:  ##############
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
            for j in range(len(df_see["center_value"].values[k])):
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
    concatenated2.save("time_img2.jpg")


imgblah = main2(5)
