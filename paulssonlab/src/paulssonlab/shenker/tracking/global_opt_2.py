import cv2
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import skimage
from PIL import Image  # , ImageDraw, ImageFilter
from skimage import morphology


def Make_d_nodes(array):
    node_arr = []
    for node in array:
        node = str(node)
        node1 = node
        node = [node + "a"]
        node1 = [node1 + "b"]
        node_arr += node + node1

    return node_arr


arr_m = [0] * 3
df_list = [pd.DataFrame()] * len(arr_m)
df_label_list = [pd.DataFrame()] * len(arr_m)

from gurobipy import *

model_tree = Model()

s = [[-1] * 1000] * len(arr_m)
a_div = [[[[-1] * 1000] * 1000] * 1000] * len(arr_m)
a_exit = [[-1] * 1000] * len(arr_m)
a_div2 = [[-1] * 1000] * len(arr_m)
a_div3 = [[-1] * 1000] * len(arr_m)
a_move = [[[-1] * 1000] * 1000] * len(arr_m)
imgs_mat = [[[]]] * len(arr_m)
th_mat = [[[]]] * len(arr_m)

labels_arr = [[]] * len(arr_m)
ass_arr = [[[]]] * len(arr_m)
center = [[]] * len(arr_m)
len_x_m = [[]] * len(arr_m)
for mn in range(len(arr_m)):
    img_ini = cv2.imread("%d.jpg" % (mn + 3))
    th = [[]] * 4
    imgs = [[]] * 4
    reta, th[0] = cv2.threshold(img_ini, 160, 255, cv2.THRESH_BINARY)
    ret, th[1] = cv2.threshold(img_ini, 180, 255, cv2.THRESH_BINARY)
    reta, th[2] = cv2.threshold(img_ini, 200, 255, cv2.THRESH_BINARY)
    reta, th[3] = cv2.threshold(img_ini, 220, 255, cv2.THRESH_BINARY)

    for i in range(len(th)):
        img = skimage.morphology.label(
            th[i], neighbors=None, background=None, return_num=False, connectivity=None
        )
        imgs[i] = img

    n_mask_img = [0] * (len(imgs))
    for i in range(len(imgs)):
        n_mask_img[i] = np.amax(imgs[i])

    def get_coordinates(img1, mask_num):
        output2 = np.where(
            img1 == mask_num
        )  # [0] gives row indices and [1] gives column indices
        mask_1 = [mask_num for j in range(len(output2[0]))]
        listOfCoordinates = list(zip(output2[0], output2[1], mask_1))
        return listOfCoordinates

    def get_array(n_mask_img, imgs):
        arr_array = [[]] * n_mask_img
        for i in range(n_mask_img):
            arr1 = get_coordinates(imgs, i + 1)
            arr_array[i] = arr1
        return arr_array

    arr = [[]] * (len(imgs))

    for i in range(0, len(imgs)):
        arr2 = get_array(n_mask_img[i], imgs[i])
        arr[i] = arr2

    # arr is the 2 X 2 matrix of image number and mask numer : arr[img_number][mask_number] last is arr[2][6] => (y,x,mask)

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

    def tree(node1, node2):
        tree1 = get_tree(n_mask_img[node1], arr[node1], n_mask_img[node2], arr[node2])

        return tree1

    def get_coord_min_max(threshold_number, mask_number):
        mask_number = mask_number - 1
        array1 = arr[threshold_number][mask_number]
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

    def get_center(threshold_number, mask_number):
        mask_number = mask_number - 1
        array1 = arr[threshold_number][mask_number]
        max_y = max(array1, key=lambda item: item[0])
        min_y = min(array1, key=lambda item: item[0])
        max_x = max(array1, key=lambda item: item[1])
        min_x = min(array1, key=lambda item: item[1])

        arr_center = [(max_y[0] + min_y[0]) / 2, (max_x[1] + min_x[1]) / 2]
        arr_center_x = arr_center[1]
        return arr_center_x

    def get_lenx(threshold_number, mask_number):
        mask_number = mask_number - 1
        array1 = arr[threshold_number][mask_number]
        max_x = max(array1, key=lambda item: item[1])
        min_x = min(array1, key=lambda item: item[1])
        len_x = max_x[1] - min_x[1]
        return len_x

    def get_img(threshold_number, mask_number):
        img3 = np.where(
            imgs[threshold_number] == mask_number, th[threshold_number], 7
        )  # this prints only the areas which are labeled as img == label #img1a==5 img1==2
        return Image.fromarray(img3).show()

    def get_cost(threshold_number, mask_number):
        mask_number = mask_number
        arr_get, arr_hist2, img2 = get_coord_min_max(threshold_number, mask_number)
        number_of_peaks = scipy.signal.find_peaks(arr_hist2[:, 0])
        length = arr_get[2] - arr_get[3]
        cost = [len(number_of_peaks[0]), length]
        if length < 10 or length == 0:
            cost[1] == 1000
        costs = cost[0] * 5 + abs((cost[1] - 100))
        return cost, costs

    costs_get_arr = []
    costs_get_arr2 = []
    center_arr = []
    lenx_arr = []
    peaks = []
    df = pd.DataFrame()

    for i in range(0, len(imgs)):
        for j in range(0 + 1, (n_mask_img[i]) + 1):
            cost_get, costs_get = get_cost(i, j)
            costs_get_arr += [costs_get]
            costs_get_arr2 += [[costs_get, i, j]]
            center_arr += [get_center(i, j)]
            lenx_arr += [get_lenx(i, j)]
            peaks += [cost_get[0]]

    df["cost"] = [x[0] for x in costs_get_arr2]
    df["threshold_number"] = [x[1] for x in costs_get_arr2]
    df["mask"] = [x[2] for x in costs_get_arr2]
    df["s_i"] = [0 for x in costs_get_arr2]
    df["a_exit"] = [0 for x in costs_get_arr2]
    df["a_move"] = [0 for x in costs_get_arr2]
    df["a_div"] = [0 for x in costs_get_arr2]

    df_list[mn] = df

    c = costs_get_arr

    df_list[mn]["center"] = center_arr
    df_list[mn]["length_x"] = lenx_arr
    df_list[mn]["peaks"] = peaks

    s[mn] = [-1] * 1000
    for i in range(0, len(df_list[mn])):
        s1 = model_tree.addVar(
            0, 1, obj=c[i], vtype=GRB.INTEGER, name="node_{0}_{1}".format(mn, i)
        )  # name="node%d,%d"%(mn,i))
        s[mn][i] = s1
        model_tree.update()

    n_constr = [0] * len(arr_m)
    n_constr_a = [0] * len(arr_m)
    n_constr_b = [0] * len(arr_m)
    n_constr_c = [0] * len(arr_m)

    n_constr[mn] = n_mask_img[len(th) - 1]
    n_constr_c[mn] = (
        n_mask_img[len(th) - 2] + n_mask_img[len(th) - 3] + n_mask_img[len(th) - 4] - 1
    )
    n_constr_a[mn] = n_mask_img[len(th) - 3] + n_mask_img[len(th) - 4] - 1
    n_constr_b[mn] = n_mask_img[len(th) - 4] - 1

    for i in range(1, n_constr[mn] + 1):
        con1 = model_tree.addConstr(
            s[mn][
                (tree(0, 3)[int(np.where([x[0] == i for x in tree(0, 3)])[0])][1]) - 1
            ]
            + s[mn][
                (n_constr_b[mn])
                + (tree(1, 3)[int(np.where([x[0] == i for x in tree(1, 3)])[0])][1])
            ]
            + s[mn][
                n_constr_a[mn]
                + (tree(2, 3)[int(np.where([x[0] == i for x in tree(2, 3)])[0])][1])
            ]
            + s[mn][n_constr_c[mn] + i]
            == 1
        )

        model_tree.update()
    # model_tree.update()

    matrix2 = []

    for i in range(len(df_list[mn])):
        row = [0] * len(df_list[mn - 1])
        for j in range(len(df_list[mn - 1])):
            xm = model_tree.addVar(
                0, 1, obj=-20, vtype=GRB.INTEGER, name="x_{0}_{1}_{2}".format(mn, i, j)
            )
            model_tree.update()
            row[j] = xm
            # row.append(xm)
        matrix2 = row  # .append(row)

        a_div2[mn][i] = matrix2
        # print(mn,i,a_div2[mn][i],'\n', len(df_list[mn-1]))

    matrix3 = []

    for i in range(len(df_list[mn])):
        row = [0] * len(df_list[mn])
        for j in range(len(df_list[mn])):
            ro2 = [0] * len(df_list[mn - 1])
            for k in range(len(df_list[mn - 1])):
                xm = model_tree.addVar(
                    0,
                    1,
                    obj=-20,
                    vtype=GRB.INTEGER,
                    name="x_{0}_{1}_{2}_{3}".format(mn, i, j, k),
                )
                model_tree.update()
                ro2[k] = xm
                # ro2.append(xm)
            row[j] = ro2
            # row.append(ro2)
        matrix3 = row  # .append(row)

        a_div3[mn][i] = matrix3
        # print(mn,i,a_div3[mn][i],'\n')#, len(df_list[mn-1]))

    for i in range(0, len(df_list[mn])):
        a1 = model_tree.addVar(
            0, 1, obj=-20, vtype=GRB.INTEGER, name="a_exit%d,%d" % (mn, i)
        )
        model_tree.update()
        a_exit[mn][i] = a1

    ######################
    if mn > 0:
        aaaab = [LinExpr()] * len(df_list[mn])

        for i in range(0, len(df_list[mn])):
            for j in range(0, len(df_list[mn - 1])):
                aaaab[i] += np.sum(a_div3[mn][i][j])

            # aaaa=    np.matrix(a_div2[mn][i]).sum()
            # con2 = model_tree.addConstr( aaaa[i]+ sum(a_move[mn][i][j] for j in range(0,len(df_list[mn-1]))) + a_exit[mn][i] - s[mn][i]==0 )
            con2 = model_tree.addConstr(
                aaaab[i] + np.sum(a_div2[mn][i]) + a_exit[mn][i] - s[mn][i] == 0
            )
            model_tree.update()
            # con3 = model_tree.addConstr(aaaab[i] <=1 )

            # con4 = model_tree.addConstr( np.sum(a_div2[mn][i]) <= 1 )
            # con4 = model_tree.addConstr( a_exit[mn][0]==0 )

            # model_tree.update()

        # if mn>0:
        for i in range(0, len(df_list[mn])):
            for j in range(0, len(df_list[mn - 1])):
                # con = model_tree.addConstr( -2*np.sum(a_div2[mn][i])+ s[mn-1][j] + s[mn][i]>=0 )
                con = model_tree.addConstr(
                    -2 * a_div2[mn][i][j] + s[mn - 1][j] + s[mn][i] >= 0
                )

                model_tree.update()
                print(
                    mn, i, j, s[mn][i], s[mn - 1][j], -2 * np.sum(a_div2[mn][i]), "\n"
                )

        # aaaab2=[LinExpr()]*len(df_list[mn])
        for i in range(0, len(df_list[mn])):
            # aaaa=    np.sum(a_div2[mn][i])
            for j in range(0, len(df_list[mn])):
                # aaaab2[i] += np.sum(a_div3[mn][i][j])
                for k in range(0, len(df_list[mn - 1])):
                    if j != k:
                        # con2 = model_tree.addConstr(-3*aaaab2[i] + s[mn][i] +s[mn][j] +s[mn-1][k] >=0 )
                        con2 = model_tree.addConstr(
                            -3 * a_div3[mn][i][j][k]
                            + s[mn][i]
                            + s[mn][j]
                            + s[mn - 1][k]
                            >= 0
                        )

                        model_tree.update()

    ######################################
    """for i in range (0,len(df_list[mn])):

            con2 = model_tree.addConstr(a_exit[mn][i]==0 )
            model_tree.update() """

    """asum=0
    bsum=0
    #for i in range (0,len(df_list[mn])):

    for j in range (0,len(df_list[mn-1])):
            bsum += s[mn-1][j]


    for i in range (0,len(df_list[mn])):
        asum += s[mn][i]
        con2 = model_tree.addConstr(a_exit[mn][i] +asum - bsum<=0 )
        model_tree.update() """

    imgs_mat[mn] = imgs
    th_mat[mn] = th
    df_list[mn]["si"] = s[mn][0 : len(df_list[mn])]

model_tree.optimize()
constrs_mat = model_tree.getConstrs()
model_tree.write("file.lp")

a_exit_arr = []
for v in model_tree.getVars():
    if v.varName[0] == "x" and v.x == 1:  # & v.x==1):
        # a_exit_arr +=[v.x]
        print("%s %g" % (v.varName, v.x))


for v in model_tree.getVars():
    if v.varName[0] == "a" and v.x == 1:  # & v.x==1):
        # a_exit_arr +=[v.x]
        print("%s %g" % (v.varName, v.x))

s_arr = []  # [0] * (len(df_list_concat)*2)

for v in model_tree.getVars():
    # print('%s %g' % (v.varName, v.x))
    if v.varName[0] == "n":  # .split(',')[-1])==mn):
        s_arr += [v.x]
        if v.x == 1:
            print("%s %g" % (v.varName, v.x))

df_list_concat = pd.concat([df_list[mn] for mn in range(0, len(arr_m))])
df_list_concat["s_i"] = s_arr
df_list_concat = df_list_concat.reset_index(drop=True)

x_df = 0
for jk in range(0, len(arr_m)):
    # x_df = 0 #len(df_list[jk])
    for i in range(0, len(df_list[jk])):
        # dfm =  df_list[jk]
        df_list[jk]["s_i"][i] = df_list_concat["s_i"][i + x_df]
        # print(df_list_concat['s_i'][i+x_df] ,'m')
    x_df += len(df_list[jk])
    # print(x_df)


"""a_exit_arr = []
for v in model_tree.getVars():
        if (v.varName[0]=='a' and v.x==1):# & v.x==1):
            a_exit_arr +=[v.x]
            print('%s %g' % (v.varName, v.x))"""
