import skimage
from PIL import Image
import cv2
import numpy as np
import matplotlib.pyplot as plt
import scipy
import pandas as pd
from skimage import morphology
import math
import gurobipy
from gurobipy import *


# FUNCTIONS USED


def Make_d_nodes(array):

    node_arr = []
    for node in array:
        node = str(node)
        node1 = node
        node = [node + "a"]
        node1 = [node1 + "b"]
        node_arr += node + node1

    return node_arr


def get_coordinates(img1, mask_num):
    output2 = np.where(img1 == mask_num)
    mask_1 = [mask_num for j in range(len(output2[0]))]
    listOfCoordinates = list(zip(output2[0], output2[1], mask_1))
    return listOfCoordinates


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
    img3 = np.where(imgs[threshold_number] == mask_number, th[threshold_number], 7)
    return Image.fromarray(img3).show()


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


arr_m = [0] * 70
df_list = [pd.DataFrame()] * len(arr_m)
s = [[-1] * 1000] * len(arr_m)
a_exit = [[-1] * 1000] * len(arr_m)
a_div2 = [[-1] * 1000] * len(arr_m)
a_div3 = [[-1] * 1000] * len(arr_m)
imgs_mat = [[[]]] * len(arr_m)
th_mat = [[[]]] * len(arr_m)

model_tree = Model()

for mn in range(len(arr_m)):

    img_ini = cv2.imread("%d.jpg" % (mn))
    th = [[]] * 4
    imgs = [[]] * 4
    reta, th[0] = cv2.threshold(img_ini, 187, 255, cv2.THRESH_BINARY)  # 160
    ret, th[1] = cv2.threshold(img_ini, 190, 255, cv2.THRESH_BINARY)  # 180
    reta, th[2] = cv2.threshold(img_ini, 205, 255, cv2.THRESH_BINARY)  # 200
    reta, th[3] = cv2.threshold(img_ini, 215, 255, cv2.THRESH_BINARY)  # 220

    for i in range(len(th)):
        img = skimage.morphology.label(
            th[i], neighbors=None, background=None, return_num=False, connectivity=None
        )
        imgs[i] = img

    n_mask_img = [0] * (len(imgs))
    for i in range(len(imgs)):
        n_mask_img[i] = np.amax(imgs[i])

    arr = [[]] * (len(imgs))

    for i in range(0, len(imgs)):
        arr2 = get_array(n_mask_img[i], imgs[i])
        arr[i] = arr2

    # arr is the 2 X 2 matrix of image number and mask numer : arr[img_number][mask_number] last is arr[2][6] => (y,x,mask)

    costs_get_arr = []
    costs_get_arr2 = []
    center_arr = []
    lenx_arr = []
    # peaks=[]
    df = pd.DataFrame()

    for i in range(0, len(imgs)):
        for j in range(0 + 1, (n_mask_img[i]) + 1):
            cost_get, costs_get = get_cost(i, j)
            costs_get_arr += [costs_get]
            costs_get_arr2 += [[costs_get, i, j]]
            center_arr += [get_center(i, j)]
            lenx_arr += [get_lenx(i, j)]
            # peaks += [cost_get[0]]

    df["cost"] = [x[0] for x in costs_get_arr2]
    df["threshold_number"] = [x[1] for x in costs_get_arr2]
    df["mask"] = [x[2] for x in costs_get_arr2]
    df["s_i"] = [0 for x in costs_get_arr2]
    df["a_exit"] = [0 for x in costs_get_arr2]
    # df['t'] = [0 for x in costs_get_arr2]

    df_list[mn] = df

    c = costs_get_arr
    df_list[mn]["center"] = center_arr
    df_list[mn]["length_x"] = lenx_arr
    # df_list[mn]['peaks'] = peaks
    df_list[mn]["t"] = mn

    s[mn] = [-1] * 1000
    for i in range(0, len(df_list[mn])):
        s1 = model_tree.addVar(
            0, 1, obj=c[i], vtype=GRB.BINARY, name="node_{0}_{1}".format(mn, i)
        )
        s[mn][i] = s1

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

    ############################## ASSIGNMENT VARIABLES ####################################################################

    matrix2 = []
    a_div2[mn - 1] = [[0] * len(df_list[mn])] * len(df_list[mn - 1])
    cost_move = [[0] * len(df_list[mn])] * len(df_list[mn - 1])
    for j in range(len(df_list[mn - 1])):
        row = [0] * len(df_list[mn])
        for i in range(len(df_list[mn])):
            cost_move[j][i] = 60 * abs(
                df_list[mn]["length_x"][i] - df_list[mn - 1]["length_x"][j]
            ) + 250 * abs(df_list[mn]["center"][i] - df_list[mn - 1]["center"][j])
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
    a_div3[mn - 1] = [[[0] * len(df_list[mn])] * len(df_list[mn])] * len(
        df_list[mn - 1]
    )
    cost_divide = [[[0] * len(df_list[mn])] * len(df_list[mn])] * len(df_list[mn - 1])
    for k in range(len(df_list[mn - 1])):
        row = [0] * len(df_list[mn])
        for i in range(len(df_list[mn])):
            ro2 = [0] * len(df_list[mn])
            for j in range(len(df_list[mn])):
                cost_divide[k][i][j] = (
                    20
                    * (
                        df_list[mn]["length_x"][i]
                        + df_list[mn]["length_x"][j]
                        - df_list[mn - 1]["length_x"][k]
                    )
                    + 300
                    * abs(df_list[mn]["center"][i] - (df_list[mn - 1]["center"][k]))
                    + 100 * abs(df_list[mn]["center"][j] - (df_list[mn]["center"][i]))
                    + 50
                    * abs(
                        df_list[mn]["center"][j]
                        + (
                            df_list[mn]["center"][i]
                            - 2 * (df_list[mn - 1]["center"][k])
                        )
                    )
                )
                xm = model_tree.addVar(
                    0,
                    1,
                    obj=cost_divide[k][i][j],
                    vtype=GRB.BINARY,
                    name="x_{0}_{1}_{2}_{3}".format(mn - 1, k, i, j),
                )
                ro2[j] = xm
            row[i] = ro2
        matrix3 = row
        a_div3[mn - 1][k] = matrix3

    cost_exit = []
    for i in range(0, len(df_list[mn - 1])):
        a1 = model_tree.addVar(
            0, 1, obj=20000, vtype=GRB.BINARY, name="a_exit%d,%d" % (mn - 1, i)
        )
        a_exit[mn - 1][i] = a1
    ####################################################################################################

    for i in range(0, len(df_list[mn - 1])):
        con2 = model_tree.addConstr(
            sum(np.sum(a_div3[mn - 1][i][j]) for j in range(0, len(df_list[mn])))
            + np.sum(a_div2[mn - 1][i])
            + a_exit[mn - 1][i]
            - s[mn - 1][i]
            == 0
        )

    if mn > 1:
        aac2a = [LinExpr(0)] * len(df_list[mn - 1])
        aac3a = [LinExpr(0)] * len(df_list[mn - 1])
        aac4a = [LinExpr(0)] * len(df_list[mn - 1])
        for i in range(0, len(df_list[mn - 1])):
            aac2a[i] = sum(a_div2[mn - 2][k][i] for k in range(0, len(df_list[mn - 2])))
            aac3a[i] = sum(
                a_div3[mn - 2][k][i][j] + a_div3[mn - 2][k][j][i]
                for j in range(0, len(df_list[mn - 1]))
                for k in range(0, len(df_list[mn - 2]))
            )
            con2a = model_tree.addConstr(
                aac2a[i]
                + aac3a[i]
                - (
                    sum(
                        np.sum(a_div3[mn - 1][i][j]) for j in range(0, len(df_list[mn]))
                    )
                    + np.sum(a_div2[mn - 1][i])
                    + a_exit[mn - 1][i]
                )
                == 0
            )

    # model_tree.update()
    imgs_mat[mn] = imgs
    th_mat[mn] = th
    print(mn)

model_tree.update()
model_tree.optimize()
model_tree.write("file.lp")


s_arr = []

for v in model_tree.getVars():
    if v.varName[0] == "n":
        s_arr += [v.x]

df_list_concat = pd.concat([df_list[mn] for mn in range(0, len(arr_m))])
df_list_concat["s_i"] = s_arr
df_list_concat = df_list_concat.reset_index(drop=True)

x_df = 0
for jk in range(0, len(arr_m)):
    for i in range(0, len(df_list[jk])):
        df_list[jk]["s_i"][i] = df_list_concat["s_i"][i + x_df]
    x_df += len(df_list[jk])


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


for i_df in range(0, len(arr_m)):
    df_a = df_list[i_df]
    df_new = df_a[df_a["s_i"] == 1]
    df_new = df_new.reset_index(drop=True)
    img_mask_video = get_img_we_want2(i_df, df_new)
    img1 = img_mask_video.save("%d.png" % i_df)


df_list_concat = df_list_concat.reset_index(drop=False)

for j in range(len(arr_m)):
    for i in range(len(df_list_concat)):
        if df_list_concat["t"][i] == j:
            df_list_concat["index"][i] = df_list_concat["index"][i] - sum(
                len(df_list[m]) for m in range(0, j)
            )


df_list_concat = df_list_concat[df_list_concat["s_i"] == 1]


asg_arr = []
for v in model_tree.getVars():
    if v.varName[0] == "x" and v.x == 1:
        asg_arr += [v.varName]
        print("%s %g" % (v.varName, v.x))

for v in model_tree.getVars():
    if v.varName[0] == "a" and v.x == 1:
        asg_arr += [v.varName]
        print("%s %g" % (v.varName, v.x))


for v in model_tree.getVars():
    if v.varName[0] == "n" and v.x == 1:
        print("%s %g" % (v.varName, v.x))

for i in range(len(df_new)):
    asg_arr += [0]


for i in range(len(asg_arr)):
    asg_arr[i] = str(asg_arr[i])
    asg_arr[i] = asg_arr[i].split("_")


df_list_concat["MoD"] = [0] * len(df_list_concat)
df_list_concat["assign"] = asg_arr
df_list_concat = df_list_concat.reset_index(drop=False)


for i in range(len(df_list_concat)):
    if len(df_list_concat["assign"][i]) == 4:
        df_list_concat["MoD"][i] = "move"
    if len(df_list_concat["assign"][i]) == 5:
        df_list_concat["MoD"][i] = "divide"
    if len(df_list_concat["assign"][i]) == 1:
        df_list_concat["MoD"][i] = "last"


import cv2
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

# video = cv2.VideoWriter('video_mask27.mp4',-1,2,(width,height)) #second last is fps should be a factor of total time points

# for i in range (0,len(arr_m)):

#    video.write(images[i])


# cv2.destroyAllWindows()
# video.release()


concatenated2 = np.concatenate([x for x in images[0 : len(arr_m)]], axis=0)


df_list_concat["centers_which"] = [0] * len(df_list_concat)
for i in range(0, len(df_list_concat)):
    if df_list_concat["MoD"][i] == "move":
        df_list_concat["centers_which"][i] = df_list_concat["assign"][i][-1:]
    if df_list_concat["MoD"][i] == "divide":
        df_list_concat["centers_which"][i] = df_list_concat["assign"][i][-2:]


df_list_concat["center_value"] = [[]] * len(df_list_concat)
for i in range(0, len(df_list_concat)):
    df_list_concat["center_value"][i] = []
    for k in range(0, len(df_list_concat)):
        if type(df_list_concat["centers_which"][i]) != int:
            for j in range(len(df_list_concat["centers_which"][i])):
                if (df_list_concat["t"][k] == df_list_concat["t"][i] + 1) and (
                    df_list_concat["index"][k]
                    == int(df_list_concat["centers_which"][i][j])
                ):
                    df_list_concat["center_value"][i] += [
                        str((df_list_concat["center"][k]))
                    ]


df_list_concat.to_csv("df_label_list_new_yay.csv")
for i in range(0, len(arr_m) - 1):
    df_see = df_list_concat[df_list_concat["t"] == i]
    df_see = df_see.reset_index(drop=True)
    print(df_see)
    for k in range(len(df_see)):
        print(k, len(df_see["center_value"][k]), df_see["center_value"][k])
        for j in range(len(df_see["center_value"][k])):
            x1 = int(df_see["center"][k])
            x2 = int(float(df_see["center_value"][k][j]))
            concatenated2 = cv2.line(
                concatenated2,
                (x1, int(38.5) + height * (i)),
                (x2, int(38.5 + (height * (i + 1)))),
                (255, 0, 0),
                5,
            )


concatenated2 = Image.fromarray(concatenated2)
concatenated2.save("test_img8dac.jpg")
