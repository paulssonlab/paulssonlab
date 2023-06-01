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


arr_m = [0] * 70
df_list = [pd.DataFrame()] * len(arr_m)
df_label_list = [pd.DataFrame()] * len(arr_m)

from gurobipy import *

model_tree = Model()
# a=0
s = [[-1] * 1000] * len(arr_m)
# constr_mat=[[model_tree.addConstr(a==0)]*100] *len(arr_m)

imgs_mat = [[[]]] * len(arr_m)
th_mat = [[[]]] * len(arr_m)

labels_arr = [[]] * len(arr_m)
ass_arr = [[[]]] * len(arr_m)
center = [[]] * len(arr_m)
len_x_m = [[]] * len(arr_m)
for mn in range(len(arr_m)):
    img_ini = cv2.imread("%d.jpg" % (mn))  ###########################################
    #

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
                        # print(l+1,i+1)
                        tree += [[l + 1, i + 1]]
        res = []
        for k in tree:
            if k not in res:
                res.append(k)
        return res  # tree

    def tree(node1, node2):
        tree1 = get_tree(n_mask_img[node1], arr[node1], n_mask_img[node2], arr[node2])

        return tree1

    """tree2= tree(0,2)
    #tree3 = tree(1,2)
    print(tree2)#tree3)"""

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

        # print(arr_hist2[:,0])

        number_of_peaks = scipy.signal.find_peaks(arr_hist2[:, 0])
        length = arr_get[2] - arr_get[3]

        # print(len(number_of_peaks[0]), length)
        cost = [len(number_of_peaks[0]), length]

        if length < 10 or length == 0:
            cost[1] == 1000
        costs = cost[0] * 5 + abs((cost[1] - 100))  # /10)
        # costs = cost[0] + abs((cost[1]//120))
        # print('a',threshold_number,mask_number)

        return cost, costs

    # img2a = get_img(0,2)
    """cost_get, costs_get = get_cost(1,2)
    print(cost_get, costs_get)"""

    costs_get_arr = []
    costs_get_arr2 = []
    df = pd.DataFrame()

    for i in range(0, len(imgs)):
        for j in range(0 + 1, (n_mask_img[i]) + 1):
            cost_get, costs_get = get_cost(i, j)

            costs_get_arr += [costs_get]
            costs_get_arr2 += [[costs_get, i, j]]

    df["cost"] = [x[0] for x in costs_get_arr2]
    df["threshold_number"] = [x[1] for x in costs_get_arr2]
    df["mask"] = [x[2] for x in costs_get_arr2]
    df["s_i"] = [0 for x in costs_get_arr2]
    df["s_i_assign"] = [0 for x in costs_get_arr2]
    df_list[mn] = df

    c = costs_get_arr

    for i in range(0, len(df_list[mn])):
        s1 = model_tree.addVar(
            0, 1, obj=c[i], vtype=GRB.INTEGER, name="node%d,%d" % (i, mn)
        )
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
        # con1 = model_tree.addConstr(s[0]+s[(tree(1,3)[int(np.where([x[0]==i for x in tree(1,3)])[0])][1])]+s[2+(tree(2,3)[int(np.where([x[0]==i for x in tree(2,3)])[0])][1])] +s[n_constr+i] == 1)
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
        # constr_mat[mn][i] = con1
        model_tree.update()
    model_tree.update()
    model_tree.optimize()

    s_arr_label = []
    for v in model_tree.getVars():
        # print('%s %g' % (v.varName, v.x),mn, )

        # v.varName.split(',')
        if int(v.varName.split(",")[-1]) == mn:
            s_arr_label += [v.x]

    constrs = model_tree.getConstrs()

    imgs_mat[mn] = imgs
    th_mat[mn] = th

    df["s_i_assign"] = s_arr_label
    count = 0
    for im in range(len(s_arr_label)):
        if s_arr_label[im] == 1:
            count += 1
    labels_arr[mn] = [i for i in range((count))]
    ass_arr[mn] = [""] * len(labels_arr[mn - 1])

    df_assign = df[df["s_i_assign"] == 1]
    df_assign = df_assign.reset_index(drop=True)

    center[mn] = []
    len_x_m[mn] = []
    center_len_x = []
    for i in range(len(labels_arr[mn])):
        center[mn] += [
            get_center(df_assign["threshold_number"][i], df_assign["mask"][i])
        ]
        len_x_m[mn] += [
            get_lenx(df_assign["threshold_number"][i], df_assign["mask"][i])
        ]
        center_len_x += [(center[mn][i], len_x_m[mn][i])]
    # center_len_x2 = center_len_x
    # len_check = len_x_m[mn]
    # center_check  = center[mn]
    center_len_x.sort()
    center[mn] = []
    len_x_m[mn] = []
    for i in range(len(labels_arr[mn])):
        center[mn] += [center_len_x[i][0]]
        len_x_m[mn] += [center_len_x[i][1]]
        # center_len_x += [(center[mn][i] ,len_x_m[mn][i]) ]

    # center_len_x.sort(key=lambda x: x[0])
    # center[mn] = center_len_x[0]
    # len_x_m[mn] = center_len_x[1]

    ass_arr[0] = "start"

    count_d2 = 0
    for i in range(len(labels_arr[mn - 1])):
        print(mn, i, "popopoo")

        # if (i+count_d2) > len(labels_arr[mn])-1:
        #   count_d2 = 0
        if (
            len_x_m[mn - 1][i] - len_x_m[mn][i + count_d2]
        ) > 53:  # for i in range (len(labels_arr[mn-1]))):
            # for i in range(len(labels_arr[mn-1])):
            ass_arr[mn][i] = "divide"
            count_d2 += 1
        # if (i+count_d2) > len(labels_arr[mn])-1:
        #   count_d2 = 0
        if len(labels_arr[mn]) < len(
            labels_arr[mn - 1]
        ):  # and (abs(len_x_m[mn][i+count_d2] - len_x_m[mn-1][i]) < 40 )):
            # for i in range(len(labels_arr[mn-1])):
            if i < len(labels_arr[mn]):
                ass_arr[mn][i] = "move"
            else:
                ass_arr[mn][i] = "exit"

        print((center[mn][0]), i + count_d2)
        print(abs(center[mn][i + count_d2] - center[mn - 1][i]))
        if abs(center[mn][i + count_d2] - center[mn - 1][i]) < 50 and (
            (len_x_m[mn][i + count_d2] >= len_x_m[mn - 1][i])
            or (len_x_m[mn - 1][i] - len_x_m[mn][i + count_d2] < 54)
        ):  # for i in range (len(labels_arr[mn-1])))) :
            # for i in range(len(labels_arr[mn-1])):
            ass_arr[mn][i] = "move"

        if (len_x_m[mn][i + count_d2] - len_x_m[mn - 1][i]) < 150 and (
            len_x_m[mn][i + count_d2] - len_x_m[mn - 1][i]
        ) > 80:  # for i in range (len(labels_arr[mn-1])))) :
            ass_arr[mn][i] = "move"

    df_label = pd.DataFrame()
    df_label["label"] = [labels_arr[mn]]
    df_label["assg"] = [ass_arr[mn]]
    df_label["center"] = [center[mn]]
    df_label["lenx"] = [len_x_m[mn]]
    df_label_list[mn] = df_label


model_tree.optimize()
constrs_mat = model_tree.getConstrs()
model_tree.write("file.lp")

s_arr = []
for v in model_tree.getVars():
    # print('%s %g' % (v.varName, v.x))
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
    print("cc", len(imgsa), len(tha))
    threshold_number = df_img["threshold_number"][0]
    mask_number = df_img["mask"][0]
    img3 = np.where(
        imgsa[threshold_number] == mask_number,
        tha[threshold_number],
        tha[threshold_number] * 0,
    )  # this prints only the areas which are labeled as img == label #img1a==5 img1==2

    if len(df_img) > 1:
        for i in range(1, len(df_img)):
            threshold_number = df_img["threshold_number"][i]
            mask_number = df_img["mask"][i]
            img3 = img3 + np.where(
                imgsa[threshold_number] == mask_number,
                tha[threshold_number],
                tha[threshold_number] * 0,
            )  # this prints only the areas which are labeled as img == label #img1a==5 img1==2

    # else:
    # print(i_df,'bb',threshold_number,mask_number)

    return Image.fromarray(img3)  # , Image.fromarray(img3).show()


for i_df in range(0, len(arr_m)):
    df_a = df_list[i_df]
    df_new = df_a[df_a["s_i"] == 1]
    df_new = df_new.reset_index(drop=True)

    img_mask_video = get_img_we_want2(i_df, df_new)
    img1 = img_mask_video.save("%d.png" % i_df)


#####df_label_list[6]['assg'][0][1] = 'move'###########################################

for ik in range(0, len(arr_m)):
    lwn = len(df_label_list[ik]["label"][0])
    count_d = 0
    for k in range(0, len(df_label_list[ik]["assg"][0])):
        print(ik, k, df_label_list[ik]["assg"][0])

        if df_label_list[ik]["assg"][0][k] == "divide":
            if k + count_d < len(df_label_list[ik]["assg"][0]) + count_d:
                df_label_list[ik]["label"][0][k + count_d] = Make_d_nodes(
                    [df_label_list[ik - 1]["label"][0][k]]
                )[0]
                df_label_list[ik]["label"][0][k + 1 + count_d] = Make_d_nodes(
                    [df_label_list[ik - 1]["label"][0][k]]
                )[1]

            elif k + count_d == len(df_label_list[ik]["assg"][0]) + count_d:
                df_label_list[ik]["label"][0][k + count_d] = Make_d_nodes(
                    [df_label_list[ik - 1]["label"][0][k]]
                )[0]

            count_d += 1

        if df_label_list[ik]["assg"][0][k] == "move":
            df_label_list[ik]["label"][0][k + count_d] = df_label_list[ik - 1]["label"][
                0
            ][k]

        if df_label_list[ik]["assg"][0][k] == "exit":
            df_label_list[ik]["label"][0] = df_label_list[ik]["label"][0][0:k]


df_label_list_concat = pd.concat([df_label_list[mn] for mn in range(0, len(arr_m))])
df_label_list_concat.to_csv("df_label_list.csv")


import glob

import cv2

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
    "video_mask_testc.mp4", -1, 2, (width, height)
)  # second last is fps should be a factor of total time points

for i in range(0, len(arr_m)):
    video.write(images[i])


cv2.destroyAllWindows()
video.release()
