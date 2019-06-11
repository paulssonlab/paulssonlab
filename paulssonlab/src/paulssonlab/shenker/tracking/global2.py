import skimage
from PIL import Image  # , ImageDraw, ImageFilter
import cv2
import numpy as np
import matplotlib.pyplot as plt
import scipy
import pandas as pd
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


arr_m = [0] * 5
df_list = [pd.DataFrame()] * len(arr_m)
df_label_list = [pd.DataFrame()] * len(arr_m)

from gurobipy import *

model_tree = Model()

s = [[-1] * 1000] * len(arr_m)
a_div = [[-1] * 1000] * len(arr_m)
a_exit = [[-1] * 1000] * len(arr_m)
a_move = [[-1] * 1000] * len(arr_m)
imgs_mat = [[[]]] * len(arr_m)
th_mat = [[[]]] * len(arr_m)

labels_arr = [[]] * len(arr_m)
ass_arr = [[[]]] * len(arr_m)
center = [[]] * len(arr_m)
len_x_m = [[]] * len(arr_m)
for mn in range(len(arr_m)):

    img_ini = cv2.imread("%d.jpg" % (mn + 34))
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
    df["a_exit"] = [0 for x in costs_get_arr2]
    df["a_move"] = [0 for x in costs_get_arr2]
    df["a_div"] = [0 for x in costs_get_arr2]
    # df['s_i_assign'] = [0 for x in costs_get_arr2]
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

    d = costs_get_arr * 4
    for i in range(0, len(df_list[mn])):
        a1 = model_tree.addVar(
            0, 1, obj=d[i], vtype=GRB.INTEGER, name="a_div%d,%d" % (i, mn)
        )
        a_div[mn][i] = a1
        model_tree.update()

    e = costs_get_arr * 3
    for i in range(0, len(df_list[mn])):
        a1 = model_tree.addVar(
            0, 1, obj=e[i], vtype=GRB.INTEGER, name="a_exit%d,%d" % (i, mn)
        )
        a_exit[mn][i] = a1
        model_tree.update()

    f = costs_get_arr * 2
    for i in range(0, len(df_list[mn])):
        a1 = model_tree.addVar(
            0, 1, obj=f[i], vtype=GRB.INTEGER, name="a_move%d,%d" % (i, mn)
        )
        a_move[mn][i] = a1
        model_tree.update()

    for i in range(0, len(df_list[mn])):

        # con2 = model_tree.addConstr(a_div[mn][i] + a_move[mn][i] + a_exit[mn][i] +s[mn][i]<= 2 )
        con2 = model_tree.addConstr(
            a_div[mn][i] + a_move[mn][i] + a_exit[mn][i] - s[mn][i] == 0
        )

        model_tree.update()

    for i in range(0, len(df_list[mn])):

        con3 = model_tree.addConstr(a_div[mn][i] + a_move[mn][i] + a_exit[mn][i] <= 1)
        model_tree.update()

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

    imgs_mat[mn] = imgs
    th_mat[mn] = th

model_tree.optimize()
constrs_mat = model_tree.getConstrs()
model_tree.write("file.lp")


s_arr = []  # [0] * (len(df_list_concat)*2)

for v in model_tree.getVars():
    # print('%s %g' % (v.varName, v.x))
    if v.varName[0] == "n":  # .split(',')[-1])==mn):
        print("%s %g" % (v.varName, v.x))
        s_arr += [v.x]

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
# df_new = df[df['s_i']==1]
# df_new = df_new.reset_index(drop=True)
# print(df_new)


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

    return Image.fromarray(img3)  # , Image.fromarray(img3).show()


# img_mask_video, img_mask_video_a= get_img_we_want2(df_new)
for i_df in range(0, len(arr_m)):
    df_a = df_list[i_df]
    df_new = df_a[df_a["s_i"] == 1]
    df_new = df_new.reset_index(drop=True)
    # print(i_df,df_new)
    img_mask_video = get_img_we_want2(i_df, df_new)
    img1 = img_mask_video.save("%d.png" % i_df)
    # print('i_df', i_df)

a_exit_arr = []
for v in model_tree.getVars():
    if v.varName[3] == "x":
        a_exit_arr += [v.x]
        # print('%s %g' % (v.varName, v.x))

a_div_arr = []
for v in model_tree.getVars():
    if v.varName[3] == "i":
        a_div_arr += [v.x]
        # print('%s %g' % (v.varName, v.x))

a_move_arr = []
for v in model_tree.getVars():
    if v.varName[2] == "m":
        a_move_arr += [v.x]
        # print('%s %g' % (v.varName, v.x),a_move_arr)


df_list_concat["a_exit"] = a_exit_arr
df_list_concat["a_div"] = a_div_arr
df_list_concat["a_move"] = a_move_arr
# df_list_concat = a_move_arr.reset_index(drop=True)


# df_list_concat = df_list_concat[df_list_concat['s_i']==1]
x_df = 0
for jk in range(0, len(arr_m)):
    for i in range(0, len(df_list[jk])):
        df_list[jk]["s_i"][i] = df_list_concat["s_i"][i + x_df]
        df_list[jk]["a_exit"][i] = df_list_concat["a_exit"][i + x_df]
        df_list[jk]["a_div"][i] = df_list_concat["a_div"][i + x_df]
        df_list[jk]["a_move"][i] = df_list_concat["a_move"][i + x_df]
    x_df += len(df_list[jk])

for jk in range(0, len(arr_m)):
    df_list[jk] = df_list[jk][df_list[jk]["s_i"] == 1]
    df_list[jk] = df_list[jk].reset_index(drop=True)

labels_arr = [[]] * (len(arr_m))
ass_arr = [[]] * (len(arr_m))
ass_arr[0] = ["start"]

for jk in range(0, len(arr_m)):
    for i in range(0, len(df_list[jk])):
        labels_arr[jk] = [i for i in range(len(df_list[jk]))]
        ass_arr[jk] = [""] * (len(df_list[jk - 1]))

for jk in range(1, len(arr_m)):
    for i in range(0, len(df_list[jk - 1])):
        print(jk, i, "popopoo")
        # ass_arr[jk][i] = 'move'
        if df_list[jk]["a_exit"][i] == 1:
            ass_arr[jk][i] = "exit"

        if df_list[jk]["a_div"][i] == 1:
            ass_arr[jk][i] = "divide"

        if df_list[jk]["a_move"][i] == 1:
            ass_arr[jk][i] = "move"


df_label_list = [pd.DataFrame()] * len(arr_m)


for mn in range(0, len(arr_m)):
    df_label = pd.DataFrame()
    df_label["label"] = [labels_arr[mn]]
    df_label["assg"] = [ass_arr[mn]]

    df_label_list[mn] = df_label


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
df_label_list_concat.to_csv("df_label_list_new.csv")


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

video = cv2.VideoWriter(
    "video_mask_testd.mp4", -1, 2, (width, height)
)  # second last is fps should be a factor of total time points

for i in range(0, len(arr_m)):

    video.write(images[i])


cv2.destroyAllWindows()
video.release()
