from operator import itemgetter

import cv2
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import skimage
from PIL import Image, ImageDraw, ImageFilter
from skimage import morphology

for m in range(0, 75):
    # m=12
    # print(i)
    # m=i
    img_ini = cv2.imread("%d.jpg" % m)

    th = [[]] * 4
    imgs = [[]] * 4

    reta, th[0] = cv2.threshold(img_ini, 150, 255, cv2.THRESH_BINARY)
    ret, th[1] = cv2.threshold(img_ini, 170, 255, cv2.THRESH_BINARY)
    reta, th[2] = cv2.threshold(img_ini, 180, 255, cv2.THRESH_BINARY)
    reta, th[3] = cv2.threshold(img_ini, 185, 255, cv2.THRESH_BINARY)

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

        costs = cost[0] + abs((cost[1] - 120) / 10)
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
            print(costs_get, i, j)

    df["cost"] = [x[0] for x in costs_get_arr2]
    df["threshold_number"] = [x[1] for x in costs_get_arr2]
    df["mask"] = [x[2] for x in costs_get_arr2]

    #####use get_cost(1,1) for mask 1 for threshold 2
    ####### gurobi for [1: 2 and 3,   2: 4 and 5,  3: 6, 7 and 8]
    from gurobipy import *

    model_tree = Model()

    c = costs_get_arr
    s = [0] * len(c)
    for i in range(0, len(c)):
        s1 = model_tree.addVar(0, 1, obj=c[i], vtype=GRB.INTEGER, name="node1")
        s[i] = s1

        model_tree.update()

    n_constr = n_mask_img[len(th) - 1]
    print(n_constr)

    n_constr_c = (
        n_mask_img[len(th) - 2] + n_mask_img[len(th) - 3] + n_mask_img[len(th) - 4] - 1
    )
    print(n_constr_c)
    n_constr_a = n_mask_img[len(th) - 3] + n_mask_img[len(th) - 4] - 1
    print(n_constr_a)
    n_constr_b = n_mask_img[len(th) - 4] - 1
    print(n_constr_b)

    constr = [model_tree.addConstr(s1 == 0)] * n_constr

    for i in range(1, n_constr + 1):
        print(i)
        # con1 = model_tree.addConstr(s[0]+s[(tree(1,3)[int(np.where([x[0]==i for x in tree(1,3)])[0])][1])]+s[2+(tree(2,3)[int(np.where([x[0]==i for x in tree(2,3)])[0])][1])] +s[n_constr+i] == 1)
        con1 = model_tree.addConstr(
            s[(tree(0, 3)[int(np.where([x[0] == i for x in tree(0, 3)])[0])][1])]
            + s[
                n_constr_b
                + (tree(1, 3)[int(np.where([x[0] == i for x in tree(1, 3)])[0])][1])
            ]
            + s[
                n_constr_a
                + (tree(2, 3)[int(np.where([x[0] == i for x in tree(2, 3)])[0])][1])
            ]
            + s[n_constr_c + i]
            == 1
        )
        # s[(tree(0,3)[int(np.where([x[0]==i for x in tree(0,3)])[0])][1])]

        constr[i - 1] = con1

    model_tree.optimize()

    s_arr = [] * len(c)
    for i in range(len(c)):
        s_a = s[i].X
        s_arr += [s_a]
        print(i, s[i].X)

    df["s_i"] = s_arr

    df_new = df[df["s_i"] == 1]
    df_new = df_new.reset_index(drop=True)
    print(df_new)

    """def get_img_we_want(threshold_number,mask_number):
        img3 = np.where(imgs[threshold_number]==mask_number, th[threshold_number],th[threshold_number]*0)   # this prints only the areas which are labeled as img == label #img1a==5 img1==2
        img4 = np.where(imgs[threshold_number]==mask_number+1, th[threshold_number],th[threshold_number]*0)
        img  = img3+img4
        return Image.fromarray(img).show()"""

    # img= get_img_we_want(1,1)

    def get_img_we_want2(df_img):
        threshold_number = df_img["threshold_number"][0]
        mask_number = df_img["mask"][0]
        img3 = np.where(
            imgs[threshold_number] == mask_number,
            th[threshold_number],
            th[threshold_number] * 0,
        )  # this prints only the areas which are labeled as img == label #img1a==5 img1==2

        for i in range(1, len(df_img)):
            threshold_number = df_img["threshold_number"][i]
            mask_number = df_img["mask"][i]
            img3 = img3 + np.where(
                imgs[threshold_number] == mask_number,
                th[threshold_number],
                th[threshold_number] * 0,
            )  # this prints only the areas which are labeled as img == label #img1a==5 img1==2
            print(threshold_number, mask_number)

        # img4 = np.where(imgs[threshold_number]==mask_number+1, th[threshold_number],th[threshold_number]*0)
        # img  = img3+img4

        return Image.fromarray(img3)  # , Image.fromarray(img3).show()

    # img_mask_video, img_mask_video_a= get_img_we_want2(df_new)
    img_mask_video = get_img_we_want2(df_new)
    img1 = img_mask_video.save("%d.png" % m)
    print("m", m)


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
    # print (img)


img0 = cv2.imread("0.png")
height, width, layers = img0.shape

video = cv2.VideoWriter(
    "video_mask_test3.mp4", -1, 2, (width, height)
)  # second last is fps should be a factor of total time points

for i in range(0, 75):
    video.write(images[i])


cv2.destroyAllWindows()
video.release()
