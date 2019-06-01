import skimage
from PIL import Image, ImageDraw, ImageFilter
from skimage import morphology
import cv2
import numpy as np
import matplotlib.pyplot as plt
import scipy

img_ini = cv2.imread("0.png", 0)
###########################################img = cv2.medianBlur(img,5)

th = [[]] * 3
imgs = [[]] * 3

ret, th[0] = cv2.threshold(img_ini, 170, 255, cv2.THRESH_BINARY)
reta, th[1] = cv2.threshold(img_ini, 180, 255, cv2.THRESH_BINARY)
reta, th[2] = cv2.threshold(img_ini, 185, 255, cv2.THRESH_BINARY)


for i in range(3):
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

    return tree


def tree(node1, node2):
    tree1 = get_tree(n_mask_img[node1], arr[node1], n_mask_img[node2], arr[node2])

    return tree1


tree2 = tree(0, 1)
tree3 = tree(1, 2)
print(tree2, tree3)


# print(img[0]==2)
# img3 = np.where(imgs[0]==2, th[0], 7)
# Image.fromarray(img3).show()
# print(img)


def get_coord_min_max(threshold_number, mask_number):
    mask_number = mask_number - 1
    array1 = arr[threshold_number][mask_number]
    max_y = max(array1, key=lambda item: item[0])  # y coordinate
    min_y = min(array1, key=lambda item: item[0])  # y coordinate
    max_x = max(array1, key=lambda item: item[1])  # y coordinate
    min_x = min(array1, key=lambda item: item[1])  # y coordinate

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
    arr_get, arr_hist2, img2 = get_coord_min_max(threshold_number, mask_number)

    # print(arr_hist2)

    number_of_peaks = scipy.signal.find_peaks(arr_hist2)
    length = arr_get[2] - arr_get[3]

    # print(len(number_of_peaks[0]), length)
    cost = [len(number_of_peaks[0]), length]

    costs = cost[0] + (cost[1] // 120)

    return cost, costs


img2a = get_img(0, 2)
cost_get, costs_get = get_cost(0, 2)
print(cost_get, costs_get)


for i in range(0, len(imgs)):
    for j in range(0, (n_mask_img[i])):
        cost_get, costs_get = get_cost(i, j)
        print(costs_get, i, j)


####### gurobi for [1: 2 and 3,   2: 4 and 5,  3: 6, 7 and 8]
from gurobipy import *

model_tree = Model()

c = [22, 3, 15, 22, 21, 2, 2, 3]
s = [0] * len(c)
for i in range(0, len(c)):
    s1 = model_tree.addVar(0, 1, obj=c[i], vtype=GRB.INTEGER, name="node1")
    s[i] = s1

    model_tree.update()
con1 = model_tree.addConstr(s[3] + s[1] + s[0] == 1)
con2 = model_tree.addConstr(s[4] + s[1] + s[0] == 1)
con3 = model_tree.addConstr(s[5] + s[2] + s[0] == 1)
con4 = model_tree.addConstr(s[6] + s[2] + s[0] == 1)
con5 = model_tree.addConstr(s[7] + s[2] + s[0] == 1)

model_tree.optimize()
print("nodes", s[0].X, s[1].X, s[2].X, s[3].X, s[4].X, s[5].X, s[6].X, s[7].X)
