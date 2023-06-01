import cv2
import numpy as np
import skimage
from PIL import Image, ImageDraw, ImageFilter
from skimage import morphology

img = cv2.imread("0.png", 0)
###########################################img = cv2.medianBlur(img,5)

th = [[]] * 3
imgs = [[]] * 3

ret, th[0] = cv2.threshold(img, 170, 255, cv2.THRESH_BINARY)
reta, th[1] = cv2.threshold(img, 180, 255, cv2.THRESH_BINARY)
reta, th[2] = cv2.threshold(img, 185, 255, cv2.THRESH_BINARY)


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

# arr is the 2 X 2 matrix of image number and mask numer : arr[img_number][mask_number] last is arr[2][6]


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


tree2 = tree(1, 2)
