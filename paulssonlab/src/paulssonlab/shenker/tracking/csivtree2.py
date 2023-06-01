import math
import random

import cv2
import numpy as np
import pandas as pd
import skimage
from PIL import Image, ImageDraw, ImageFilter

bright_halo = "rgb(68,68,68)"
bright_blur = "rgb(128,128,128)"
bright = "rgb(225,225,225)"
pixelsize = 50
g_blur_r = 1
t = 70
list2 = [2, 0, 1, -2, -1]
pixelsize_list = [100, 50, 80, 30, 40]
bright_halo_list = ["rgb(68,68,68)", "rgb(48,48,48)", "rgb(88,88,88)"]
bright_blur_list = ["rgb(0,0,0)", "rgb(60,60,60)", "rgb(40,40,40)"]
bright_list = ["rgb(215,215,215)", "rgb(255,255,255)", "rgb(240,240,240)"]


# dr = ImageDraw.Draw(img1)
#       lines = dr.line((1800, 0, 1800, 150), fill='red')


def make_ellipse(image, x, y, dx, dy, next_1, color):
    dr = ImageDraw.Draw(image)
    ell = dr.ellipse((x, y, dx, dy), fill=color)
    ell2 = dr.ellipse((x + next_1, y, dx + next_1, dy), color)
    img1 = Image.new("RGB", (next_1, dy - y + 1), color=color)
    image.paste(img1, ((x + int((dx - x) / 2), int(y))))
    return ell, ell2, img1


def halo(image, x, y, dx, dy, color):
    dr = ImageDraw.Draw(image)
    ell = dr.ellipse((x, y, dx, dy), fill=color)
    return ell


def draw_ellipse(length):
    bright_halo = random.choice(bright_halo_list)
    bright_blur = random.choice(bright_blur_list)
    bright = random.choice(bright_list)
    img = Image.new(
        "RGBA", (800, 160), "rgba(0,0,0,0)"
    )  # we do get a background less image
    ell1 = make_ellipse(img, x, y, dx, dy, length, bright_blur)
    ell1 = make_ellipse(img, x + 20, y + 4, dx - 20, dy - 4, length, bright)  # x7-7
    return img


def image_edit(img, angle):
    pixelsize = random.choice(pixelsize_list)
    img = img.rotate(angle)
    imgSmall = img.resize((pixelsize, pixelsize), resample=Image.BILINEAR)
    result = imgSmall.resize(img.size, Image.NEAREST)
    result = result.filter(ImageFilter.GaussianBlur(g_blur_r))

    return result


def image_edit2(img):
    img = img.rotate(0)
    imgSmall = img.resize((150, 150), resample=Image.BILINEAR)
    result = imgSmall.resize(img.size, Image.NEAREST)
    result = result.filter(ImageFilter.GaussianBlur(3))

    return result


def position_ellipse(imga, length_mat, x_pos_mat, y_pos_mat, label):
    for i in range(len(length_mat)):
        img = draw_ellipse(int(length_mat[i]))
        result = image_edit(img, random.choice(list2))
        imga.paste(result, (x_pos_mat[i], y_pos_mat[i]), result)

    imga = image_edit2(imga)  ###

    imga = imga.save("%d.jpg" % (label - 0))
    return imga


def pos_mat_gen(len_mat_in):
    pos_mat2 = [0] * len(len_mat_in)
    pos_mat2[0] = 0
    for i in range(1, len(len_mat_in)):
        pos_mat2[i] = int(
            (len_mat_in[i - 1] + pos_mat2[i - 1]) + (160 * (len_mat_in[i - 1]) / 450)
        )  ################

    # print (pos_mat2)
    return pos_mat2


df = pd.DataFrame()

m1_mat = [np.array([None] * 100000)] * t
m1_mat[0] = np.array([82, 88, 74])
m1_mat[0] = np.array(m1_mat[0][m1_mat[0] != None])

df["t"] = [0] * t
df["len_mat"] = [m1_mat[0]] * t
df["right_edge"] = [m1_mat[0]] * t
df["center"] = [m1_mat[0]] * t
sum_len = [[]] * t

for i in range(1, t):
    df["t"][i] = i
    x = 125
    y = 22
    dx = x + 75
    dy = y + 33

    sum_len[i - 1] = [0] * len(m1_mat[i - 1])

    for j in range(len(m1_mat[i - 1])):
        if m1_mat[i - 1][j] <= 180:
            m1_mat[i][j] = int(m1_mat[i - 1][j] + random.uniform(10, 20))

        elif m1_mat[i - 1][j] > 180:
            m1_mat[i][j] = int(m1_mat[i - 1][j] / 2)

        sum_len[i - 1][j] = sum(m1_mat[i - 1][k] for k in range(0, j + 1))

    count = 0
    for j in range(len(m1_mat[i - 1])):
        if type(m1_mat[i][j + count]) == int and type(m1_mat[i - 1][j]) == int:
            if (
                (int(m1_mat[i - 1][j]) == 2 * int(m1_mat[i][j + count]))
                or (int(m1_mat[i - 1][j]) == (2 * (m1_mat[i][j + count]) - 1))
                or (int(m1_mat[i - 1][j]) == (2 * (m1_mat[i][j + count]) + 1))
            ):
                m1_mat[i] = list(m1_mat[i])
                m1_mat[i].insert(j + count, int(m1_mat[i - 1][j] / 2))
                m1_mat[i] = np.asarray(m1_mat[i])

                count += 1

    arr = []
    for j in range(len(m1_mat[i - 1])):
        if (sum_len[i - 1][j]) > 1200:
            arr += [j]

    # print(i-1,arr)
    # print(m1_mat[i-1] )
    m1_mat[i - 1] = np.delete(m1_mat[i - 1], arr)
    # print(m1_mat[i-1] ,'>')
    """for k in arr:

        print(sum_len[i-1][k])"""

    sum_len[i - 1] = [0] * len(m1_mat[i - 1])
    for j in range(len(m1_mat[i - 1])):
        sum_len[i - 1][j] = sum(m1_mat[i - 1][k] for k in range(0, j + 1))

    df["right_edge"][i - 1] = sum_len[i - 1]
    df["center"][i - 1] = [0] * len(m1_mat[i - 1])

    m1_mat[i] = np.array(m1_mat[i][m1_mat[i] != None])
    df["len_mat"][i - 1] = m1_mat[i - 1]

    for j in range(len(m1_mat[i - 1])):
        df["center"][i - 1][j] = np.asarray(df["right_edge"][i - 1][j]) - (
            np.asarray(df["len_mat"][i - 1][j]) / 2
        )


df["assign"] = [[""] * 50] * t

mat = []
for i in range(1, t):
    mat = df[df["t"] == (i - 1)]["len_mat"].values
    len_m = len(mat[0])
    mat2 = [""] * len_m

    for j in range(len_m):
        # print( (len(df['len_mat'][i-1]) , len(df['len_mat'][i])  ))
        count22 = 0
        count33 = 0
        if len(df["len_mat"][i - 1]) <= len(df["len_mat"][i]):
            if df["len_mat"][i - 1][j] <= 180:  # df['len_mat'][i][j+count22]:
                mat2[j] = "move"
            else:
                mat2[j] = "divide"
                count22 += 1

        else:
            if j < len(df["len_mat"][i]):
                if df["len_mat"][i - 1][j] <= 180:  # df['len_mat'][i][j+count33]:
                    mat2[j] = "move"
                else:
                    mat2[j] = "divide"
                    count33 += 1
            elif j >= len(df["len_mat"][i]):
                mat2[j] = "exit"

    df["assign"][i - 1] = mat2


df["from"] = [[None] * 50] * t
for i in range(0, t):
    df["from"][i] = df["center"][i]


df["re2"] = [[None] * 50] * t
for i in range(0, t):
    df["re2"][i] = df["right_edge"][i]


df["to"] = [[None] * 50] * t
df["to_length"] = [[None] * 50] * t
for i in range(0, t - 1):
    df["to"][i] = df["center"][i + 1]
    df["to_length"][i] = df["len_mat"][i + 1]


for i in range(0, t):
    count2 = 0
    for j in range(len(df["right_edge"][i])):
        if df["assign"][i][j] == "divide":
            df["from"][i] = np.insert(
                df["from"][i], j + count2, df["from"][i][j + count2]
            )
            df["re2"][i] = np.insert(df["re2"][i], j + count2, df["re2"][i][j + count2])

            count2 += 1


"""
df['to2'] = df['to']
for i in range(0, t-1):
    for j in range(len(df['to'][i])):
        df['to2'][i][j] = np.asarray(df['to'][i][j]) + (0.5 * np.asarray(df['to_length'][i][j] ))

df['from2'] = df['from']
for i in range(0, t-1):
    for j in range(len(df['from'][i])):
        df['from2'][i][j] = np.asarray(df['from'][i][j]) + (0.5 * np.asarray(df['to_length'][i][j] ))
"""


m1 = []
for i in range(t):
    # img1 = Image.new("RGB", (2740, 150), "rgba(18,18,18,100)")
    img1 = Image.new("RGB", (1440, 150), "rgba(18,18,18,100)")
    m1 = m1_mat[i]
    img = position_ellipse(img1, m1, pos_mat_gen(m1), [22] * len(m1), i)

import glob

import cv2

images = []

filenames = [img for img in glob.glob("*.jpg")]
filenames.sort()

int_files = [int(filenames[i][:-4]) for i in range(len(filenames))]
int_files.sort()
final_files = ["{}.jpg".format((str(i))) for i in int_files]


images = []
for img in final_files:
    n = cv2.imread(img)
    images.append(n)
    # print (img)


img0 = cv2.imread("0.jpg")
height, width, layers = img0.shape

video = cv2.VideoWriter(
    "video_intensity12.mp4", -1, 5, (width, height)
)  # second last is fps should be a factor of total time points

for i in range(0, t):
    video.write(images[i])


cv2.destroyAllWindows()
video.release()

####
concatenated2 = np.concatenate([x for x in images[0 : t - 1]], axis=0)


df.to_csv("df3.csv")
for m in range(0, t - 1):
    for i in range(len(df["to"][m])):
        # print(m,i,len(df['from2'][m]), len(df['to2'][m]))
        # x1= df['from'][m][i]*1.4+135#    +150
        # x2=df['to'][m][i]*1.4 +135#    +150

        if df["re2"][m][i] < 500:
            x1 = df["from"][m][i] * 1.4 + 135  #    +150
            x2 = df["to"][m][i] * 1.4 + 135

        else:
            x1 = df["from"][m][i] * 1.4 + 135  #    +150
            x2 = df["to"][m][i] * 1.4 + 135

        # print(x1,x2)
        concatenated2 = cv2.line(
            concatenated2,
            (int(x1), int(38.5) + height * (m)),
            (int(x2), int(38.5 + (height * (m + 1)))),
            (255, 0, 0),
            5,
        )

concatenated2 = Image.fromarray(concatenated2)
concatenated2.save("mask_img.jpg")

df = df.drop(columns=["from", "to", "to_length", "re2"])
df.to_csv("df2.csv")
