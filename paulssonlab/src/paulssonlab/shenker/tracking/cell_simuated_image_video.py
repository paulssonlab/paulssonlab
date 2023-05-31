import math
import random

import numpy as np
import skimage
from PIL import Image, ImageDraw, ImageFilter

bright_halo = "rgb(68,68,68)"
bright_blur = "rgb(128,128,128)"
bright = "rgb(225,225,225)"
pixelsize = 50
g_blur_r = 2
t = 75  # 150
# list2=[2,0,1]#,5,8,10]
list2 = [0]
pixelsize_list = [50, 80, 30, 40]


bright_halo_list = ["rgb(68,68,68)", "rgb(48,48,48)", "rgb(88,88,88)"]
# bright_blur_list= ['rgb(128,128,128)','rgb(148,148,148)','rgb(108,108,108)']
# bright_list = ['rgb(225,225,225)','rgb(250,250,250)','rgb(180,180,180)']

# bright_halo_list = ['rgb(225,225,225)']
bright_blur_list = ["rgb(30,30,30)"]
# bright_list = ['rgb(225,225,0)','rgb(225,0,225)','rgb(225,225,225)']
bright_list = ["rgb(225,225,225)"]


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
    # ell0 = make_ellipse(img,x-12,y-12,dx+12,dy+12,length,bright_halo)
    ell1 = make_ellipse(img, x, y, dx, dy, length, bright_blur)
    ell1 = make_ellipse(img, x + 7, y + 5, dx - 7, dy - 5, length, bright)

    return img


def divide(length, x_there, x_here, i):
    length = int(length / 2)
    length1 = length

    oo = 0

    if length < 80:
        oo = -int((75 + (length)) / 4)  # -int((length/2))########2
    else:
        oo = -int((75 + (length)) / 8)  ######10

    """x_lol=0
            if i>74:
                    m = -0.25+0.25-0.2
                    x_lol=0"""

    if length < 100:
        length = int(length * 2.3)
        img = draw_ellipse(length)
        result = image_edit(img, random.choice(list2))
        img1.paste(result, (x_there, 0), result)
        img = draw_ellipse(length)
        result = image_edit(img, random.choice(list2))
        img1.paste(result, (x_here + (oo) + int(length * 0.2), 0), result)
        return img1

    if length > 90:
        # length = int(length*2)#####################
        # length = length1
        m = -0.25  #########################################

        oo1 = 0
        if i > 74:
            oo1 = 0
            m = -0.25 + 0.25 - 0.2 + 0.5
            oo1 == oo * 2

        x_here = x_here + oo
        img2 = divide(
            length1, x_there, x_there + (int(75 + (length / 2))) - oo1, i
        )  #############oo
        img2 = divide(
            length1,
            x_here + (int((75 + (length / 2)) * (0.5 + m))),
            x_here + (int((75 + (length / 2)) * (1.5 + m))),
            i,
        )

        return img2


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
    result = result.filter(ImageFilter.GaussianBlur(2))

    return result


for i in range(0, t):
    img1 = Image.new("RGB", (2740, 150), "rgba(18,18,18,100)")

    length = 50 + i * 10  # in
    x = 125
    y = 22
    dx = x + 75
    dy = y + 33
    # divide_length = int(99+(2*2))
    if length < 80:
        img = draw_ellipse(int(length * 1.5))
        result = image_edit(img, random.choice(list2))
        img1.paste(result, (0, 0), result)
        img1 = image_edit2(img1)  ###

        """img1 = img1.convert('LA')
        img1_ary = np.asarray(img1)[:,:,0]
        mask = img1_ary > skimage.filters.threshold_otsu(img1_ary)
        Image.fromarray(mask*255.).show()"""

        img1 = img1.save("%d.jpg" % i)

    if length < 100 and length > 70:
        img = draw_ellipse(int(length * 2.6))
        result = image_edit(img, random.choice(list2))
        img1.paste(result, (0, 0), result)
        img1 = image_edit2(img1)  ###

        """img1 = img1.convert('LA')
        img1_ary = np.asarray(img1)[:,:,0]
        mask = img1_ary > skimage.filters.threshold_otsu(img1_ary)
        Image.fromarray(mask*255.).show()"""

        img1 = img1.save("%d.jpg" % i)

    elif length > 90:
        length = int(length / 2)
        length1 = length

        img = draw_ellipse(int(length * 2.6))
        result = image_edit(img, random.choice(list2))

        if length < 100:
            img1.paste(result, (0, 0), result)

        elif length > 90:
            m = 1
            if i > 34 and i < 76 - 1:
                m = 1
            if i > 75 - 1:
                m = 1.5 - 0.25 - 0.5
            img1 = divide(length, 0, int(((75 + int(length * 1.9)) / 2) * m), i)

        length = length1

        oo = 0

        if length < 80:
            print("m", i, length)
            oo = -int((75 + (length)) * (0.5 + ((6 - i) / 6)))
        else:
            oo = 0  # -int((75+(length))/2)
            print(i, "j")

        img = draw_ellipse(int(length * 2.6))
        result = image_edit(img, random.choice(list2))

        if length < 100:
            img1.paste(result, (int(75 + int(length * 2.4)) + oo, 0), result)

        elif length > 90:
            print(i, "i")

            length = int(length * 1.6)
            m = 0
            n = 0

            if i == 15:
                m = -0.5
            if i < 35 and i > 15:
                m = -0.5
            if i > 34 and i < 76 - 1:
                m = 0.75 - 0.5 - 0.1 - 0.5 - 0.2  ###########################

            if i > 75 - 1:
                m = (
                    1.5 - 0.25 - 0.3 - 0.1 - 0.5
                )  ######------------------------------0.5
                n = 0  # 0.25

            img1 = divide(
                length1,
                int((75 * 2 + length) + int((int(75 + (length / 2))) * m))
                + oo
                + int(length * 0.4),
                int(75 * 2 + length)
                + int((int(75 + (length / 2))) * (1 + m + n))
                + oo
                + int(length * 0.4),
                i,
            )  #######
        img1 = image_edit2(img1)
        img1 = img1.save("%d.jpg" % i)
        """
        img1 = img1.convert('LA')
        img1_ary = np.asarray(img1)[:,:,0]
        mask = img1_ary > skimage.filters.threshold_otsu(img1_ary)
        #Image.fromarray(mask*255.).show()"""


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


img0 = cv2.imread("20.jpg")
height, width, layers = img0.shape

video = cv2.VideoWriter(
    "video_test_5.mp4", -1, 2, (width, height)
)  # second last is fps should be a factor of total time points

for i in range(0, t):
    video.write(images[i])


cv2.destroyAllWindows()
video.release()
