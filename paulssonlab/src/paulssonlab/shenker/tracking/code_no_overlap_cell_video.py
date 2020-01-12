from PIL import Image, ImageDraw, ImageFilter
import numpy as np
import skimage
import random
import math


bright_halo = "rgb(68,68,68)"
bright_blur = "rgb(128,128,128)"
bright = "rgb(225,225,225)"
pixelsize = 50
g_blur_r = 2
t = 150
list2 = [2, 3, 4, 5, 0, 1]
# list2=[0]
pixelsize_list = [50, 80, 30, 40]


# bright_halo_list = ['rgb(68,68,68)','rgb(48,48,48)','rgb(88,88,88)']
# bright_blur_list= ['rgb(128,128,128)','rgb(148,148,148)','rgb(108,108,108)']
# bright_list = ['rgb(225,225,225)','rgb(240,240,240)','rgb(200,200,200)']

bright_halo_list = ["rgb(225,225,225)"]
bright_blur_list = ["rgb(225,225,225)"]
# bright_list = ['rgb(225,225,0)','rgb(225,0,225)','rgb(225,225,225)']
bright_list = ["rgb(225,225,225)"]


def make_ellipse(image, x, y, dx, dy, next_1, color):

    dr = ImageDraw.Draw(image)
    ell = dr.ellipse((x, y, dx, dy), fill=color)
    ell2 = dr.ellipse((x + next_1, y, dx + next_1, dy), color)
    img1 = Image.new("RGB", (next_1, dy - y + 1), color=color)
    image.paste(img1, ((x + int((dx - x) / 2), int(y))))
    return ell, ell2, img


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


def divide(length, x_there, x_here):  # x_here is x_there + divide_length

    length = int(length / 2)

    # x_there = int(x_there)
    # x_here = x_here+int(75*int(i/8))#/40))

    img = draw_ellipse(length)
    result = image_edit(img, random.choice(list2))
    if length < 100:
        img1.paste(result, (x_there, 0), result)
        img1.paste(result, (x_here, 0), result)  # +(2*i)#i%18)*4
        return img1
    if length > 90:
        # print('i',i)
        m = 0

        # if (i>78):
        #   m = (i/75)
        # x_here=x_here+(int(75+length)*2)#(int((75+(length/2))*(i/35)))
        # x_here_add =1000# int(75*int(i/40)) #between the 2 cell sets
        img2 = divide(length, x_there, x_there + (int(75 + (length / 2))))
        img2 = divide(
            length,
            x_here + (int((75 + (length / 2)) * (0.5 + m))),
            x_here + (int((75 + (length / 2)) * (1.5 + m))),
        )

        # img2 = divide(length,x_there, x_there+int(divide_length))
        # img2 = divide(length,x_here+(2*(i-int(35*(i/40))))+4*i, x_here+int(divide_length)+(2*(i-int(35*(i/40))))+4*i)
        # img2 = divide(length,x_here+(2*(i))+4*i, x_here+divide_length+(2*(i))+4*i)
        ##img2.show()
        return img2


# (i-(40*int(i/40))) = i


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
    divide_length = int(99 + (2 * 2))

    if length < 100:

        img = draw_ellipse(length)
        result = image_edit(img, random.choice(list2))
        img1.paste(result, (0, 0), result)
        img1 = image_edit2(img1)  ###
        img1 = img1.save("%d.png" % i)

    elif length > 90:

        length = int(length / 2)
        length1 = length

        img = draw_ellipse(length)
        result = image_edit(img, random.choice(list2))

        if length < 100:
            img1.paste(result, (0, 0), result)

        elif length > 90:
            img1 = divide(length, 0, (int(75 + (length / 2))))

        # i = i - (int(i/10)*10)

        length = length1

        # if i>15:
        #    length = length+ int(i/15)*int(length * (i%15)/15)

        img = draw_ellipse(length)
        result = image_edit(img, random.choice(list2))

        if length < 100:
            img1.paste(result, (int(75 + (length)), 0), result)  # +((i-1)*5)

        elif length > 90:
            print(i, "i")

            if i < 70:
                m = i / 35
            elif i > 69 and i < 135:
                m = i / 69
                m = int(m)
            elif i > 134:
                m = i / 134
                m = int(m)
            # img1 = divide(length,int((75*2+length))+ (int(75+(length/2))), int((75*2+length))+ (int(75+(length/2))*2))#######
            img1 = divide(
                length,
                int((75 * 2 + length) + ((int(75 + (length / 2))) * int(m))),
                int((75 * 2 + length)) + ((int(75 + (length / 2))) * (1 + int(m))),
            )  #######

        img1 = image_edit2(img1)  ###
        img1 = img1.save("%d.png" % i)


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
    # print (img)


img0 = cv2.imread("0.png")
height, width, layers = img0.shape

video = cv2.VideoWriter("video5zzzzzzzzzzzj.mp4", -1, 5, (width, height))

for i in range(0, t):

    video.write(images[i])


cv2.destroyAllWindows()
video.release()
