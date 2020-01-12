from PIL import Image, ImageDraw, ImageFilter
import numpy as np
import skimage
import random


bright_halo = "rgb(68,68,68)"
bright_blur = "rgb(128,128,128)"
bright = "rgb(225,225,225)"
pixelsize = 50
g_blur_r = 2
t = 150
list2 = [2, 3, 4, 8]


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

    img = Image.new(
        "RGBA", (800, 160), "rgba(0,0,0,0)"
    )  # we do get a background less image
    ell0 = make_ellipse(img, x - 12, y - 12, dx + 12, dy + 12, length, bright_halo)
    ell1 = make_ellipse(img, x, y, dx, dy, length, bright_blur)
    ell1 = make_ellipse(img, x + 7, y + 5, dx - 7, dy - 5, length, bright)

    return img


def divide(length, x_there, x_here):  # x_here is x_there + divide_length

    length = int(length / 2)

    img = draw_ellipse(length)
    result = image_edit(img)
    if length < 100:
        img1.paste(result, (x_there, 0), result)
        img1.paste(result, (x_here + (2 * i), 0), result)  # +(2*i)
        return img1
    if length > 90:

        img2 = divide(length, x_there, x_there + divide_length)
        img2 = divide(length, x_here + (4 * i), x_here + divide_length + (4 * i))
        # img2.show()
        return img2


def image_edit(img):

    img = img.rotate(random.choice(list2))
    imgSmall = img.resize((pixelsize, pixelsize), resample=Image.BILINEAR)
    result = imgSmall.resize(img.size, Image.NEAREST)
    result = result.filter(ImageFilter.GaussianBlur(g_blur_r))

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
        result = image_edit(img)
        img1.paste(result, (0, 0), result)
        img1 = img1.save("%d.png" % i)

    elif length > 90:

        length = int(length / 2)
        length1 = length

        img = draw_ellipse(length)
        result = image_edit(img)

        if length < 100:
            img1.paste(result, (0, 0), result)

        elif length > 90:
            img1 = divide(length, 0, divide_length)

        length = length1

        img = draw_ellipse(length)
        result = image_edit(img)

        if length < 100:
            img1.paste(result, (int(99 + (i * 5)), 0), result)  #

        elif length > 90:
            img1 = divide(
                length, int(99 + (i * 10)), (int(99 + (i * 10))) + divide_length
            )  #

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
    print(img)


img0 = cv2.imread("0.png")
height, width, layers = img0.shape

video = cv2.VideoWriter("video5zx.mp4", -1, 10, (width, height))

for i in range(0, t):

    video.write(images[i])


cv2.destroyAllWindows()
video.release()
