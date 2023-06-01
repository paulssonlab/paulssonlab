import random

import numpy as np
import skimage
from PIL import Image, ImageDraw, ImageFilter

# img1 = Image.new("RGB",(1740,150),'rgba(18,18,18,100)')
# img2 = Image.new("RGBA",(800,80),'rgba(18,18,18,100)')
# img2 = img2.filter(ImageFilter.GaussianBlur(3))

# img1.paste(img2, (5,15))
angle = 0.1


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


def divide(length, x_there, x_here):
    length = int(length / 2)
    img = Image.new(
        "RGBA", (800, 160), "rgba(0,0,0,0)"
    )  # we do get a background less image
    ell0 = make_ellipse(img, x - 12, y - 12, dx + 12, dy + 12, length, bright_halo)
    ell1 = make_ellipse(img, x, y, dx, dy, length, bright_blur)
    ell1 = make_ellipse(img, x + 7, y + 5, dx - 7, dy - 5, length, bright)

    # img=img.rotate(angle)
    img = img.rotate(random.choice(list2))
    # result = img
    imgSmall = img.resize((pixelsize, pixelsize), resample=Image.BILINEAR)
    result = imgSmall.resize(img.size, Image.NEAREST)
    result = result.filter(ImageFilter.GaussianBlur(g_blur_r))

    img1.paste(result, (x_there, 0), result)

    img = Image.new(
        "RGBA", (800, 160), "rgba(0,0,0,0)"
    )  # we do get a background less image
    ell0 = make_ellipse(img, x - 12, y - 12, dx + 12, dy + 12, length, bright_halo)
    ell1 = make_ellipse(img, x, y, dx, dy, length, bright_blur)
    ell1 = make_ellipse(img, x + 7, y + 5, dx - 7, dy - 5, length, bright)

    # img=img.rotate(angle)
    img = img.rotate(random.choice(list2))
    # result = img
    imgSmall = img.resize((pixelsize, pixelsize), resample=Image.BILINEAR)
    result = imgSmall.resize(img.size, Image.NEAREST)
    result = result.filter(ImageFilter.GaussianBlur(g_blur_r))

    img1.paste(result, (x_here, 0), result)
    return img1


bright_halo = "rgb(68,68,68)"
bright_blur = "rgb(128,128,128)"
bright = "rgb(225,225,225)"
pixelsize = 50
g_blur_r = 2
t = 25
list2 = [2, 3, 4, 8]


for i in range(0, t):
    img1 = Image.new("RGB", (2740, 150), "rgba(18,18,18,100)")

    angle = 0.1
    print(i)
    # img = Image.new("RGBA",(800,160),'rgba(0,0,0,0)') # we do get a background less image
    length = 50 + i * 10  # in
    x = 125
    y = 22
    dx = x + 75
    dy = y + 33

    if length < 100:
        img = Image.new(
            "RGBA", (800, 160), "rgba(0,0,0,0)"
        )  # we do get a background less image
        ell0 = make_ellipse(img, x - 12, y - 12, dx + 12, dy + 12, length, bright_halo)
        ell1 = make_ellipse(img, x, y, dx, dy, length, bright_blur)
        ell1 = make_ellipse(img, x + 7, y + 5, dx - 7, dy - 5, length, bright)

        img = img.rotate(random.choice(list2))
        imgSmall = img.resize((pixelsize, pixelsize), resample=Image.BILINEAR)
        result = imgSmall.resize(img.size, Image.NEAREST)
        result = result.filter(ImageFilter.GaussianBlur(g_blur_r))

        img1.paste(result, (0, 0), result)

        img1 = img1.save("%d.png" % i)

    elif length > 90:
        length = int(length / 2)
        length1 = length

        img = Image.new(
            "RGBA", (800, 160), "rgba(0,0,0,0)"
        )  # we do get a background less image
        ell0 = make_ellipse(img, x - 12, y - 12, dx + 12, dy + 12, length, bright_halo)
        ell1 = make_ellipse(img, x, y, dx, dy, length, bright_blur)
        ell1 = make_ellipse(img, x + 7, y + 5, dx - 7, dy - 5, length, bright)

        img = img.rotate(random.choice(list2))
        imgSmall = img.resize((pixelsize, pixelsize), resample=Image.BILINEAR)
        result = imgSmall.resize(img.size, Image.NEAREST)
        result = result.filter(ImageFilter.GaussianBlur(g_blur_r))

        if length < 100:
            img1.paste(result, (0, 0), result)
        ##

        elif length > 90:
            img1 = divide(length, 0, int(99 + (i * 2)))

        length = length1
        img = Image.new(
            "RGBA", (800, 160), "rgba(0,0,0,0)"
        )  # we do get a background less image
        ell0 = make_ellipse(img, x - 12, y - 12, dx + 12, dy + 12, length, bright_halo)
        ell1 = make_ellipse(img, x, y, dx, dy, length, bright_blur)
        ell1 = make_ellipse(img, x + 7, y + 5, dx - 7, dy - 5, length, bright)

        img = img.rotate(random.choice(list2))
        imgSmall = img.resize((pixelsize, pixelsize), resample=Image.BILINEAR)
        result = imgSmall.resize(img.size, Image.NEAREST)
        result = result.filter(ImageFilter.GaussianBlur(g_blur_r))

        if length < 100:
            img1.paste(result, (int(99 + (i * 12)), 0), result)

        elif length > 90:
            img1 = divide(
                length, int(99 + (i * 12)), (int(99 + (i * 12))) + (int(99 + (i * 2)))
            )

        img1 = img1.save("%d.png" % i)


import cv2

# import cv2
img0 = cv2.imread("0.png")
img1 = cv2.imread("1.png")
img2 = cv2.imread("2.png")
img3 = cv2.imread("3.png")
img4 = cv2.imread("4.png")
img5 = cv2.imread("5.png")
img6 = cv2.imread("6.png")
img7 = cv2.imread("7.png")
img8 = cv2.imread("8.png")
img9 = cv2.imread("9.png")
img10 = cv2.imread("10.png")
img11 = cv2.imread("11.png")
img12 = cv2.imread("12.png")
img13 = cv2.imread("13.png")
img14 = cv2.imread("14.png")
img15 = cv2.imread("15.png")
img16 = cv2.imread("16.png")
img17 = cv2.imread("17.png")
img18 = cv2.imread("18.png")
img19 = cv2.imread("19.png")
img20 = cv2.imread("20.png")
img21 = cv2.imread("21.png")
img22 = cv2.imread("22.png")
img23 = cv2.imread("23.png")
img24 = cv2.imread("24.png")

height, width, layers = img1.shape

video = cv2.VideoWriter("video5o.mp4", -1, 5, (width, height))

video.write(img0)
video.write(img1)
video.write(img2)
video.write(img3)
video.write(img4)
video.write(img5)
video.write(img6)
video.write(img7)
video.write(img8)
video.write(img9)
video.write(img10)
video.write(img11)
video.write(img12)
video.write(img13)
video.write(img14)
video.write(img15)
video.write(img16)
video.write(img17)
video.write(img18)
video.write(img19)
video.write(img20)
video.write(img21)
video.write(img22)
video.write(img23)
video.write(img24)

cv2.destroyAllWindows()
video.release()
