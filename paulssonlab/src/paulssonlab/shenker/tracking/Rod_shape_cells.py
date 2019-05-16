from PIL import Image, ImageDraw, ImageFilter


img1 = Image.new("RGB", (840, 1400), "grey")


def make_ellipse(image, x, y, dx, dy, next_1):

    dr = ImageDraw.Draw(image)
    ell = dr.ellipse((x, y, dx, dy), fill="white")
    ell2 = dr.ellipse((x + next_1, y, dx + next_1, dy), "white")

    img1 = Image.new("RGB", (next_1, dy - y + 1), color="white")

    image.paste(img1, ((x + int((dx - x) / 2), int(y))))
    #
    return ell, ell2, img


img = Image.new("RGB", (800, 60), "black")
ell1 = make_ellipse(img, 5, 10, 60, 48, 50)
ell2 = make_ellipse(img, 120, 10, 170, 48, 50)
ell2 = make_ellipse(img, 230, 10, 290, 48, 80)
ell2 = make_ellipse(img, 380, 10, 430, 48, 40)
ell2 = make_ellipse(img, 450, 10, 500, 48, 40)
ell2 = make_ellipse(img, 560, 10, 610, 48, 40)
ell2 = make_ellipse(img, 647, 10, 700, 48, 40)

img = img.rotate(0)
img1.paste(img, (20, 10))


img2 = Image.new("RGB", (800, 60), "black")
ell1 = make_ellipse(img2, 5, 10, 60, 48, 50)
ell2 = make_ellipse(img2, 120, 10, 170, 48, 50)
ell2 = make_ellipse(img2, 280, 10, 340, 48, 80)
ell2 = make_ellipse(img2, 400, 10, 430, 48, 40)
ell2 = make_ellipse(img2, 450, 10, 500, 48, 40)
ell2 = make_ellipse(img2, 560, 10, 610, 48, 40)
ell2 = make_ellipse(img2, 647, 10, 700, 48, 40)

img2 = img2.rotate(1)
img1.paste(img2, (20, 80))


img3 = Image.new("RGB", (260, 60), "black")

ell2 = make_ellipse(img3, 2, 10, 52, 48, 40)
ell2 = make_ellipse(img3, 70, 10, 120, 48, 40)
ell2 = make_ellipse(img3, 150, 10, 210, 48, 40)

img3 = img3.filter(ImageFilter.GaussianBlur(2))

img3 = img3.rotate(2)
img1.paste(img3, (20, 160))


img3 = Image.new("RGB", (800 - 260, 60), "black")

ell2 = make_ellipse(img3, 2, 10, 52, 48, 40)
ell2 = make_ellipse(img3, 70, 10, 120, 48, 40)
ell2 = make_ellipse(img3, 180, 10, 250, 48, 40)
ell4 = make_ellipse(img3, 380, 10, 450, 48, 40)

img3 = img3.filter(ImageFilter.GaussianBlur(2))

img3 = img3.rotate(-2)
img1.paste(img3, (20 + 260, 160))


# img1= img1.filter(ImageFilter.GaussianBlur(4))
img1.show()
