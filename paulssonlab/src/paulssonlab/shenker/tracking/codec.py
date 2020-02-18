from PIL import Image, ImageDraw, ImageFilter


img1 = Image.new("RGB", (840, 400), "grey")
img2 = Image.new("RGBA", (800, 80), "rgba(18,18,18,100)")
img2 = img2.filter(ImageFilter.GaussianBlur(3))

img1.paste(img2, (5, 15))


def make_ellipse(image, x, y, dx, dy, next_1, color):

    dr = ImageDraw.Draw(image)

    # ell3 = dr.ellipse((x-10,y-10,dx+next_1+10,dy+10),fill ='red')
    ell = dr.ellipse((x, y, dx, dy), fill=color)
    ell2 = dr.ellipse((x + next_1, y, dx + next_1, dy), color)

    img1 = Image.new("RGB", (next_1, dy - y + 1), color=color)

    image.paste(img1, ((x + int((dx - x) / 2), int(y))))
    #
    return ell, ell2, img  # ,ell3


def halo(image, x, y, dx, dy, color):

    dr = ImageDraw.Draw(image)
    ell = dr.ellipse((x, y, dx, dy), fill=color)

    return ell


img = Image.new(
    "RGBA", (800, 160), "rgba(0,0,0,0)"
)  # we do get a background less image

ell0 = make_ellipse(img, 0, 0, 84 + 20, 55, 50, "rgb(38,38,38)")
ell1 = make_ellipse(img, 5 + 10, 2 + 10, 79 + 10, 33 + 10, 50, "rgb(128,128,128)")
ell1 = make_ellipse(img, 17 + 10, 7 + 10, 66 + 10, 28 + 10, 50, "rgb(240,240,240)")
# img = img.filter(ImageFilter.GaussianBlur(3))
img = img.rotate(3.23)
imga = img.save("did_id.png")
imgc = Image.open("did_id.png")

imgSmall = imgc.resize((50, 50), resample=Image.BILINEAR)

result = imgSmall.resize(imgc.size, Image.NEAREST)
result = result.filter(ImageFilter.GaussianBlur(2))
# imc = im.filter(ImageFilter.BoxBlur(5))
img1.paste(result, (8, 3), result)


img = Image.new(
    "RGBA", (800, 160), "rgba(0,0,0,0)"
)  # we do get a background less image
ell0 = make_ellipse(img, 0 + 140, 0, 84 + 20 + 140, 55, 50, "rgb(28,28,28)")
ell1 = make_ellipse(img, 5 + 140, 2 + 13, 79 + 140, 33 + 13, 50, "rgb(128,128,128)")
ell1 = make_ellipse(img, 17 + 140, 7 + 13, 66 + 140, 28 + 13, 50, "rgb(180,180,180)")
img = img.rotate(3.23)
imga = img.save("did_id.png")
imgc = Image.open("did_id.png")

imgSmall = imgc.resize((50, 50), resample=Image.BILINEAR)

result = imgSmall.resize(imgc.size, Image.NEAREST)
result = result.filter(ImageFilter.GaussianBlur(2))
img1.paste(result, (0, 10), result)


img = Image.new(
    "RGBA", (800, 160), "rgba(0,0,0,0)"
)  # we do get a background less image
ell0 = make_ellipse(img, 0 + 140, 0, 84 + 20 + 140, 55, 50, "rgb(28,28,28)")
ell1 = make_ellipse(img, 5 + 140, 2 + 8, 79 + 140, 33 + 8, 50, "rgb(128,128,128)")
ell1 = make_ellipse(img, 17 + 140, 7 + 8, 66 + 140, 28 + 8, 50, "rgb(220,220,220)")

img = img.rotate(7)
imga = img.save("did_id.png")
imgc = Image.open("did_id.png")

imgSmall = imgc.resize((50, 50), resample=Image.BILINEAR)

result = imgSmall.resize(imgc.size, Image.NEAREST)
result = result.filter(ImageFilter.GaussianBlur(3))
img1.paste(result, (120, 3), result)


img = Image.new(
    "RGBA", (800, 160), "rgba(0,0,0,0)"
)  # we do get a background less image
ell0 = make_ellipse(img, 0 + 140, 0 + 25, 84 + 20 + 140, 55 + 25, 50, "rgb(28,28,28)")
ell1 = make_ellipse(img, 5 + 140, 2 + 35, 79 + 140, 33 + 35, 50, "rgb(128,128,128)")
ell1 = make_ellipse(img, 17 + 140, 7 + 35, 66 + 140, 28 + 35, 50, "rgb(220,220,220)")

img = img.rotate(-6)
imga = img.save("did_id.png")
imgc = Image.open("did_id.png")

imgSmall = imgc.resize((50, 50), resample=Image.BILINEAR)

result = imgSmall.resize(imgc.size, Image.NEAREST)
result = result.filter(ImageFilter.GaussianBlur(3))
img1.paste(result, (240, 25), result)


"""

img2 = Image.new("RGB",(800,60),'black')
ell1 = make_ellipse(img2,5,10,60,48,50)
ell2 = make_ellipse(img2,120,10,170,48,50)
ell2 = make_ellipse(img2,280,10,340,48,80)
ell2 = make_ellipse(img2,400,10,430,48,40)
ell2 = make_ellipse(img2,450,10,500,48,40)
ell2 = make_ellipse(img2,560,10,610,48,40)
ell2 = make_ellipse(img2,647,10,700,48,40)

img2=img2.rotate(1)
img1.paste(img2, (20, 80))




img3 = Image.new("RGB",(260,60),'black')

ell2 = make_ellipse(img3,2,10,52,48,40)
ell2 = make_ellipse(img3,70,10,120,48,40)
ell2 = make_ellipse(img3,150,10,210,48,40)

img3 = img3.filter(ImageFilter.GaussianBlur(2))

img3=img3.rotate(2)
img1.paste(img3, (20, 160))



img3 = Image.new("RGB",(800-260,60),'black')

ell2 = make_ellipse(img3,2,10,52,48,40)
ell2 = make_ellipse(img3,70,10,120,48,40)
ell2 = make_ellipse(img3,180,10,250,48,40)
ell4 = make_ellipse(img3,380,10,450,48,40)

img3 = img3.filter(ImageFilter.GaussianBlur(2))

img3=img3.rotate(-2)
img1.paste(img3, (20+260, 160))


"""
# img1= img1.filter(ImageFilter.Blur(1))
img1.show()
