from PIL import Image, ImageDraw


def make_rectangle(ID, height, width):
    return ((ID * width) + "\n") * height


print(
    make_rectangle("0", 2, 122),
    "0111111111100222222222200333333333300444444444400555555555500666666666600777777777700888888888800999999999900aaaaaaaaaa00",
    "00111111111100222222222200333333333300444444444400555555555500666666666600777777777700888888888800999999999900aaaaaaaaaa00",
    "00111111111100222222222200333333333300444444444400555555555500666666666600777777777700888888888800999999999900aaaaaaaaaa00",
    "00111111111100222222222200333333333300444444444400555555555500666666666600777777777700888888888800999999999900aaaaaaaaaa00",
    make_rectangle("0", 2, 122),
)


img = Image.new("RGB", (1525, 100), color="black")

img2 = Image.new("RGB", (125, 50), "white")
draw = ImageDraw.Draw(img2)
r, g, b = 0, 0, 0
dr = (14 - r) / 5.0
dg = (211 - g) / 5.0
db = (32 - b) / 5.0
for i in range(500):
    r, g, b = r + dr, g + dg, b + db
    draw.line((i, 0, i, 500), fill=(int(r), int(g), int(b)))

# img2.show()

img.paste(img2, (25, 25))
img.paste(img2, (175, 25))
img.paste(img2, (325, 25))
img.paste(img2, (325 + 150, 25))
img.paste(img2, (625, 25))
img.paste(img2, (625 + 150, 25))
img.paste(img2, (925, 25))
img.paste(img2, (925 + 150, 25))
img.paste(img2, (1225, 25))
img.paste(img2, (1225 + 150, 25))


from PIL import ImageFilter

im1 = img.filter(ImageFilter.BoxBlur(90))
img.show()


from random import randint as rint

"""
img = Image.new('RGB', (1525, 100), color = 'black')

img2 = Image.new('RGB', (125,50), color = 'white')


img.paste(img2, (25, 25)) 
img.paste(img2, (175, 25)) 
img.paste(img2, (325, 25)) 
img.paste(img2, (325+150, 25)) 
img.paste(img2, (625, 25)) 
img.paste(img2, (625+150, 25)) 
img.paste(img2, (925, 25)) 
img.paste(img2, (925+150, 25)) 
img.paste(img2, (1225, 25)) 
img.paste(img2, (1225+150, 25)) 

#img.show()


from PIL import ImageFilter

im1 = img.filter(ImageFilter.BoxBlur(9))
img.show()
"""
