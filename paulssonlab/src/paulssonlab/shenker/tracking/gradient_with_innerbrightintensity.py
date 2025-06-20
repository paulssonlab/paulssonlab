import math

import cv2
import numpy as np
from PIL import Image, ImageDraw, ImageFilter

# ============================================================================


def ellipse_bbox(h, k, a, b, theta):
    ux = a * math.cos(theta)
    uy = a * math.sin(theta)
    vx = b * math.cos(theta + math.pi / 2)
    vy = b * math.sin(theta + math.pi / 2)
    box_halfwidth = np.ceil(math.sqrt(ux**2 + vx**2))
    box_halfheight = np.ceil(math.sqrt(uy**2 + vy**2))
    return (
        (int(h - box_halfwidth), int(k - box_halfheight)),
        (int(h + box_halfwidth), int(k + box_halfheight)),
    )


# ----------------------------------------------------------------------------


# Rotated elliptical gradient - faster, vectorized numpy approach
def make_gradient_v2(width, height, h, k, a, b, theta):
    # Precalculate constants
    st, ct = math.sin(theta), math.cos(theta)
    aa, bb = a**2, b**2

    # Generate (x,y) coordinate arrays
    y, x = np.mgrid[-k : height - k, -h : width - h]
    # Calculate the weight for each pixel
    weights = (((x * ct + y * st) ** 2) / aa) + (((x * st - y * ct) ** 2) / bb)

    return np.clip(1.0 - weights, 0, 1)


# ============================================================================


def draw_image(a, b, theta, inner_scale, save_intermediate=False):
    # Calculate the image size needed to draw this and center the ellipse
    _, (h, k) = ellipse_bbox(0, 0, a, b, theta)  # Ellipse center
    h += 2  # Add small margin
    k += 2  # Add small margin
    width, height = (h * 2 + 1, k * 2 + 1)  # Canvas size

    # Parameters defining the two ellipses for OpenCV (a RotatedRect structure)
    ellipse_outer = ((h, k), (a * 2, b * 2), math.degrees(theta))
    ellipse_inner = (
        (h, k),
        (a * 2 * inner_scale, b * 2 * inner_scale),
        math.degrees(theta),
    )

    # Generate the transparency layer -- the outer ellipse filled and anti-aliased
    transparency = np.zeros((height, width), np.uint8)
    cv2.ellipse(transparency, ellipse_outer, 255, -1, cv2.LINE_AA)

    # Generate the gradient and scale it to 8bit grayscale range
    intensity = np.uint8(make_gradient_v2(width, height, h, k, a, b, theta) * 255)

    # Draw the inter ellipse filled and anti-aliased
    cv2.ellipse(intensity, ellipse_inner, 255, -1, cv2.LINE_AA)

    # Turn it into a BGRA image
    result = cv2.merge([intensity, intensity, intensity, transparency])
    return result


# ============================================================================

"""a, b = (100.0, 20.0) # Semi-major and semi-minor axis
theta = math.radians(40.0) # Ellipse rotation (radians)
inner_scale = 0.6 # Scale of the inner full-white ellipse"""

# cv2.imwrite("eligrad2.png", draw_image(a, b, theta, inner_scale, True))


img1 = Image.new("RGB", (1220, 420), "grey")


img = Image.new("RGB", (1200, 80), "black")


cv2.imwrite("eligrad2a.png", draw_image(80, 24, -0.1, 0.6, True))
ima = Image.open("eligrad2a.png")
ima = ima.filter(ImageFilter.GaussianBlur(3))
img.paste(ima, (0, 10), ima)

cv2.imwrite("eligrad2b.png", draw_image(80, 24, 0.1, 0.6, True))
imb = Image.open("eligrad2b.png")
imb = imb.filter(ImageFilter.GaussianBlur(3))
img.paste(imb, (140, 10), imb)

cv2.imwrite("eligrad2c.png", draw_image(140, 24, 0, 0.6, True))
imc = Image.open("eligrad2c.png")
imc = imc.filter(ImageFilter.GaussianBlur(3))
img.paste(imc, (300, 10), imc)

cv2.imwrite("eligrad2d.png", draw_image(140, 24, -0.07, 0.5, True))
imd = Image.open("eligrad2d.png")
imd = imd.filter(ImageFilter.GaussianBlur(3))
img.paste(imd, (545, 10), imd)

cv2.imwrite("eligrad2e.png", draw_image(140, 24, 0.02, 0.65, True))
ime = Image.open("eligrad2e.png")
ime = ime.filter(ImageFilter.GaussianBlur(3))
img.paste(ime, (805, 10), ime)


cv2.imwrite("eligrad2b.png", draw_image(80, 24, 0.1, 0.65, True))
imb = Image.open("eligrad2b.png")
imb = imb.filter(ImageFilter.GaussianBlur(3))
img.paste(imb, (1100, 10), imb)

img1.paste(img, (10, 10))


img = Image.new("RGB", (1200, 80), "black")


cv2.imwrite("eligrad2c.png", draw_image(140, 24, 0, 0.65, True))
imc = Image.open("eligrad2c.png")
imc = imc.filter(ImageFilter.GaussianBlur(3))
img.paste(imc, (0, 10), imc)

cv2.imwrite("eligrad2c.png", draw_image(140, 24, 0, 0.65, True))
imc = Image.open("eligrad2c.png")
imc = imc.filter(ImageFilter.GaussianBlur(3))
img.paste(imc, (240, 10), imc)

cv2.imwrite("eligrad2c.png", draw_image(140, 24, 0, 0.65, True))
imc = Image.open("eligrad2c.png")
imc = imc.filter(ImageFilter.GaussianBlur(3))
img.paste(imc, (480, 10), imc)

cv2.imwrite("eligrad2c.png", draw_image(140, 24, 0, 0.65, True))
imc = Image.open("eligrad2c.png")
imc = imc.filter(ImageFilter.GaussianBlur(3))
img.paste(imc, (480 + 240, 10), imc)

cv2.imwrite("eligrad2c.png", draw_image(140, 24, 0, 0.65, True))
imc = Image.open("eligrad2c.png")
imc = imc.filter(ImageFilter.GaussianBlur(3))
img.paste(imc, (480 + 240 + 240, 10), imc)

img1.paste(img, (10, 100))


img1.show()
