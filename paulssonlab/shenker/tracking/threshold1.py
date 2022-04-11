import cv2
import numpy as np
from matplotlib import pyplot as plt

img = cv2.imread("0.png", 0)
img = cv2.medianBlur(img, 5)

ret, th1 = cv2.threshold(img, 0, 255, cv2.THRESH_BINARY)


titles = ["Original Image", "Global Thresholding (v = 127)"]
images = [img, th1]

"""for i in range(0,2):
    plt.subplot(1,2,i+1),plt.imshow(images[i],'gray')
    plt.title(titles[i])
    plt.xticks([]),plt.yticks([])
plt.show()"""

# th1.show()
Image.fromarray(th1).show()
Image.fromarray(th1).save("1.png")

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

video = cv2.VideoWriter(
    "video_th2.mp4", -1, 2, (width, height)
)  # second last is fps should be a factor of total time points

for i in range(0, 23):

    video.write(images[i])


cv2.destroyAllWindows()
video.release()
