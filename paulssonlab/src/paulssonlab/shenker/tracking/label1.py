import cv2
import numpy as np
import skimage
from matplotlib import pyplot as plt
from PIL import Image, ImageDraw, ImageFilter

# https://scipy-lectures.org/packages/scikit-image/index.html
from skimage import morphology

"""img = cv2.imread('0.png',0)
img = cv2.medianBlur(img,5)

ret,th1 = cv2.threshold(img,180,255,cv2.THRESH_BINARY)
images = [img, th1]


img1 = skimage.morphology.label(th1, neighbors=None, background=None, return_num=False, connectivity=None)
print(img1) # these are the labels


np.savetxt('file1.txt',img1)

#img2 = np.where(img1>0, th1, th1*0)
img2 = np.where(th1>0, th1,220)  # changing the last parameter changes the background color
img3 = np.where(img1==5, th1,0)   # this prints only the areas which are labeled as img== label

np.savetxt('file0.txt',th1)

#Image.fromarray(th1).show()
##Image.fromarray(img2).show()
Image.fromarray(img3).show()

np.savetxt('file2.txt',img2)"""

"""img = cv2.imread('0.png',0)
img = cv2.medianBlur(img,5)

ret,th1 = cv2.threshold(img,170,255,cv2.THRESH_BINARY)
images = [img, th1]


img1 = skimage.morphology.label(th1, neighbors=None, background=None, return_num=False, connectivity=None)
print(img1) # these are the labels


np.savetxt('file1.txt',img1)

#img2 = np.where(img1>0, th1, th1*0)
img2 = np.where(th1>0, th1,220)  # changing the last parameter changes the background color
img3 = np.where(img1==2, th1,0)   # this prints only the areas which are labeled as img== label

np.savetxt('file0.txt',th1)

#Image.fromarray(th1).show()
##Image.fromarray(img2).show()
Image.fromarray(img3).show()

np.savetxt('file2.txt',img2)"""

####################################################

"""img = cv2.imread('0.png',0)
img = cv2.medianBlur(img,5)

ret,th1 = cv2.threshold(img,170,255,cv2.THRESH_BINARY)
images = [img, th1]

reta,th1a = cv2.threshold(img,180,255,cv2.THRESH_BINARY)
images = [img, th1a]


img1 = skimage.morphology.label(th1, neighbors=None, background=None, return_num=False, connectivity=None)
#print(img1) # these are the labels

img1a = skimage.morphology.label(th1a, neighbors=None, background=None, return_num=False, connectivity=None)
#print(img1a) # these are the labels

np.savetxt('file1.txt',img1)
n_mask_img1 = np.amax(img1)
np.savetxt('file1a.txt',img1a)
n_mask_img1a = np.amax(img1a)


#output = img1[img1==0]
output2 = np.where(img1==2) #[0] gives row indices and [1] gives column indices
mask_1 = [2 for i in range (len(output2[0]))]
listOfCoordinates= list(zip(output2[0], output2[1], mask_1))
#print(output2,listOfCoordinates)
for cord in listOfCoordinates:
    print(cord)

np.savetxt('file_coordinates_listOfCoordinates.txt',listOfCoordinates)


output2a = np.where(img1a==5) #[0] gives row indices and [1] gives column indices
mask_1a = [2 for i in range (len(output2a[0]))]
listOfCoordinatesa= list(zip(output2a[0], output2a[1], mask_1a))
#print(output2a,listOfCoordinatesa)
for cord in listOfCoordinatesa:
    print(cord)

np.savetxt('file_coordinates_listOfCoordinatesa.txt',listOfCoordinatesa)










th1_bg = th1[img1==0]



#np.savetxt('file_coordinates_bg.txt',output2)
#img2 = np.where(img1>0, th1, th1*0)
img2 = np.where(th1>0, th1,220)  # changing the last parameter changes the background color
img3 = np.where(img1a==5, th1,7)   # this prints only the areas which are labeled as img == label #img1a==5 img1==2

np.savetxt('file0.txt',th1)

#Image.fromarray(th1).show()
##Image.fromarray(img2).show()
###Image.fromarray(img3).show()

np.savetxt('file2.txt',img2)"""


img = cv2.imread("0.png", 0)
img = cv2.medianBlur(img, 5)

ret, th1 = cv2.threshold(img, 170, 255, cv2.THRESH_BINARY)
images = [img, th1]


reta, th1a = cv2.threshold(img, 180, 255, cv2.THRESH_BINARY)
imagesa = [img, th1a]


img1 = skimage.morphology.label(
    th1, neighbors=None, background=None, return_num=False, connectivity=None
)
img1a = skimage.morphology.label(
    th1a, neighbors=None, background=None, return_num=False, connectivity=None
)


np.savetxt("file1.txt", img1)
n_mask_img1 = np.amax(img1)
np.savetxt("file1a.txt", img1a)
n_mask_img1a = np.amax(img1a)


"""for i in range (1,n_mask_img1+1):
    output2 = np.where(img1==i)
    mask_1 = [i for j in range (len(output2[0]))]
    listOfCoordinates= list(zip(output2[0], output2[1], mask_1))"""


def get_coordinates(img1, mask_num):
    output2 = np.where(
        img1 == mask_num
    )  # [0] gives row indices and [1] gives column indices
    mask_1 = [mask_num for j in range(len(output2[0]))]
    listOfCoordinates = list(zip(output2[0], output2[1], mask_1))

    return listOfCoordinates


arr1 = get_coordinates(img1, 1)
arr2 = get_coordinates(img1, 2)

arr = arr1 + arr2

np.savetxt("file_coordinates_listOfCoordinates.txt", arr)

arr1a = get_coordinates(img1a, 1)
arr2a = get_coordinates(img1a, 2)
arr3a = get_coordinates(img1a, 3)
arr4a = get_coordinates(img1a, 4)
arr5a = get_coordinates(img1a, 5)

arra = arr1a + arr2a + arr3a + arr4a + arr5a

np.savetxt("file_coordinates_listOfCoordinatesa.txt", arra)

"""arr_new = [0] * len(arra)

for i in range (len(arra)):
    for j in range(len(arr)):
        if arra[i] == arr[j]:
            arr_new[i] = arr[i]
            #print(arr[i],i)

arr_new = arr_new[arr_new!=0]"""

for i in range(len(arr1)):
    if arr1a[0][0] == arr1[i][0] and arr1a[0][1] == arr1[i][1]:
        print(i, "1,1")
for i in range(len(arr2)):
    if arr1a[0][0] == arr2[i][0] and arr1a[0][1] == arr2[i][1]:
        print(i, "1,2")

for i in range(0, len(arr1)):
    if arr1[i][0] == arr2a[0][0] and arr1[i][1] == arr2a[0][1]:
        print(i, "2,1")
for i in range(len(arr2)):
    if arr2a[0][0] == arr2[i][0] and arr2a[0][1] == arr2[i][1]:
        print(i, "2,2")


for i in range(len(arr1)):
    if arr3a[0][0] == arr1[i][0] and arr3a[0][1] == arr1[i][1]:
        print(i, "3,1")
for i in range(len(arr2)):
    if arr3a[0][0] == arr2[i][0] and arr3a[0][1] == arr2[i][1]:
        print(i, "3,2")

for i in range(len(arr1)):
    if arr4a[0][0] == arr1[i][0] and arr4a[0][1] == arr1[i][1]:
        print(i, "4,1")
for i in range(len(arr2)):
    if arr4a[0][0] == arr2[i][0] and arr4a[0][1] == arr2[i][1]:
        print(i, "4,2")

for i in range(len(arr1)):
    if arr5a[0][0] == arr1[i][0] and arr5a[0][1] == arr1[i][1]:
        print(i, "5,1")
for i in range(len(arr2)):
    if arr5a[0][0] == arr2[i][0] and arr5a[0][1] == arr2[i][1]:
        print(i, "5,2")
