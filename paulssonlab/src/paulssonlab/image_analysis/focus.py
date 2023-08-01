import numpy as np
import skimage.filters


# FROM: Pertuz, S., Puig, D., & Garcia, M. A. (2013). Analysis of focus measure operators for shape-from-focus. Pattern Recognition, 46(5), 1415-1432.
def lap4(img, ksize=3):
    return np.var(skimage.filters.laplace(img, ksize=ksize))


# FROM: Pertuz, S., Puig, D., & Garcia, M. A. (2013). Analysis of focus measure operators for shape-from-focus. Pattern Recognition, 46(5), 1415-1432.
def gra6(img):
    G_x = skimage.filters.sobel_h(img)
    G_y = skimage.filters.sobel_v(img)
    return np.sum(G_x**2 + G_y**2)


# FROM: Pertuz, S., Puig, D., & Garcia, M. A. (2013). Analysis of focus measure operators for shape-from-focus. Pattern Recognition, 46(5), 1415-1432.
def gra7(img):
    G = skimage.filters.sobel(img)
    return np.var(G - G.mean())


# FROM: Pertuz, S., Puig, D., & Garcia, M. A. (2013). Analysis of focus measure operators for shape-from-focus. Pattern Recognition, 46(5), 1415-1432.
def sta3(img):
    return np.var(img - img.mean())
