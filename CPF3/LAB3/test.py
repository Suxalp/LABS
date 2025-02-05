import numpy as np
import cv2
import conv2d

# 1. GREYSCALE
image = cv2.imread("image.png")
cv2.imshow("Original", image)
image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
cv2.imshow("Original", image)

# 2. GAUSSIAN BLUR
blurred = np.hstack([
cv2.GaussianBlur(image, (3, 3), 0),
cv2.GaussianBlur(image, (5, 5), 0),
cv2.GaussianBlur(image, (7, 7), 0)])
cv2.imshow("Gaussian", blurred)
cv2.waitKey(0)

Kx = np.array([[-1, 0, 1], [-2, 0, 2], [-1, 0, 1]], np.float32)
Ky = np.array([[1, 2, 1], [0, 0, 0], [-1, -2, -1]], np.float32)

    # Convolve kernels with image
Ix = conv2d(img, Kx)
Iy = conv2d(img, Ky)

    # Calculate gradient magnitude and direction
G = np.hypot(Ix, Iy)
    # Normalize gradient values to 0-255
G = G / G.max() * 255
theta = np.arctan2(Iy, Ix)
   