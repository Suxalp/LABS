import numpy as np
import argparse
import cv2


ap = argparse.ArgumentParser()
ap.add_argument("-i", "--image", required=True, help="Path to the image")
args = vars(ap.parse_args())

# Load the Image and Convert to Grayscale
image = cv2.imread(args["image"])
gray_image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
cv2.imshow("Grey Scale", gray_image)

# Gaussian Blur (Pre-processing)
blurred_image = cv2.GaussianBlur(gray_image, (3, 3), 0)


# Edge Detection Methods
#   Canny Edge
canny_edges = cv2.Canny(blurred_image, 50, 100)
cv2.imshow("Canny Edge Detection", canny_edges)

#   Laplacian
laplacian_edges = cv2.Laplacian(blurred_image, cv2.CV_64F)
laplacian_edges = np.uint8(np.absolute(laplacian_edges))
cv2.imshow("Laplacian Edge Detection", laplacian_edges)

#   Sobel
sobel_x = cv2.Sobel(blurred_image, cv2.CV_64F, 1, 0, ksize=3)
sobel_y = cv2.Sobel(blurred_image, cv2.CV_64F, 0, 1, ksize=3)
sobel_x = np.uint8(np.absolute(sobel_x))
sobel_y = np.uint8(np.absolute(sobel_y))
sobel_combined = cv2.bitwise_or(sobel_x, sobel_y)
cv2.imshow("Sobel Edge Detection", sobel_combined)

# Wait for a key press to close all windows
cv2.waitKey(0)
cv2.destroyAllWindows()
