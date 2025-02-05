import numpy as np
import cv2


# Load the image
image = cv2.imread('image.png')
if image is None:
    print(f"Error: Could not load image at {image_path}")
else:
    cv2.imshow("Original", image)

    # Apply bilateral filter with different parameter sets
    blurred = np.hstack([
        cv2.bilateralFilter(image, d=9, sigmaColor=75, sigmaSpace=75),
        cv2.bilateralFilter(image, d=9, sigmaColor=150, sigmaSpace=150),
        cv2.bilateralFilter(image, d=9, sigmaColor=500, sigmaSpace=500)
    ])

    # Display the result
    cv2.imshow("Edge Preserving Blur", blurred)
    cv2.waitKey(0)
    cv2.destroyAllWindows()
