import cv2
import numpy as np

# Load image
image = cv2.imread('image5.png', cv2.IMREAD_COLOR)

# Set border size and kernel size
bor = 30  # border size
ker = 50  # blur intensity

# Get size of the image
rows, cols, _ = image.shape

# zero padding border
mask = np.zeros((rows, cols))

# mask border
mask[:bor, :] = 1  # Top
mask[-bor:, :] = 1  # Bottom
mask[:, :bor] = 1  # Left
mask[:, -bor:] = 1  # Right

# Apply blur to the entire image
blurred_image = cv2.blur(image, (ker, ker))

# Convert the mask to 3-channel by stacking it for each color channel (BGR)
mask3 = np.stack([mask]*3, axis=-1)

# Combine blurred border with the original center
result_image = np.where(mask3 == 1, blurred_image, image)

# Display images
cv2.imshow('Original Image', image)
cv2.imshow('Blurred', blurred_image)
cv2.imshow('Blurred Border Image', result_image)
cv2.waitKey(0)
cv2.destroyAllWindows()
