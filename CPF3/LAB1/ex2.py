import cv2
# Load an image
image = cv2.imread('input.png')
# Specify the new width and height (change these values as needed)
new_width = 300
new_height = 200
# Resize the image
resized_image = cv2.resize(image, (new_width, new_height))
# Save or display the resized image
cv2.imwrite('resized_image.png', resized_image)
cv2.imshow('Resized Image', resized_image)
cv2.waitKey(0)
cv2.destroyAllWindows()