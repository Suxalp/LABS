import cv2
path= 'input.png'
image = cv2.imread(path, cv2.IMREAD_GRAYSCALE)
print(image.shape)
cv2.imshow("Image", image)
cv2.waitKey(0)