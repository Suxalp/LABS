import numpy as np
import cv2
canvas = np.zeros((400, 400, 3), dtype="uint8")
green = (0, 255, 0)

cv2.line(canvas, (0, 0), (400, 400), green, 5)
cv2.imshow("Canvas", canvas)


red = (0, 0, 255)
cv2.rectangle(canvas, (50, 50), (150, 150), red,5)
cv2.imshow("Canvas", canvas)


blue = (255, 0, 0)

for r in range(0, 175, 25):
    cv2.circle(canvas, (200, 200), 50, blue, -1)

cv2.imshow("Canvas", canvas)
cv2.waitKey(0)





