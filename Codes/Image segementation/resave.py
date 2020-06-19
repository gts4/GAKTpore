from skimage.io import imread, imshow, imsave
import numpy as np
import os
os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = pow(2,40).__str__()
import cv2
img_name = "D:/Image location/image_001.tif"
save_name = "D:/Image location/image_001-resave.tif"
img = imread(img_name)
img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
imsave(save_name,img)