from skimage.io import imread, imsave
import numpy as np
import pyfftw, cv2
import multiprocessing as mp
from pyfftw.interfaces.scipy_fftpack import ifft2, fft2, fftshift, ifftshift
pyfftw.config.NUM_THREADS = mp.cpu_count()
img_name = "D:/Image location/image_001.tif"
save_name = "D:/Image location/image_001-fft.tif"
image = imread(img_name)
if len(image.shape) > 2:
    image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
#image = cv2.resize(image,(image.shape[1]//2,image.shape[0]//2))
print("Calculating Gaussian Filter")
r_in = 5 #change r_in and r_out to filter out the correct amount of information
#r_out = 0.1*np.sqrt(image.shape[1]*image.shape[0])
r_out = 40
print("R out:",r_out)
filt = 2**np.ceil(np.log2(max(image.shape[0],image.shape[1])))
x, y = np.meshgrid(np.linspace(np.round(-0.5*image.shape[1]),np.round(0.5*image.shape[1]),image.shape[1]),np.linspace(np.round(-0.5*image.shape[0]),np.round(0.5*image.shape[0]),image.shape[0]))
factors = np.multiply(x,x)+np.multiply(y,y)
del x,y
fact_remove = np.where(factors < 1)
scaleLarge = (2*r_out/filt)**2
scaleSmall = (2*r_in/filt)**2
factors = (1-np.exp(-factors*scaleLarge))*np.exp(-factors*scaleSmall)
if (len(fact_remove[0])):
    factors[fact_remove] = 0
factors = factors.astype(np.float32)
print("Performing FFT")
image=fftshift(fft2(image))
print("Multiplying Gaussian Filter with FFT")
image=ifftshift(np.multiply(image,factors))
del factors
print("Performing inverse FFT")
image = ifft2(image)
image = np.absolute(image - np.amin(image))
contrast_enhance = 255/(6*np.std(image))
image *= contrast_enhance
image += 125 - np.mean(image)
image = image.round()
image[np.where(image<0)]=0
image[np.where(image>255)]=255
image = image.astype(np.uint8)
print("Saving image")
imsave(save_name,image)
#cv2.imwrite(save_name,image)