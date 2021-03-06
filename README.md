# GAKTpore: Stereological Characterisation Methods

## Installation
To install the package from PyPi:

```bash
$ pip install GAKTpore
```

## Quick Use

GAKTpore is distributed with a function utilising the GAKTpore 'AnalysePores' class to output the relevant data files into a save folder and generate a colour map. This outputs the data used for the initial GAKTpore paper. Test-case files are included on Github, with images.

The syntax of the function are as follows:

### GAKTpore.Analyse.**run**

#### **Parameters**

##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; IMAGE_NAME: str,
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The directory path for the image to be analysed. The image is expected to be grayscale. The image formats supported can be found on the OpenCV [imread](https://docs.opencv.org/3.4/d4/da8/group__imgcodecs.html#imread) reference.

##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; SAVE_FOLDER:str,
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The directory where the output files will be saved to. There are 3 output files by default:
* A file containing the standard evaluation parameters (such as mean diameter, number of samples etc.) per radial step.
* A file containing the standard evaluation parameters (such as pore area, cirularity etc.) for each pore, allowing further analysis.
* The Area fraction image generated using the 'jet' colourmap via Matplotlib.

The filename will be determined by the name of the image file and input parameters (Threshold, upscale multiplier, contour filtering method (FFT or Savgol) and the data) as &lt;Image name&gt;-&lt;Threshold&gt;-&lt;Upscale Multiplier&gt;-&lt;type of data&gt;_&lt;filtering method&gt;.file-type, for example 'Test1-127-1x_FFT.csv'.

##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; THRES: int,
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The binarisation integer, expected to be between 0-255.

##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; SCALE: float,
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Distance per pixel, presumed to be taken from the 'scale' bar of a microstructure image.

##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; UPSCALE_MULTIPLIER: int, *Optional*. Default Value = 1
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The number to upscale the image being processed. Will multiply the image resulution by the value given.

##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; W_BG: bool, *Optional*. Default Value = True, 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Whether the background of the image is white (True) or black (False).

##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; FFT: bool, *Optional*. Default Value = True, 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Whether to use FFT bandpass to smooth contours. Setting this to False will use the Savgol Filter from Scipy instead (Not validated yet, but much faster).

##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; parallel: bool, *Optional*. Default Value = True, 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Whether to use multiple cpu cores to process the image in parallel.

##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; cpu_count: int, *Optional*. Default Value = -1,
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The number of parallel computations used when a multiprocessing function is used. -1 to use the number of available cores.

##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; npy_save: bool, *Optional*. Default Value = False,
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Whether to save the numpy array of the territory map. This map is a mask containing integers, where the integers refer to the pore number the territory belongs to.

##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; plt_save:bool, *Optional*. Default Value = False,
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Whether to use matplotlib.save to save the area fraction image. Matplotlib.save is slow and will crash when using extremely high resolution images (> 25k x 25k), but is otherwise stable. If set to false, will use cv2.imsave to save the area fraction image.

##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; draw_contours:bool, *Optional*. Default Value=True,
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Whether to draw the pores onto the map in black.

##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; vmin:float, *Optional*. Default Value=0,
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; For the area fraction colour map, sets the minimum value to correspond with starting colour of the colour map. 

##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; vmax:float, *Optional*. Default Value=1,
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; For the area fraction colour map, sets the maximum value to correspond with final colour of the colour map. 

##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; radii_n:int, *Optional*. Default Value=10
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Number of segmented steps to use between the centre of the image and the maximum radius.

## Usage

GAKTpore provides a class 'AnalysePores' with multiple analytical tools. 

### GAKTpore.**AnalysePores**
The initialisation for the class works as a simple binarisation and pore detection tool utilising the OpenCV implementation of findContours.

#### **Parameters:** 

##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; img: np.array,
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 2D Grayscale Image
##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; threshold_value: int, *Optional*. Default Value: 125
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Threshold value for binarising the image
##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; scale = float, *Optional*. Default Value: 1
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Distance per pixel, presumed to be taken from the 'scale' bar of a microstructure image.
##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; G = bool, *Optional*. Default Value:False
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Whether to apply a Gaussian filter of sigma=2
##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; white_background= bool, *Optional*. Default Value: True
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Whether the background of the image is white (True) or black (False)
##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; cpu_count= int, *Optional*. Default Value:-1
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The number of parallel computations used when a multiprocessing function is used. -1 to use the number of available cores.

### GAKTpore.AnalysePores.**process**
The next major function is process (**process_parallel** for the multiprocessing version). This function computes the properties of the pores (Area, Circularity etc.).
Note that this must be run before the territory areas can be calculated.

#### **Parameters:** 
##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; FFT: bool, *Optional*. Default Value: True
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Whether to use FFT bandpass to smooth contours. Setting this to False will use the Savgol Filter from Scipy instead (Not validated yet, but much faster).

### GAKTpore.AnalysePores.**process_free_area**
This function (and its parallel counterpart, **process_free_area_parallel**) calculates the territory area for each pore by computing the closest pore contour for each pixel of the image provided.

#### **Parameters:** 
##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; zoom: int, *Optional*. Default Value: 1
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Increase the resolution of the map used for computing the territory area. Example: zoom=2 will use double the resolution of the input image to calculate the territory area.

### GAKTpore.AnalysePores.**process_homogeneity_colour_map**
This function (and its parallel counterpart, **process_homogeneity_colour_map_parallel**) generates a colour map using the area fractions (Pore area divided by Territory area).
#### **Parameters:** 
##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; mapper: matplotlib.colors.LinearSegmentedColormap, *Optional*. Default Value: matplotlib.cm.get_cmap("jet")
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Sets the colourmap to be used when colouring the image. Uses the colour "jet" by default. Provided in this layout to support custom matplotlib maps.
##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; vmin: float, *Optional*. Default Value: 0
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Sets the minimum value to correspond with starting colour of the colour map. 
##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; vmax: float, *Optional*. Default Value: 1
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Sets the maximum value to correspond with final colour of the colour map. 
##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; draw_contours: bool, *Optional*. Default Value: False
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Whether to draw the pores onto the map.
  
#### **Returns:**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Colour map of the same resolution as the one in *process_free_area*.

### GAKTpore.AnalysePores.**process_radial_contour**:
  Computes the number of pores and the porosity percentage in segmented steps from the centre of the image.
#### **Parameters:** 

##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; radii_n: int, *Optional*. Default Value: 10
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Number of segmented steps to use between the centre and the maximum radius.

##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; radius_centre: np.array, *Optional*. Default Value: None
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The pixel position to use as the centre for the segmented circle. If not supplied, simply takes the centre of the image.

##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; radius_max: float, *Optional*. Default Value: None
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The maximum radius to do the calculations for. If not supplied, simply calculates the distance to edge of the image from the centre.

##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; draw_contours: bool, *Optional*. Default Value: False
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Whether to draw the pores onto the map.

## Citation

When using this package please cite:

*   Sheppard, G. and Tassenberg, K. et al. _GAKTpore: Stereological Characterisation Methods for Porous Foams in Biomedical Applications_, Materials 2021, MDPI.
