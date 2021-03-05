# GAKTpore: Stereological Characterisation Methods

## Installation
To install the package from PyPi:

```bash
$ pip install GAKTPore
```

## Usage

GAKTPore provides a class 'AnalysePores' with multiple analytical tools. 

### GAKTPore.**AnalysePores**
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

### GAKTPore.AnalysePores.**process**
The next major function is process (**process_parallel** for the multiprocessing version). This function computes the properties of the pores (Area, Circularity etc.).
Note that this must be run before the territory areas can be calculated.

#### **Parameters:** 
##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; FFT: bool, *Optional*. Default Value: True
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Whether to use FFT bandpass to smooth contours. Setting this to False will use the Savgol Filter from Scipy instead (Not validated yet, but much faster).

### GAKTPore.AnalysePores.**process_free_area**
This function (and its parallel counterpart, **process_free_area_parallel**) calculates the territory area for each pore by computing the closest pore contour for each pixel of the image provided.

#### **Parameters:** 
##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; zoom: int, *Optional*. Default Value: 1
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Increase the resolution of the map used for computing the territory area. Example: zoom=2 will use double the resolution of the input image to calculate the territory area.

### GAKTPore.AnalysePores.**process_homogeneity_colour_map**
This function (and its parallel counterpart, **process_homogeneity_colour_map_parallel**) generates a colour map using the area fractions (Pore area divided by Territory area).
#### **Parameters:** 
##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; mapper: matplotlib.colors.LinearSegmentedColormap, *Optional*. Default Value: matplotlib.cm.get_cmap("jet")
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Sets the colourmap to be used when colouring the image. Uses the colour "jet" by default. Provided in this layout to support custom matplotlib maps.
##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; vmin: float, *Optional*. Default Value: 0
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Sets the minimum value to correspond with starting colour of the colourmap. 
##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; vmax: float, *Optional*. Default Value: 1
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Sets the maximum value to correspond with final colour of the colourmap. 
##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; draw_contours: bool, *Optional*. Default Value: False
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Whether to draw the pores onto the map.
  
#### **Returns:**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Colour map of the same resolution as the one in *process_free_area*.

### GAKTPore.AnalysePores.**process_radial_contour**:
  Computes the number of pores and the porosity percentage in segmented steps from the centre of the image.
#### **Parameters:** 

##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; radii_n: int, *Optional*. Default Value: 10
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Number of segmented steps to use between the centre and the maximum radius

##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; radius_centre: np.array, *Optional*. Default Value: None
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The pixel position to use as the centre for the segmented circle. If not supplied, simply takes the centre of the image.

##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; radius_max: float, *Optional*. Default Value: None
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The maximum radius to do the calculations for. If not supplied, simply calculates the distance to edge of the image from the centre.

##### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; draw_contours: bool, *Optional*. Default Value: False
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Whether to draw the pores onto the map.

### Testing Script

A script utilising these functions to output the relevant data into a csv file and generate a colour map is included. This is what was used for the initial GAKTpore paper.

## Citation

When using this package please cite:

*   Sheppard, G. and Tassenberg, K. et al. _GAKTpore: Stereological Characterisation Methods for Porous Foams in Biomedical Applications_, Materials 2021, MDPI.
