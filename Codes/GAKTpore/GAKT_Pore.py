#Import all the important packages we're going to need
import cv2, psutil, sys
import numpy as np
from scipy.ndimage import gaussian_filter, gaussian_filter1d
from scipy.spatial.distance import cdist
from scipy.spatial import cKDTree
from functools import partial
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar
import pyfftw
from pyfftw.interfaces.scipy_fftpack import ifft, fft
import matplotlib.path as mpltPath
import matplotlib as mpl
import multiprocessing as mp
from itertools import groupby
from skimage.measure import regionprops
from skimage.color import label2rgb
from skimage.draw import polygon
from types import ModuleType, FunctionType
from gc import get_referents
#AnalysePores.conts #Contours
#AnalysePores.l_rad #List of all the largest radius through pore relating to contour indices (e.g. Value for position 0 is the largest diameter, relating to contour 0)
#AnalysePores.area_rad #Radius of circle is calculated assuming the area is that of a perfect circle (sqrt(area/pi))
#AnalysePores.circ #List of all the circularity relating to contour indices
#AnalysePores.waviness #List of all the waviness relating to contour indices
#AnalysePores.aspect_ratios #List of all the aspect ratios relating to contour indices
#AnalysePores.area #List of all the areas relating to contour indices
#AnalysePores.FFT_conts # FFT'd Contours
#AnalysePores.total_porosity #Total porosity
#AnalysePores.homogeneity_map #Map of contour space - each x,y coordinate links to a contour indice
#AnalysePores.homogeneity_zoom # How much is the map enlarged by compared to original image - by default it is set to 1.
#AnalysePores.homogeneity_areas #Free space area
#AnalysePores.radii_n #Number of rings
#AnalysePores.radius_centre #Centre of the rings
#AnalysePores.radius_max #Final radius
#AnalysePores.radius_contours #The contour indices in each ring. For example: AnalysePores.area[AnalysePores.radius_contours[0]] will give you the areas in the first ring.
#AnalysePores.radius_step_len #The length of each step in radii. For example: the second ring will be within 1*step_len and 2*step_len.
#AnalysePores.radius_porosity #Porosity for rings calculated from binary image.
#AnalysePores.radius_porosity_num #Number of black pixels used to calculate means
#AnalysePores.radius_porosity_num #Number of pixels used to calculate means
#AnalysePores.homogeneity_colour_map #A RGB colour map based on volume fraction of pore. Must be made first using process_homogeneity_colour_map(). draw_contours=True option can be used to display the contours used to make the map.


class AnalysePores():
	def __init__(self, img, threshold_value = 125, scale = 1, G = False, white_background=True):
		self.scale = scale
		bimg = img.copy()
		if len(img.shape) > 2:
			bimg = cv2.cvtColor(bimg, cv2.COLOR_BGR2GRAY)
		if white_background:
			bimg = (bimg <= threshold_value).astype(np.uint8)*255
		else:
			bimg = (bimg >= threshold_value).astype(np.uint8)*255
		self.binary_img = bimg
		if G:
			bimg = (gaussian_filter(bimg,sigma=2) > threshold_value).astype(np.uint8)*255
		bimg_analysis = np.zeros(bimg.shape+np.array([10,10]),dtype=np.uint8)
		if not white_background:
			bimg_analysis += 255
		bimg_analysis[5:-5,5:-5] = bimg
		contours, hierarchy = cv2.findContours(bimg_analysis, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)
		main_contours = []
		main_hierarchy = hierarchy[0,:,:]
		main_index = []
		bimg_offset = np.array([5,5])
		img_size = img.shape[0]*img.shape[1]
		for i in range(len(contours)):
			if cv2.contourArea(contours[i]) > cv2.arcLength(contours[i], True) and cv2.contourArea(contours[i]) < 0.8*img_size:
				main_index += [i]
				main_contours += [contours[i]-bimg_offset]
				#and cv2.contourArea(contours[i]) < 0.85*img_size:

		main_contours = np.array(main_contours)
		main_index = np.array(main_index)
		main_binary = []
		for n in main_contours:
			x = n.T[0][0]
			y = n.T[1][0]
			dx = np.roll(x,-1)-x
			dy = np.roll(y,-1)-y
			main_binary += [np.sum(y*dx - x*dy) >= 0]
		main_binary = np.array(main_binary)
		
		pores = np.where(main_binary)[0]
		main_pores = main_index[pores]
		holes = np.where(main_binary == False)[0]
		parents = main_hierarchy[main_index[holes]][:,3]
		
		remove_list = np.where(np.isin(parents,main_pores) == False)[0]
		if len(remove_list) > 0:
			main_contours = np.delete(main_contours,holes[remove_list])
			main_index = np.delete(main_index,holes[remove_list])
			main_binary = np.delete(main_binary,holes[remove_list])
			pores = np.where(main_binary)[0]
			holes = np.where(main_binary == False)[0]
			parents = main_hierarchy[main_index[holes]][:,3]

		cont_pore_index = np.zeros(len(main_contours),dtype=int)
		cont_pore_index[main_binary] = np.arange(len(pores))
		main_pore_index = []
		for p in pores:
			main_pore_index += [[p]]
		main_pores = main_index[pores]
		add_list = np.isin(parents,main_pores)

		the_tree = cKDTree(np.array([main_pores,np.zeros(len(main_pores))]).T)
		_, nn = the_tree.query(np.array([parents[add_list],np.zeros(len(parents[add_list]))]).T)
		add_list = holes[add_list]
		cont_pore_index[add_list] = nn
		for i in range(len(add_list)):
			main_pore_index[nn[i]] += [add_list[i]]

		self.conts = main_contours
		self.main_pore_index = np.array(main_pore_index)
		self.cont_pore_index = np.array(cont_pore_index)
		self.conts_centre = []
		for c in main_contours:
			self.conts_centre += [contourCentre(c)]
		self.conts_centre = np.array(self.conts_centre)
		self.pore_centre = []
		self.l_rad = []
		self.area_rad = []
		self.cont_circ = []
		self.cont_waviness = []
		self.cont_aspect_ratios = []
		self.circ = []
		self.waviness = []
		self.aspect_ratios = []
		self.cont_area = []
		self.area = []
		self.FFT_conts = []
		self.total_porosity = -1
		self.homogeneity_map = []
		self.homogeneity_areas =[]
		self.homogeneity_zoom = 1
		self.homogeneity_colour_map = []
		self.homogeneity_colour_mapper = -1
		self.radii_n = -1
		self.radius_centre = -1
		self.radius_max = -1
		self.radius_contours = -1
		self.radius_step_len = -1
		self.radius_steps = -1
		self.radius_porosity = -1
		self.radius_porosity_num = -1
		self.radius_porosity_total_num = -1
	def process(self):
		x,y = np.meshgrid(np.linspace(0,1,30),np.linspace(0,1,30))
		x = np.concatenate(x)
		y = np.concatenate(y)
		grid = np.vstack((x, y)).T
		del x,y
		for cont in self.conts:
			props = calc_contour_properties(cont)
			self.cont_circ += [props[0]]
			self.cont_waviness += [props[1]]
			self.cont_aspect_ratios += [props[2]]
			self.cont_area += [props[3]]
			self.FFT_conts += [props[4]]
		self.total_porosity = calc_poro(self.binary_img)
		self.cont_circ = np.array(self.cont_circ)
		self.cont_waviness = np.array(self.cont_waviness)
		self.cont_aspect_ratios = np.array(self.cont_aspect_ratios)
		self.cont_area = np.array(self.cont_area)*self.scale**2
		self.FFT_conts = np.array(self.FFT_conts)
		for p in self.main_pore_index:
			self.l_rad += [calc_largest_rad(self.FFT_conts[p],g=grid)]
			self.circ += [self.cont_circ[p[0]]]
			self.waviness += [self.cont_waviness[p[0]]]
			self.aspect_ratios += [self.cont_aspect_ratios[p[0]]]
			self.area += [self.cont_area[p[0]]]
			self.pore_centre += [self.conts_centre[p[0]]]
			if len(p) > 1:
				for i in range(1,len(p)):
					self.area[-1] -= self.cont_area[p[i]]
			self.area_rad += [np.sqrt(self.area[-1]/np.pi)]
		self.l_rad = np.array(self.l_rad)*self.scale
		self.area_rad = np.array(self.area_rad)*self.scale
		self.circ = np.array(self.circ)
		self.waviness = np.array(self.waviness)
		self.aspect_ratios = np.array(self.aspect_ratios)
		self.area = np.array(self.area)*self.scale**2
		self.pore_centre = np.array(self.pore_centre)
		
	def process_parallel(self):
		cpu_count = mp.cpu_count()
		pool = mp.Pool(processes=cpu_count)
		x,y = np.meshgrid(np.linspace(0,1,30),np.linspace(0,1,30))
		x = np.concatenate(x)
		y = np.concatenate(y)
		grid = np.vstack((x, y)).T
		del x,y
		data = np.array(pool.map(calc_contour_properties, iter(self.conts)))
		for props in data:
			self.cont_circ += [props[0]]
			self.cont_waviness += [props[1]]
			self.cont_aspect_ratios += [props[2]]
			self.cont_area += [props[3]]
			self.FFT_conts += [props[4]]
		self.total_porosity = calc_poro(self.binary_img)
		self.cont_circ = np.array(self.cont_circ)
		self.cont_waviness = np.array(self.cont_waviness)
		self.cont_aspect_ratios = np.array(self.cont_aspect_ratios)
		self.cont_area = np.array(self.cont_area)*self.scale**2
		self.FFT_conts = np.array(self.FFT_conts)
		FFT_Pores = []
		for p in self.main_pore_index:
			FFT_Pores += [self.FFT_conts[p]]
			self.circ += [self.cont_circ[p[0]]]
			self.waviness += [self.cont_waviness[p[0]]]
			self.aspect_ratios += [self.cont_aspect_ratios[p[0]]]
			self.area += [self.cont_area[p[0]]]
			self.pore_centre += [self.conts_centre[p[0]]]
			if len(p) > 1:
				for i in range(1,len(p)):
					self.area[-1] -= self.cont_area[p[i]]
			self.area_rad += [np.sqrt(self.area[-1]/np.pi)]
			#self.circ += [self.cont_circ[p[0]]*self.area[-1]/self.cont_area[p[0]]]
		self.l_rad = np.array(pool.map(partial(calc_largest_rad,g=grid), iter(FFT_Pores)))*self.scale
		self.circ = np.array(self.circ)
		self.waviness = np.array(self.waviness)
		self.aspect_ratios = np.array(self.aspect_ratios)
		self.area = np.array(self.area)
		self.area_rad = np.array(self.area_rad)
		self.pore_centre = np.array(self.pore_centre)
		pool.close()
	def process_free_area(self,zoom=1):
		assert len(self.FFT_conts) > 0,"No processed contours found!"
		lens = np.cumsum(np.array([len(c) for c in self.FFT_conts]))-1
		the_tree = cKDTree(np.concatenate(self.FFT_conts),leafsize=1)
		new_len_x = (self.binary_img.shape[1]-1)*zoom+1
		new_len_y = (self.binary_img.shape[0]-1)*zoom+1
		pts = np.meshgrid(np.linspace(0,self.binary_img.shape[1]-1,new_len_x),np.linspace(0,self.binary_img.shape[0]-1,new_len_y))
		pts = np.concatenate(np.array(pts).T)
		_, nn = the_tree.query(pts)
		the_tree = cKDTree(np.array([lens,np.zeros(len(lens))]).T,leafsize=1)
		_, nn_cc = the_tree.query(np.array([nn,np.zeros(len(nn))]).T)
		nn_cc = np.array((the_tree.data[nn_cc][:,0]<nn).astype(int)+nn_cc)
		nn_map = np.zeros((new_len_y,new_len_x),dtype=np.uint32)
		nn_map[tuple((pts.T*zoom).astype(int))[::-1]] = self.cont_pore_index[nn_cc]
		self.homogeneity_areas = []
		self.homogeneity_map = nn_map.copy()
		for region in regionprops(nn_map+1):
			self.homogeneity_areas += [region.area]
		self.homogeneity_areas = np.array(self.homogeneity_areas)*self.scale**2/(zoom)**2
		self.homogeneity_zoom = zoom
	def process_free_area_parallel(self,zoom=1):
		assert len(self.FFT_conts) > 0,"No processed contours found!"
		lens = np.cumsum(np.array([len(c) for c in self.FFT_conts]))-1
		the_tree = cKDTree(np.concatenate(self.FFT_conts),leafsize=1)
		cpu_count = min(mp.cpu_count(),int(psutil.virtual_memory().available/getsize(the_tree)))
		pool = mp.Pool(processes=cpu_count)
		new_len_x = (self.binary_img.shape[1]-1)*zoom+1
		new_len_y = (self.binary_img.shape[0]-1)*zoom+1
		pts = np.meshgrid(np.linspace(0,self.binary_img.shape[1]-1,new_len_x),np.linspace(0,self.binary_img.shape[0]-1,new_len_y))
		pts = np.concatenate(np.array(pts).T)
		nn = np.concatenate(pool.map(partial(nn_proc,the_tree=the_tree), iter([pts[i:i + int(len(pts)//cpu_count)] for i in range(0, len(pts), int(len(pts)//cpu_count))])))
		chunks = int(len(nn)//cpu_count)
		chk_list = [nn[i:i + chunks] for i in range(0, len(nn), chunks)]
		the_tree = cKDTree(np.array([lens,np.zeros(len(lens))]).T,leafsize=1)
		nn_cc = np.concatenate(pool.map(partial(chk_proc,the_tree = the_tree), iter(chk_list)))
		nn_map = np.zeros((new_len_y,new_len_x),dtype=np.uint32)
		nn_map[tuple((pts.T*zoom).astype(int))[::-1]] = self.cont_pore_index[nn_cc]
		#arange = np.arange(len(lens))
		self.homogeneity_areas = []
		for region in regionprops(nn_map+1):
			self.homogeneity_areas += [region.area]
		#self.homogeneity_areas = np.concatenate(pool.map(partial(area_proc,nn_map=nn_map), iter([arange[i:i + int(len(arange)//cpu_count)] for i in range(0, len(arange), int(len(arange)//cpu_count))])))
		self.homogeneity_map = nn_map
		self.homogeneity_areas = np.array(self.homogeneity_areas)*self.scale**2/zoom**2
		self.homogeneity_zoom = zoom
		pool.close()
	def process_radial_contour(self,radii_n=10,radius_centre = False,radius_max=False):
		if type(radius_centre) == bool:
			radius_centre = np.array(self.binary_img.shape[::-1])//2
		if type(radius_max) == bool:
			radius_max = min(self.binary_img.shape[::-1]-np.array(radius_centre))
		step_len = radius_max/radii_n
		step = np.arange(0,radii_n)*step_len
		contour_dists = np.sum((self.pore_centre - radius_centre)**2,axis=1)
		x,y = np.meshgrid(np.linspace(0,self.binary_img.shape[1]-1,self.binary_img.shape[1]),np.linspace(0,self.binary_img.shape[0]-1,self.binary_img.shape[0]))
		x = (x-radius_centre[0])**2
		y = (y-radius_centre[1])**2
		rad_map = x+y
		del x,y
		contour_step = []
		radius_porosity = []
		radius_porosity_num = []
		radius_porosity_total_num = []
		for i in range(1,len(step)):
			contour_step += [np.where((contour_dists < step[i]**2)&(contour_dists >= step[i-1]**2))[0]]
			vals = 1-self.binary_img[(rad_map < step[i]**2)&(rad_map >= step[i-1]**2)]/255
			radius_porosity += [np.mean(vals)]
			radius_porosity_num += [np.sum(vals)]
			radius_porosity_total_num += [len(vals)]
		contour_step += [np.where((contour_dists >= step[-1]**2))[0]]
		vals = 1-self.binary_img[(rad_map >= step[-1]**2)]/255
		radius_porosity += [np.mean(vals)]
		radius_porosity_num += [np.sum(vals)]
		radius_porosity_total_num += [len(vals)]
		self.radii_n = radii_n
		self.radius_centre = radius_centre
		self.radius_max = radius_max
		self.radius_contours = contour_step
		self.radius_step_len = step_len*self.scale
		self.radius_steps = step*self.scale
		self.radius_porosity = radius_porosity
		self.radius_porosity_num = radius_porosity_num
		self.radius_porosity_total_num = radius_porosity_total_num
	def process_homogeneity_colour_map(self,mapper=False,vmin=0,vmax=1,draw_contours=False):
		norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
		if type(mapper) == bool:
			mapper = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.get_cmap("terrain"))
		vol_frac = self.area/self.homogeneity_areas
		colour_array = []
		for i in range(len(self.area)):
			colour_array += [np.round(np.array(mapper.to_rgba(vol_frac[i])[:3])*255).astype(int)]
		colour_array = np.array(colour_array)
		vol_frac_map = label2rgb(self.homogeneity_map+1,colors=colour_array).astype(int)
		if draw_contours:
			cv2.fillPoly(vol_frac_map, pts=self.conts, color=[0,0,0])
		self.homogeneity_colour_map = vol_frac_map
		self.homogeneity_colour_mapper = mapper
	def process_homogeneity_colour_map_parallel(self,mapper=False,vmin=0,vmax=1,draw_contours=False):
		norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
		if type(mapper) == bool:
			mapper = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.get_cmap("jet"))
		vol_frac = self.area/self.homogeneity_areas
		cpu_count = mp.cpu_count()
		pool = mp.Pool(processes=cpu_count)
		arange = np.arange(len(self.area))
		colour_array = np.concatenate(pool.map(partial(colour_proc,mapper=mapper,vol_frac=vol_frac), iter([arange[i:i + int(len(arange)//cpu_count)] for i in range(0, len(arange), int(len(arange)//cpu_count))])))
		vol_frac_map = label2rgb(self.homogeneity_map+1,colors=colour_array).astype(int)
		if draw_contours:
			cv2.fillPoly(vol_frac_map, pts=self.conts, color=[0,0,0])
		self.homogeneity_colour_map = vol_frac_map
		self.homogeneity_colour_mapper = mapper
		pool.close()
		
def cont_loc_proc(conts):
	rcs = []
	for i in conts:
		r = np.concatenate(i)
		c = r[:,1]
		r = r[:,0]
		rcs += [np.array(polygon(c, r)).T]
	return np.concatenate(rcs)

def colour_proc(conts,mapper,vol_frac):
	colour_array = []
	for i in conts:
		colour_array += [np.round(np.array(mapper.to_rgba(vol_frac[i])[:3])*255).astype(int)]
	return colour_array

def area_proc(conts,nn_map):
	homogeneity_areas = []
	for n in conts:
		map_c = (nn_map == n)
		wh = np.where(map_c)
		ymin, xmin = np.amin(wh,axis=1)
		ymax, xmax = np.amax(wh,axis=1)
		map_c = map_c[max(ymin-1,0):ymax+1,max(xmin-1,0):xmax+1].astype(np.uint8)
		tc, _ = cv2.findContours(map_c, cv2.RETR_LIST, cv2.CHAIN_APPROX_NONE)
		lens = [len(n) for n in tc]
		homogeneity_areas += [cv2.contourArea(tc[np.argmax(lens)])]
	return np.array(homogeneity_areas)

def chk_proc(chk,the_tree):
	_, nn = the_tree.query(np.array([chk,np.zeros(len(chk))]).T)
	return (the_tree.data[nn][:,0]<chk).astype(int)+nn

def nn_proc(pts,the_tree):
	_, nn = the_tree.query(pts)
	return nn

#Return: Circularity, Waviness, Aspect Ratio, Accurate Area, X-Y Points of contour
def calc_contour_properties(cont):
	x = np.concatenate(cont[:,:,0])
	x = np.append(x,x[0])
	y = np.concatenate(cont[:,:,1])
	y = np.append(y,y[0])
	points = np.array([x,y]).T
	distance = np.cumsum( np.sqrt(np.sum( np.diff(points, axis=0)**2, axis=1 )) )
	distance = np.insert(distance, 0, 0)/distance[-1]
	if len(x) < 1000:
		weight = len(x)/1000
		new_s = np.linspace(0, 1, 1000)
		f = interp1d(distance, points, kind="linear", axis=0)
		x,y = f(new_s).T
		del f
	else:
		weight=1
	del points
	x = fft(x)
	y = fft(y)
	grp_len = np.array([len(list(j)) for i, j in groupby(np.concatenate(cont[:,:,0])-np.roll(np.concatenate(cont[:,:,0]),-1))])
	grp_len = np.array([sum(list(j)) for i, j in groupby(grp_len)])
	grp_len_y = np.array([len(list(j)) for i, j in groupby(np.concatenate(cont[:,:,1])-np.roll(np.concatenate(cont[:,:,1]),-1))])
	grp_len_y = np.array([sum(list(j)) for i, j in groupby(grp_len_y)])
	orig_x = np.concatenate(cont[:,:,0])
	orig_y = np.concatenate(cont[:,:,1])
	d_x = np.array([np.array(list(j)) for i, j in groupby(gaussian_filter1d(np.roll(orig_x,-1)-orig_x,sigma=1))])
	d_x = np.array([float(1 == i[0])*len(i) for i in d_x])
	d_x = d_x[d_x > 2]
	d_x_neg = np.array([np.array(list(j)) for i, j in groupby(gaussian_filter1d(np.roll(orig_x,-1)-orig_x,sigma=1))])
	d_x_neg = np.array([float(-1 == i[0])*len(i) for i in d_x_neg])
	d_x_neg = d_x_neg[d_x_neg > 2]
	d_x_0 = np.array([np.array(list(j)) for i, j in groupby(gaussian_filter1d(np.roll(orig_x,-1)-orig_x,sigma=1))])
	d_x_0 = np.array([(1-abs(i[0]))*len(i) for i in d_x_0])
	d_x_0 = d_x_0[d_x_0 > 2]
	comb_x = np.concatenate([d_x_0,d_x,d_x_neg]) 
	d_y = np.array([np.array(list(j)) for i, j in groupby(gaussian_filter1d(np.roll(orig_y,-1)-orig_y,sigma=1))])
	d_y = np.array([float(1 == i[0])*len(i) for i in d_y])
	d_y = d_y[d_y > 2]
	d_y_neg = np.array([np.array(list(j)) for i, j in groupby(gaussian_filter1d(np.roll(orig_y,-1)-orig_y,sigma=1))])
	d_y_neg = np.array([float(-1 == i[0])*len(i) for i in d_y_neg])
	d_y_neg = d_y_neg[d_y_neg > 2]
	d_y_0 = np.array([np.array(list(j)) for i, j in groupby(gaussian_filter1d(np.roll(orig_y,-1)-orig_y,sigma=1))])
	d_y_0 = np.array([(1-abs(i[0]))*len(i) for i in d_y_0])
	d_y_0 = d_y_0[d_y_0 > 2]
	comb_y = np.concatenate([d_y_0,d_y,d_y_neg]) 
	if len(comb_x) > 2:
		low_pass_x = max(np.mean(np.sort(comb_x)[:-len(comb_x)//2])/2,3)
	else:
		low_pass_x = 3
	fft_wavelength = np.fft.fftfreq(len(x),weight)
	fft_wavelength[0] = 1
	fft_wavelength = abs(1/fft_wavelength)
	fft_wavelength[0] = np.inf
	low_pass_x = np.sum(fft_wavelength >= low_pass_x)//2
	if len(comb_y) > 2:
		low_pass_y = max(np.mean(np.sort(comb_y)[:-len(comb_y)//2])/2,3)
	else:
		low_pass_y = 3
	low_pass_y = np.sum(fft_wavelength >= low_pass_y)//2
	x[low_pass_x:-low_pass_x] = 0
	y[low_pass_y:-low_pass_y] = 0
	freq = np.fft.fftfreq(len(x),weight)*2j*np.pi
	dx = ifft(x*freq).real
	dy = ifft(y*freq).real
	x = ifft(x).real
	y = ifft(y).real
	dists = np.cumsum( np.insert(np.sqrt(dx**2+dy**2),0,0)[:-1] )
	dists = dists/dists[-1]
	f = interp1d(dists, np.array([x,y]).T, kind="linear", axis=0)
	x,y = f(distance).T
	x = gaussian_filter1d(x[1:],sigma=1)
	y = gaussian_filter1d(y[1:],sigma=1)
	ratio_x = (max(np.concatenate(cont[:,:,0]))-min(np.concatenate(cont[:,:,0])))/(np.amax(x)-np.amin(x))
	ratio_y = (max(np.concatenate(cont[:,:,1]))-min(np.concatenate(cont[:,:,1])))/(np.amax(y)-np.amin(y))
	x = ratio_x*(x-np.amin(x))+min(np.concatenate(cont[:,:,0]))
	y = ratio_y*(y-np.amin(y))+min(np.concatenate(cont[:,:,1]))
	dx = (np.roll(x,-1)-np.roll(x,1))/2
	dy = (np.roll(y,-1)-np.roll(y,1))/2
	area = np.abs(0.5*np.sum(y*dx - x*dy))
	circ = 4*np.pi*area/np.sum( np.sqrt(np.sum( np.array([dx**2,dy**2]).T, axis=1 )) )**2
	dists = np.sqrt(dx**2+dy**2)
	dx = dx/dists
	dy = dy/dists
	wavs = ((dx*np.roll(dy,-1) - np.roll(dx,-1)*dy) >= -0.03489949670250097) #np.sin(2*np.pi/180) // 2 degrees
	wavs = np.sum( dists[wavs] )/np.sum( dists )
	return circ, wavs, aspect_ratio(cont), area, np.array([x,y]).T # Largest Diameter through pore, Diameter from Area, Circularity, Waviness, Aspect Ratio, Accurate Area, X-Y Points of contour

def calc_largest_rad(fft_c,g):
	x = fft_c[0][:,0]
	y = fft_c[0][:,1]
	xmin = np.amin(x)
	xmax = np.amax(x)
	ymin = np.amin(y)
	ymax = np.amax(y)
	grid = g*(xmax-xmin,ymax-ymin)+(xmin,ymin)
	coords = np.vstack((x, y)).T
	total_coords = np.concatenate(fft_c)
	poly = mpltPath.Path(coords)
	get_pts = poly.contains_points(grid)
	grid = grid[get_pts]
	for i in range(1,len(fft_c)):
		poly_o = mpltPath.Path(fft_c[i])
		get_pts = poly_o.contains_points(grid)
		grid = grid[(get_pts == False)]
	the_tree = cKDTree(total_coords,leafsize=1)
	ds, _ = the_tree.query(grid)
	ds = grid[np.argmax(ds)]
	x_frac = (xmax-xmin)/len(np.unique(grid[:,0]))
	y_frac = (ymax-ymin)/len(np.unique(grid[:,1]))
	xmin = ds[0]-x_frac
	xmax = ds[0]+x_frac
	ymin = ds[1]-y_frac
	ymax = ds[1]+y_frac
	grid = g*(xmax-xmin,ymax-ymin)+(xmin,ymin)
	get_pts = poly.contains_points(grid)
	grid = grid[np.where(get_pts)[0]]
	for i in range(1,len(fft_c)):
		poly_o = mpltPath.Path(fft_c[i])
		get_pts = poly_o.contains_points(grid)
		grid = grid[(get_pts == False)]
	the_tree = cKDTree(total_coords,leafsize=1)
	ds, _ = the_tree.query(grid)
	largest_rad = np.amax(ds)
	return largest_rad


def aspect_ratio(c):
	centre = contourCentre(c)
	new_c = np.array([np.concatenate(c[:,:,0]),np.concatenate(c[:,:,1])]).T-centre
	def min_boundary(ang):
		rotated_c = new_c[:,0]*np.cos(ang)-new_c[:,1]*np.sin(ang)
		return max(rotated_c)-min(rotated_c)
	min_ang = minimize_scalar(min_boundary,bounds=(-np.pi,np.pi)).x
	rotated_c = [new_c[:,0]*np.cos(min_ang)-new_c[:,1]*np.sin(min_ang),new_c[:,0]*np.sin(min_ang)+new_c[:,1]*np.cos(min_ang)]
	return (max(rotated_c[0])-min(rotated_c[0]))/(max(rotated_c[1])-min(rotated_c[1])) 

def calc_poro(img):
	total_area = np.sum(img/(255))
	return 100*total_area/(img.shape[0]*img.shape[1])

def contourCentre(c):
	M = cv2.moments(c)
	cX = int(M["m10"] / M["m00"])
	cY = int(M["m01"] / M["m00"])
	return [cX,cY]

#This equation is found online to calculate the contourCentre
#We then define a function to find Contours within our specified radii
def contourRadius(Contour_centres,centre,r_in,r_out):
	p = np.sqrt(np.sum((Contour_centres-centre)**2,axis=1))
	return np.where((p>=r_in)&(p<=r_out))[0]

#This calculates the distance from the centre and returns
#the indicies of the contours that are within range.
#We define a function to give all the areas of a given contour list
# for use later on.
def contourAreas(Contours):    
	areas = []
	for c in Contours:
		areas += [cv2.contourArea(c)]
	return np.array(areas)

#Example code to loop through radii and plot mean area
def calc_mean_area(main_contours):
    if len(main_contours)>0:
        areas = contourAreas(main_contours)
        mean_area = np.mean(areas)
        std_area = np.std(areas)
    else:
        mean_area = 0
        std_area = 0
    return np.array(std_area),np.array(mean_area)

def accurate_area(c):
    x = np.concatenate(c[:,:,0])
    y = np.concatenate(c[:,:,1])
    #x = np.append(x,x[0])
    #y = np.append(y,y[0])
    points = np.array([x,y]).T
    distance = np.cumsum( np.sqrt(np.sum( np.diff(points, axis=0)**2, axis=1 )) )
    distance = np.insert(distance, 0, 0)/distance[-1]
    new_s = np.linspace(0, 1, (len(x)-1)*100)
    f = interp1d(distance, points, kind="linear", axis=0)
    x,y = f(new_s).T
    del f,points,distance
    x = fft(x)
    y = fft(y)
    low_pass = int(len(x)*0.001)
    x[low_pass:-low_pass] = 0
    y[low_pass:-low_pass] = 0
    fftd = np.arange(0,len(x))*2j*np.pi/len(x)
    pos = len(x)//2-1
    fftd[-pos:] = -fftd[1:pos+1][::-1]
    dx = -ifft(x*fftd).real
    dy = -ifft(y*fftd).real
    x = ifft(x).real
    y = ifft(y).real
    xmin = np.amin(x)
    xmax = np.amax(x)
    ymin = np.amin(y)
    ymax = np.amax(y)
    ratio_x = (max(np.concatenate(c[:,:,0]))-min(np.concatenate(c[:,:,0])))/(xmax-xmin)
    ratio_y = (max(np.concatenate(c[:,:,1]))-min(np.concatenate(c[:,:,1])))/(ymax-ymin)
    new_x = (x-xmin)*ratio_x + min(np.concatenate(c[:,:,0]))
    new_y = (y-ymin)*ratio_y + min(np.concatenate(c[:,:,1]))
    return np.abs(0.5*np.sum(new_y*dx - new_x*dy))

BLACKLIST = type, ModuleType, FunctionType
def getsize(obj):
    if isinstance(obj, BLACKLIST):
        raise TypeError('getsize() does not take argument of type: '+ str(type(obj)))
    seen_ids = set()
    size = 0
    objects = [obj]
    while objects:
        need_referents = []
        for obj in objects:
            if not isinstance(obj, BLACKLIST) and id(obj) not in seen_ids:
                seen_ids.add(id(obj))
                size += sys.getsizeof(obj)
                need_referents.append(obj)
        objects = get_referents(*need_referents)
    return size
