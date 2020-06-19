import os, time, cv2
os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = pow(2,40).__str__()
import numpy as np
from skimage.io import imread, imsave
import matplotlib.pyplot as plt
import GAKT_Pore as GT


if __name__ == '__main__':
	IMAGE_NAME = "C:/University SSD Drive/file/location/of/image/to_be_analysed.tif"
	img = cv2.imread(IMAGE_NAME,0)
	F_NAME = os.path.splitext(IMAGE_NAME)[0]
	SAVE_NAME = F_NAME+".csv"
	HIST_SAVE_NAME = F_NAME+"-hist.csv"
	SAVE_FIG_NAME = F_NAME+"-homogeneity.tif"
	print("Starting program.")
	Pore_analysis = GT.AnalysePores(img,threshold_value=100,scale=1,G=False,white_background=True)
	print("Starting contour FFT processing.")
	t_pore_parallel = time.time()
	Pore_analysis.process_parallel()
	t_cont = time.time()
	print("Finished contour processing.",t_cont-t_pore_parallel,"seconds.")
	print("Starting Homogeneity.")
	Pore_analysis.process_free_area_parallel()
	t_homo = time.time()
	print("Finished Homogeneity processing.",t_homo-t_cont,"seconds.")
	print("Starting Homogeneity map creation.")
	Pore_analysis.process_homogeneity_colour_map_parallel(draw_contours=True)
	t_homo_map = time.time()
	print("Finished Homogeneity map processing.",t_homo_map-t_homo,"seconds.")
	print("Starting radial processing.")
	Pore_analysis.process_radial_contour(radii_n=20)
	t_rad = time.time()
	print("Finished radial processing.",t_rad-t_homo_map,"seconds.")
	t_pore_parallel = time.time()-t_pore_parallel
	print("Total Parallel Time:",t_pore_parallel)

	#Pore_analysis.conts #Contours
	#Pore_analysis.l_rad #List of all the largest diameter through pore relating to contour indices (e.g. Value for position 0 is the largest diameter, relating to contour 0)
	#Pore_analysis.area_rad #Radius of circle is calculated assuming the area is that of a perfect circle (sqrt(area/pi))
	#Pore_analysis.circ #List of all the circularity relating to contour indices
	#Pore_analysis.waviness #List of all the waviness relating to contour indices
	#Pore_analysis.aspect_ratios #List of all the aspect ratios relating to contour indices
	#Pore_analysis.area #List of all the areas relating to contour indices
	#Pore_analysis.FFT_conts # FFT'd Contours
	#Pore_analysis.total_porosity #Total porosity
	#Pore_analysis.homogeneity_map #Map of contour space - each x,y coordinate links to a contour indice
	#Pore_analysis.homogeneity_zoom # How much is the map enlarged by compared to original image - by default it is set to 1.
	#Pore_analysis.homogeneity_areas #Free space area
	#Pore_analysis.radii_n #Number of rings
	#Pore_analysis.radius_centre #Centre of the rings
	#Pore_analysis.radius_max #Final radius
	#Pore_analysis.radius_contours #The contour indices in each ring. For example: Pore_analysis.area[Pore_analysis.radius_contours[0]] will give you the areas in the first ring.
	#Pore_analysis.radius_step_len #The length of each step in radii. For example: the second ring will be within 1*step_len and 2*step_len.
	#Pore_analysis.radius_steps #Array of lengths between each ring.
	#Pore_analysis.radius_porosity #Porosity for rings calculated from binary image.
	#Pore_analysis.radius_porosity_num #Number of black pixels used to calculate means
	#Pore_analysis.radius_porosity_num #Number of pixels used to calculate means
	#Pore_analysis.homogeneity_colour_map #A RGB colour map based on volume fraction of pore. Must be made first using process_homogeneity_colour_map().
	#Pore_analysis.homogeneity_colour_mapper #Matplotlib mapper that can be used to show the colourbar.
	#assert False, "Stop"
	combined_array = np.array([ 
		Pore_analysis.radius_steps,
		Pore_analysis.radius_porosity,
		np.array([np.mean(Pore_analysis.l_rad[rc]) for rc in Pore_analysis.radius_contours]),
		np.array([np.max(Pore_analysis.l_rad[rc]) for rc in Pore_analysis.radius_contours]),
		np.array([np.std(Pore_analysis.l_rad[rc]) for rc in Pore_analysis.radius_contours]),
		np.array([len(rc) for rc in Pore_analysis.radius_contours]),
		np.array([np.mean(Pore_analysis.area[rc]) for rc in Pore_analysis.radius_contours]),
		np.array([np.std(Pore_analysis.area[rc]) for rc in Pore_analysis.radius_contours]),
		np.array([np.mean(Pore_analysis.circ[rc]) for rc in Pore_analysis.radius_contours]),
		np.array([np.std(Pore_analysis.circ[rc]) for rc in Pore_analysis.radius_contours]),
		np.array([np.mean(Pore_analysis.waviness[rc]) for rc in Pore_analysis.radius_contours]),
		np.array([np.std(Pore_analysis.waviness[rc]) for rc in Pore_analysis.radius_contours]),
		np.array([np.mean(Pore_analysis.aspect_ratios[rc]) for rc in Pore_analysis.radius_contours]),
		np.array([np.std(Pore_analysis.aspect_ratios[rc]) for rc in Pore_analysis.radius_contours]),
		np.array([np.mean(Pore_analysis.area[rc]/Pore_analysis.homogeneity_areas[rc]) for rc in Pore_analysis.radius_contours]),
		np.array([np.std(Pore_analysis.area[rc]/Pore_analysis.homogeneity_areas[rc]) for rc in Pore_analysis.radius_contours]),
		np.array([np.mean(Pore_analysis.area[rc]/(2*Pore_analysis.l_rad[rc])**2) for rc in Pore_analysis.radius_contours]),
		np.array([np.std(Pore_analysis.area[rc]/(2*Pore_analysis.l_rad[rc])**2) for rc in Pore_analysis.radius_contours])
		]).T
	with open(SAVE_NAME,"w") as f_save:
		f_save.write("Step (um),Porosity (%),Mean Diameter (um),Max Diameter(um), Std Diameter (um), Number of samples, Mean Area (um^2), Std Area (um^2), Circularity, std Circ, Waviness, std Wav, Aspect Ratio, std Aspect Ratio, Mean Area fraction, std Area fraction, mean Elongation, std Elongation")
		for a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r in combined_array:
			f_save.write("\n"+str(a)+","+str(b)+","+str(c*2)+","+str(d*2)+","+str(e*2)+","+str(f)+","+str(g)+","+str(h)+","+str(i)+","+str(j)+","+str(k)+","+str(l)+","+str(m)+","+str(n)+","+str(o)+","+str(p)+","+str(q)+","+str(r))
	combined_array = np.array([Pore_analysis.l_rad*2,Pore_analysis.area_rad*2,Pore_analysis.circ,Pore_analysis.waviness,Pore_analysis.aspect_ratios,Pore_analysis.area,Pore_analysis.area/Pore_analysis.homogeneity_areas]).T
	with open(HIST_SAVE_NAME,"w") as f_save:
		f_save.write("LSTP(um), Circular diameter, Circularity, Waviness, Aspect Ratio, Pore Area, Local Volume fraction %, Elongation (based on LSTP)")
		for i, j, k, l, m, n, o in combined_array:
			f_save.write("\n"+str(i)+","+str(j)+","+str(k)+","+str(l)+","+str(m)+","+str(n)+","+str(o)+","+str(n/i**2))
	f = plt.figure(figsize=[6.4,4.8])
	plt.imshow(Pore_analysis.homogeneity_colour_map)
	plt.colorbar(Pore_analysis.homogeneity_colour_mapper)
	plt.axis('off')
	bbox = f.axes[1].get_window_extent().transformed(f.dpi_scale_trans.inverted())
	f.set_size_inches(4.8*img.shape[1]/img.shape[0]+bbox.width*2.5, 4.8)
	fig_bbox = f.axes[0].get_tightbbox(f.canvas.get_renderer())
	plt.savefig(SAVE_FIG_NAME,dpi=(img.shape[0]*f.dpi/fig_bbox.height),bbox_inches='tight')
	plt.show()

def disp_pore(n):
	plt.plot(*np.concatenate(Pore_analysis.conts[n]).T)
	plt.plot(*Pore_analysis.FFT_conts[n].T)
	plt.show()