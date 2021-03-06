import os, time, cv2
os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = pow(2,40).__str__()
import numpy as np
import matplotlib.pyplot as plt
import GAKTpore as GT

def run(IMAGE_NAME:str, SAVE_FOLDER:str, THRES:int, SCALE:float, UPSCALE_MULTIPLIER:int = 1, W_BG:bool = True, FFT:bool = True, parallel=True, cpu_count=-1, npy_save:bool = False, plt_save:bool = False, draw_contours:bool=True, vmin:float=0, vmax:float=1, radii_n:int=10):
	img = cv2.imread(IMAGE_NAME,0)
	if UPSCALE_MULTIPLIER != 1:
		img = cv2.resize(img,(int(img.shape[1]*UPSCALE_MULTIPLIER), int(img.shape[0]*UPSCALE_MULTIPLIER)), interpolation=cv2.INTER_CUBIC)
	F_NAME = os.path.join(SAVE_FOLDER,os.path.splitext(os.path.split(IMAGE_NAME)[1])[0])
	if FFT:
		SAVE_NAME = F_NAME+"-"+str(THRES)+"-"+str(UPSCALE_MULTIPLIER)+"x_FFT.csv"
		HIST_SAVE_NAME = F_NAME+"-"+str(THRES)+"-"+str(UPSCALE_MULTIPLIER)+"x-hist_FFT.csv"
		SAVE_FIG_NAME = F_NAME+"-"+str(THRES)+"-"+str(UPSCALE_MULTIPLIER)+"x-homogeneity_FFT.png"
	else:
		SAVE_NAME = F_NAME+"-"+str(THRES)+"-"+str(UPSCALE_MULTIPLIER)+"x_savgol.csv"
		HIST_SAVE_NAME = F_NAME+"-"+str(THRES)+"-"+str(UPSCALE_MULTIPLIER)+"x-hist_savgol.csv"
		SAVE_FIG_NAME = F_NAME+"-"+str(THRES)+"-"+str(UPSCALE_MULTIPLIER)+"x-homogeneity_savgol.png"
	print("Starting program. Processing: "+os.path.basename(IMAGE_NAME))
	Pore_analysis = GT.AnalysePores(img,cpu_count=cpu_count,threshold_value=THRES,scale=SCALE/UPSCALE_MULTIPLIER,G=False,white_background=W_BG)
	if FFT:
		print("Starting contour FFT processing.")
	else:
		print("Starting contour Savgol processing.")
	if parallel:
		t_pore_parallel = time.time()
		Pore_analysis.process_parallel(FFT=FFT)
		t_cont = time.time()
		print("Finished contour processing.",t_cont-t_pore_parallel,"seconds.")
		print("Starting Homogeneity.")
		Pore_analysis.process_free_area_parallel()
		t_homo = time.time()
		print("Finished Homogeneity processing.",t_homo-t_cont,"seconds.")
		print("Starting Homogeneity map creation.")
		Pore_analysis.process_homogeneity_colour_map_parallel(draw_contours=draw_contours,vmax=vmax,vmin=vmin)
		t_homo_map = time.time()
		print("Finished Homogeneity map processing.",t_homo_map-t_homo,"seconds.")
		print("Starting radial processing.")
		Pore_analysis.process_radial_contour(radii_n=radii_n)
		t_rad = time.time()
		print("Finished radial processing.",t_rad-t_homo_map,"seconds.")
		t_pore_parallel = time.time()-t_pore_parallel
		print("Total Parallel Time:",t_pore_parallel)
	else:
		t_pore = time.time()
		Pore_analysis.process(FFT=FFT)
		t_cont = time.time()
		print("Finished contour processing.",t_cont-t_pore,"seconds.")
		print("Starting Homogeneity.")
		Pore_analysis.process_free_area()
		t_homo = time.time()
		print("Finished Homogeneity processing.",t_homo-t_cont,"seconds.")
		print("Starting Homogeneity map creation.")
		Pore_analysis.process_homogeneity_colour_map(draw_contours=draw_contours,vmax=vmax,vmin=vmin)
		t_homo_map = time.time()
		print("Finished Homogeneity map processing.",t_homo_map-t_homo,"seconds.")
		print("Starting radial processing.")
		Pore_analysis.process_radial_contour(radii_n=radii_n)
		t_rad = time.time()
		print("Finished radial processing.",t_rad-t_homo_map,"seconds.")
		t_pore_parallel = time.time()-t_pore_parallel
		print("Total Time:",t_pore)
	if npy_save:
		if FFT:
			with open(F_NAME+"_FFT_homogeneity.npy",'wb') as f:
				np.save(f,Pore_analysis.homogeneity_map)
		else:
			with open(F_NAME+"_savgol_homogeneity.npy",'wb') as f:
				np.save(f,Pore_analysis.homogeneity_map)
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
		f_save.write("LSTP(um), Circular diameter, Circularity, Waviness, Aspect Ratio, Pore Area, Local Area fraction, Elongation (based on LSTP)")
		for i, j, k, l, m, n, o in combined_array:
			f_save.write("\n"+str(i)+","+str(j)+","+str(k)+","+str(l)+","+str(m)+","+str(n)+","+str(o)+","+str(n/(np.pi*(i/2)**2))) # corrected elongation on 29/12/2020
	if plt_save:
		f = plt.figure(figsize=[6.4,4.8])
		plt.imshow(Pore_analysis.homogeneity_colour_map)
		plt.colorbar(Pore_analysis.homogeneity_colour_mapper)
		plt.axis('off')
		bbox = f.axes[1].get_window_extent().transformed(f.dpi_scale_trans.inverted())
		f.set_size_inches(4.8*img.shape[1]/img.shape[0]+bbox.width*2.5, 4.8)
		fig_bbox = f.axes[0].get_tightbbox(f.canvas.get_renderer())
		plt.savefig(SAVE_FIG_NAME,dpi=(img.shape[0]*f.dpi/fig_bbox.height),bbox_inches='tight')
		plt.close(f)
	else:
		a = np.array([Pore_analysis.homogeneity_colour_mapper.get_clim()])
		f = plt.figure(figsize=(0.58,3.77))
		i = plt.imshow(a, cmap=Pore_analysis.homogeneity_colour_mapper.cmap)
		plt.gca().set_visible(False)
		f.tight_layout(pad=0)
		cax = plt.axes([0.1, 0.02, 0.35, 0.93])
		plt.colorbar(cax=cax)
		new_dpi = img.shape[0]/3.77
		f.dpi = new_dpi
		f.canvas.draw()
		max_len = max([len(i.get_text()) for i in cax.get_yticklabels()])
		f.set_size_inches((0.58+(max_len-1)*0.12,3.77))
		f.canvas.draw()
		data = np.array(f.canvas.renderer.buffer_rgba())[:,:,0:3]
		plt.close('all')
		cmap_shape = Pore_analysis.homogeneity_colour_map.shape
		img_2 = np.ones((cmap_shape[0],cmap_shape[1]+data.shape[1],cmap_shape[2]),dtype=np.uint8)*255
		img_2[:cmap_shape[0],:cmap_shape[1],:] = Pore_analysis.homogeneity_colour_map
		img_2[:min(data.shape[0],img_2.shape[0]),-data.shape[1]:,:] = data[:min(data.shape[0],img_2.shape[0]),:,:]
		cv2.imwrite(SAVE_FIG_NAME,img_2[:, :, ::-1])