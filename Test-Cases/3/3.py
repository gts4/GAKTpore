import GAKTpore, os

#The if-statement

if __name__ == "__main__":
	output_directory = "./output/"
	if not os.path.exists(output_directory):
		os.mkdir(output_directory)
	GAKTpore.analysis.run("./3-fft.tif",output_directory,80,200/260) #For the scale, 200um was measured to be 262 pixels.