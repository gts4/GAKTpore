# User Guide to set up a basic Python environment for use with GAKTpore

The purpose of this guide is to aid in the set up of Python for the basic use of GAKTpore on Windows.

Python can be downloaded and installed from https://www.python.org/downloads/. Make sure to install the 64-bit version.

![Python 3.9 Installation Dialog](https://github.com/gts4/GAKTpore/raw/master/GAKTpore%20user-guide%20images/Python-install.png)

Note: Be sure to tick "Add Python to PATH". This ensures that Python can be used from the Command Prompt or Powershell without specifying the full path to python.exe each time. 

Press "Install now" to install Python.

We recommend a lightweight editor called [Visual Studio Code](https://code.visualstudio.com/download) for running code.

Once it is installed, open the program and click File -> Open Folder. Navigate to a folder of your choice, meant to store the code for running GAKTpore.

For this example, a folder named "GAKTpore testing folder" is used.
![GAKTpore testing folder](https://github.com/gts4/GAKTpore/raw/master/GAKTpore%20user-guide%20images/GAKTpore-testing-folder.png)

Click the new file button, located on the right side of the left panel, on the same row as the name of the chosen folder. 

For the name of the new file, use 'main.py' (without the quotation marks).

After which, please click on the file on the left panel to display it on right panel and then copy the following code into the file:

```Python
import GAKTpore, os

if __name__ == "__main__":
  image_path = "./1-fft.tif"
  output_directory = "./output/"
  if not os.path.exists(output_directory):
    os.mkdir(output_directory)
  GAKTpore.analysis.run(image_path,output_directory,85,200/260) 
```

Click File->Save or use the shortcut CTRL+S to save the file.

Before continuing, let us analyse the code:

* The first line will import the packages GAKTpore and os. os will allow us to check if folders exist and create folders (as well as many other file related functions).
* The 3rd line presents an if statement to check whether it is being run by the main instance of python. The technicalities are complex but this line is necessary to make use of parallel processing.
* The 4th line defines the path of to your image file. The './' selects the current folder ("GAKTpore testing folder" in this example) and the name of the image file being used is [1-fft.tif](https://github.com/gts4/GAKTpore/raw/master/Test-Cases/1/1-fft.tif). This is the first example from the provided [Test Cases](https://github.com/gts4/GAKTpore/tree/master/Test-Cases).
* The 5th line defines the path to your output directory. For this example, it just set to put them into a folder named "output" located in the current directory.
* The 6th and 7th line check whether the output folder provided exists, and if not it will create the folder instead.
* The 8th line will run GAKTpore.

![GAKTpore parameter numbers](https://github.com/gts4/GAKTpore/raw/master/GAKTpore%20user-guide%20images/GAKTpore-main-py.png)

The image above conveys the positions of the first four parameters. Their purpose is listed below:
1. The path to the image designated for analysis. This is defined on the 4th line.
2. The path to the output folder. This is defined on the 5th line.
3. The binary threshold - This number should be selected very carefully! Use a different piece of software such as ImageJ or GIMP to find the number to use for binarisation.
4. The scale - This is distance divided by pixel length. For the example of 1-fft.png, the 200um scalebar was measured to be 260 pixels long. Therefore the scale is 200/260.

Once these are configured, you are ready to run GAKTpore. But first, you must install GAKTpore package.

Click Terminal -> New Terminal. A new panel should appear on the bottom (called the terminal), allowing you to type text into it. 

Copy and paste the text below into the terminal.

```powershell
pip install GAKTpore
```

![GAKTpore Installation](https://github.com/gts4/GAKTpore/raw/master/GAKTpore%20user-guide%20images/GAKTpore-install.png)

This will install GAKTpore onto your system. After it is downloaded, you are ready to run your code.

*Note: Before running the code, you should have downloaded [1-fft.tif](https://github.com/gts4/GAKTpore/raw/master/Test-Cases/1/1-fft.tif) and placed it in your chosen directory*

Run the following code in the terminal.

```powershell
python ./main.py
```

![Running GAKTpore](https://github.com/gts4/GAKTpore/raw/master/GAKTpore%20user-guide%20images/GAKTpore-running.png)

You should have a new folder called output, with 3 files inside them. 2 csv files and 1 png file. 

The csv files can be analysed in a spreadsheet software like Microsoft Excel.

The png can be viewed in any image viewer. 

You may have noticed that the output image does not have good contrast. This is because the highest AF (Area Fraction) value on the colour bar is 1 by default. The maximum measured AF is no-where near 1. 
The file ending in "hist_FFT.csv" has all the information about the pores and their AF, and it can be shown that the maximum AF is 0.411.
By adding a 5th parameter to the function that runs GAKTpore (on line 8) as such:
```Python
  GAKTpore.analysis.run(image_path,output_directory,85,200/260,vmax=0.411) 
```
*Note: The "vmax=" determines the maximum number on the colour bar. A parameter with "vmin=" will set the minimum number on the colour bar (0 by default).*

By changing the line and rerunning the code, you can see a much better improvement in the contrast.

## Tip for selecting a different image file path

When selecting a file from a different folder, you can click on the 'copy path' button on the top selection of Windows Explorer as displayed below using a red arrow:

![Copy Path](https://github.com/gts4/GAKTpore/raw/master/GAKTpore%20user-guide%20images/Copy-path.png)

However when you replace the path into your main.py file, make sure you have an r before the quotation marks with your directory path inside.
*Example path to a different image file*

```Python
  image_path = r"D:\Python\Different Folder\2-fft.tif"
```

The r is **super important**! Python normally interprets the "\\" as an [escape sequence](https://en.wikipedia.org/wiki/Escape_sequence) and so your path will not be interpreted correctly. 
To get around this, you can either have r in front of the quotation marks (as above) or replace each instance of "\\" with "\\\" or "/".
