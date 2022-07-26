import os
import fnmatch
from PIL import Image
import numpy as np

def listNameOfFiles(directory: str, extension="tif") -> list:
	"""
	Fetch files name. 
	Does not consider .DS_Store files. 
	Returns a list of the names of the files. 
	"""
	foundFiles = []
	for file in os.listdir(directory):
		if fnmatch.fnmatch(file, f'*.{extension}'):
			if file == ".DS_Store":
				pass
			else:
				foundFiles.append(file)

	foundFiles.sort()
	return foundFiles

def read_file(filePath, imageType, mirror=False):
	"""
	Reads the .tif file and convert them in a np.array or PIL image. 
	If mirror == True, mirrors the PIL image. 
	Returns the image in the right format. 
	"""
	if imageType == "numpy":
		image = tiff.imread(filePath)
	if imageType == "PIL":
		image = Image.open(filePath)
	if mirror == True:
		image = ImageOps.mirror(image)

	return image