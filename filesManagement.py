import os
import fnmatch
from PIL import Image, ImageOps
import numpy as np
import tifffile as tiff


def create_new_directory(directory: str, newFileName: str):
	"""
	Create new folder. 
	Returns the path of the new folder. 
	"""
	newDir = directory+"/"+newFileName
	if not os.path.exists(newDir):
		os.makedirs(newDir)
	return newDir
	
def list_name_of_files(directory: str, extension="tif") -> list:
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

def read_file(filePath, imageType, mirror=False, flip=False):
	"""
	Reads the .tif file and convert them in a np.array or PIL image. 
	If mirror == True, mirrors the PIL image. 
	Returns the image in the right format. 
	"""
	if imageType == "numpy":
		image = tiff.imread(filePath)
	if imageType == "PIL":
		image = Image.open(filePath)
	if mirror:
		image = ImageOps.mirror(image)
	if flip:
		image = ImageOps.flip(image)

	return image