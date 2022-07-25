import os
import fnmatch
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
from skimage.registration import phase_cross_correlation
import tifffile as tiff

def calculate_shift_PCC(image1, image2):
	"""
	Takes two in np.arrays and calculates the spacial shift of the second np.array according 
	to the first np.array with phase cross-correlation.
	Returns a list corresponding to the shift [weight,height], where a positive height 
	corresponds to a shift to the bottom and a positive weight corresponds to a shift to the 
	right.
	"""
	reverseShift, error, diffphase = phase_cross_correlation(image1, image2)
	#print(f'Shift, Error, diffphase : {shift, error, diffphase}')
	shift = [reverseShift[1], reverseShift[0]]

	return shift

def create_tile_image(tile:list, shift:list, imageSize=[1024,512]):
	"""
	Calculates the size of the tile image. 
	Creates a 8-bit black PIL image of the size of the tile.  
	Returns a 8-bit black PIL image. 
	"""
	weight = imageSize[0] + (abs(shift[0]) * (tile[0]-1))
	height = imageSize[1] + (abs(shift[1]) * (tile[1]-1))
	blackImageSize = [weight, height]
	
	newImage = Image.new(mode="L", size=blackImageSize)

	return newImage

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

def merge_images_sidebyside(image1, image2):
    """
    Merge two images into one, displayed side by side.
    Returns the merged image object.
    """

    (width1, height1) = image1.size
    (width2, height2) = image2.size

    result_width = width1 + width2
    result_height = max(height1, height2)

    result = Image.new('RGB', (result_width, result_height))
    result.paste(im=image1, box=(0, 0))
    result.paste(im=image2, box=(width1, 0))
    
    return result

def read_file(filePath, imageType):
	"""
	Reads the .tif file and convert them in a np.array or PIL image. 
	Returns the image in the right format. 
	"""
	if imageType == "numpy":
		image = tiff.imread(filePath)
	if imageType == "PIL":
		image = Image.open(filePath)

	return image





