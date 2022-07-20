import os
import fnmatch
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
from skimage.registration import phase_cross_correlation
import tifffile as tiff

def calculate_shift(image1, image2):
	shift, error, diffphase = phase_cross_correlation(image1, image2)
	print(f'Shift, Error, diffphase : {shift, error, diffphase}')

	return shift, error, diffphase

#def create_new_image(filepath, tileDimensions):
#	image = Image.open(filepath)

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
    """Merge two images into one, displayed side by side
    :param file1: path to first image file
    :param file2: path to second image file
    :return: the merged Image object
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
	Reads the .tif file and convert them in a np.array. 
	Returns the file as a np.array. 
	"""
	if imageType == "numpy":
		image = tiff.imread(filePath)
	if imageType == "PIL":
		image = Image.open(filePath)

	return image





