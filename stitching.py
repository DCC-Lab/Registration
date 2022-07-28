from PIL import ImageOps, Image
import matplotlib.pyplot as plt
import numpy as np
from skimage.registration import phase_cross_correlation
import tifffile as tiff
import scipy.signal

import filesManagement as fman
from imageTreatment import *
#from typing import *

class Stitching(ImageTreatment):
	def __init__(self, sourceDir:str, tileD:list, imageSize:list, vShift:list, hShift:list, intensityCorrection:bool = False):
		super().__init__(sourceDir=sourceDir)
		if intensityCorrection == True:
			# the directory we are interested in is actually the one with intensity correction images. 
			self.directory, self.files = self.correct_intensity_envelop()
		else:
			self.directory = self.sourceDir

		self.tileD = tileD
		self.vShift = vShift
		self.hShift = hShift
		self.imageSize = imageSize
		self.intensityCorrection = intensityCorrection


	def calculate_coordinates_firstImage(self, tile):
		"""
		The first image is the only one that needs to consider the vertical and the horizontal shifts to be positioned properly. After positioning this first image, all the other images can be positioned according to one single image. 
		Calculates the vertical (vShift) and horizontal (hShift) shifts with phase cross-correlation. 
		According to the shifts, calculates the coordinates of the top-left pixel of the first image.
		Returns the coordinates in [-x, y], since we mirror the images to stitch.  
		"""
		hShift = self.calculate_shift_PCC(index1=0, index2=1)
		vShift = self.calculate_shift_PCC(index1=0, index2=tileD[0]) 

		# if an x value is negative, it means the neighbouring image goes to the left, so the first image must be pushed to the right. 
		if hShift[0] < 0:
			x = (tileD[1] - 1) * abs(hShift[0])
		elif vShift[0] < 0:
			x = (tileD[1] - 1) * abs(vShift[0])
		# if both x values are negative, the first image must be at the right extremity.
		elif vShift[0] and hShift[0] < 0:
			x = tile.size[0] - self.imageSize[0]
		else:
			x = 0

		# if an y value is negative, it means the neighbouring image goes upwards, so the first image must be positioned downwards. 
		if hShift[1] < 0:
			y = (tileD[0] - 1) * abs(hShift[1])
		elif vShift[1] < 0: 
			y = (tileD[0] - 1) * abs(vShift[1])

		# if both y values are negative, the first image must be at the bottom extremity.
		elif vShift[1] and hShift[1] < 0:
			x = tile.size[1] - self.imageSize[1]
		else:
			y = 0

		return [-x, y]


	def calculate_shift_PCC(self, index1:int, index2:int) -> list:
		"""
		Input the indexes of two images in a set.
		Calculates the spatial shift between two images using the phase cross-correlation.
		Returns a list corresponding to the shift [width,height], where a positive height 
		corresponds to a shift to the bottom and a positive weight corresponds to a shift to the right.
		"""
		image1 = fman.read_file(filePath=self.directory + "/" + self.files[index1], imageType="numpy")
		image2 = fman.read_file(filePath=self.directory + "/" + self.files[index2], imageType="numpy")
		print(self.files[index1], self.files[index2])

		reverseShift, error, disphase = phase_cross_correlation(reference_image=image1, moving_image=image2)
		shift = [reverseShift[1], reverseShift[0]]
	
		return shift

	def calculate_shift_convolution(self, index1:int, index2:int) -> list:
		allImages = fman.list_name_of_files(directory=self.directory)

		image1 = fman.read_file(filePath=self.directory + "/" + allImages[index1], imageType="numpy")
		image2 = fman.read_file(filePath=self.directory + "/" + allImages[index2], imageType="numpy")

		shift = scipy.signal.fftconvolve(image1, image2[::-1,::-1], mode='same')
		itself = scipy.signal.fftconvolve(image1, image1[::-1,::-1], mode='same')

		maxPeakShift = np.unravel_index(np.argmax(shift), shift.shape)
		maxPeakItself = np.unravel_index(np.argmax(itself), itself.shape)

		return

	def create_tile_image(self):
		"""
		Calculates the size of the tile image. 
		Creates a 8-bit black PIL image of the size of the tile.  
		Returns a 8-bit black PIL image. 
		"""
		width = self.imageSize[0] + (abs(self.hShift[0]) * (self.tileD[0]-1)) + abs(self.vShift[0])
		height = self.imageSize[1] + (abs(self.vShift[1]) * (self.tileD[1]-1)) + abs(self.hShift[1])
		#print(f"BACKGROUND IMAGE WEIGHT {width} AND HEIGHT {height}")

		newImage = Image.new(mode="L", size=[width, height])
	
		return newImage
	
	def merge_images_sidebyside(self, index1:int, index2:int):
	    """
	    Input the indexes of two images in a set. 
	    Merge two images into one, displayed side by side.
	    Returns the merged image object.
	    """
	    image1 = fman.read_file(filePath=self.directory + "/" + self.files[index1], imageType="PIL", mirror=True)
	    image2 = fman.read_file(filePath=self.directory + "/" + self.files[index2], imageType="PIL", mirror=True)
	
	    (width1, height1) = image1.size
	    (width2, height2) = image2.size
	
	    result_width = width1 + width2
	    result_height = max(height1, height2)
	
	    result = Image.new('RGB', (result_width, result_height))
	    result.paste(im=image1, box=(0, 0))
	    result.paste(im=image2, box=(width1, 0))
	    
	    return result
	
	def stitching_scrapbooking_allImages(self):
		""" 
		Creates the background tile image of the right size. 
		For all images of the list of files : 
			Calculates the [x,y] position (pixels) where the top-left pixel of the image has to be pasted. 
			Opens the image in PIL. 
			Pastes the image on the background tile image. 
		Returns the tile image with all the images pasted on it. 
		"""
		tile = self.create_tile_image()

		firstImageCoordinates = self.calculate_coordinates_firstImage(tile=tile)
		image = fman.read_file(filePath=self.directory + "/" + self.files[0], imageType="PIL", mirror=True)
		tile.paste(image, firstImageCoordinates)

		position = [1,0] # [width,height]
		i = 1

		# keeps track of the coordinates of the first image of the row (vCoordinates) and of the previous image (hCoordinates)
		vCoordinates = firstImageCoordinates
		hCoordinates = firstImagesCoordinates

		while position[1] < self.tileD[1]: # colonnes, y
			while position[0] < self.tileD[0]: # rangées, x
				# if first image of the row, use the image on top to calculate the shift
				if position[0] == 0:
					shift = self.calculate_shift_PCC(index1=i-tileD[0], index2=i)
					coordinates = [a + b for a, b in zip(vCoordinates, shift)]
					hCoordinates = coordinates
					vCoordinates = coordinates
				# if not first image of the row, use the previous image to calcualte the shift
				else:
					shift = self.calculate_shift_PCC(index1=i-1, index2=i)
					coordinates = [a + b for a, b in zip(hCoordinates, shift)]
					hCoordinates = coordinates

				image = fman.read_file(filePath=self.directory + "/" + self.files[i], imageType="PIL", mirror=True)
				tile.paste(image, coordinates)
				print(f"shift and coordinates [x,y] of image : {shift} and {coordinates}")
	
				coordinates[0] += 1
				i += 1
			coordinates[1] += 1
			coordinates[0] = 0
	
		return tile
	
	















