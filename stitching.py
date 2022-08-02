from PIL import ImageOps, Image
import matplotlib.pyplot as plt
import numpy as np
from skimage.registration import phase_cross_correlation
import tifffile as tiff
from scipy.signal import fftconvolve

import filesManagement as fman
from imageTreatment import *
#from typing import *

class Stitching(ImageTreatment):
	def __init__(self, sourceDir:str, tileD:list, imageSize:list, isIntensityCorrection:bool = False, isMirrored:bool = False):
		super().__init__(sourceDir=sourceDir)
		if isIntensityCorrection == True:
			# if true, the directory we are interested in is actually the one with intensity-corrected images. 
			self.directory, self.files = self.correct_intensity_envelop()
		else:
			self.directory = self.sourceDir

		self.tileD = tileD
		self.imageSize = imageSize
		self.isMirrored = isMirrored

		# vertical (vShift) and horizontal (hShift)n shifts between the first image and its neighbours. 
		#self.hShift = self.calculate_shift_PCC(index1=0, index2=1)
		#self.vShift = self.estimate_vertical_shift(index=0)

	def calculate_coordinates_firstImage(self, image):
		"""
		The first image is the only one that needs to consider the vertical and the horizontal shifts to be positioned properly. After positioning this first image, all the other images can be positioned according to one single image. 
		Calculates the vertical (vShift) and horizontal (hShift) shifts with phase cross-correlation. 
		According to the shifts, calculates the coordinates of the top-left pixel of the first image.
		Verifies if the user wants to mirror the images. If so, the sign of the x coordinate changes. 
		Returns the coordinates in [x, y]. 
		"""

		print(f"HSHIFT : {self.hShift}")
		# if an x value is negative, it means the neighbouring image goes to the left, so the first image must be pushed to the right. 
		if self.vShift[0] and self.hShift[0] < 0:
			x = image.size[0] - self.imageSize[0]
		elif self.hShift[0] < 0:
			x = (self.tileD[1] - 1) * abs(self.hShift[0])
		elif self.vShift[0] < 0:
			x = (self.tileD[1] - 1) * abs(self.vShift[0])
		else:
			x = 0

		# if an y value is negative, it means the neighbouring image goes upwards, so the first image must be positioned downwards. 
		if self.vShift[1] and self.hShift[1] < 0:
			y = image.size[1] - self.imageSize[1]
		elif self.hShift[1] < 0:
			y = (self.tileD[0] - 1) * abs(self.hShift[1])
		elif self.vShift[1] < 0: 
			y = (self.tileD[0] - 1) * abs(self.vShift[1])
		else:
			y = 0

		if self.isMirrored == True:
			x = -x

		return [int(x), int(y)]

	def calculate_shift_PCC(self, index1, index2) -> list:
		"""
		Input the indexes of two images in a set.
		If index1 and index2 are int, opens the two images are numpy images. Else, it must be an image, pillow or numpy, that is input in the function. 
		Calculates the spatial shift between two images using the phase cross-correlation.
		Verifies if the user wants to mirror the images. If so, the sign of the x coordinate changes.
		Returns a list corresponding to the shift [width,height], where a positive height corresponds to a shift to the bottom and a positive weight corresponds to a shift to the right.
		"""
		if type(index1) and type(index2) == int:
			npimage1 = fman.read_file(filePath=self.directory + "/" + self.files[index1], imageType="numpy")
			npimage2 = fman.read_file(filePath=self.directory + "/" + self.files[index2], imageType="numpy")
		else: 
			npimage1 = np.asarray(index1)
			npimage2 = np.asarray(index2)

		reverseShift, error, disphase = phase_cross_correlation(reference_image=npimage1, moving_image=npimage2)
		
		if self.isMirrored == True and type(index1) == int:
			shift = [-int(reverseShift[1]), int(reverseShift[0])]
		else:
			shift = [int(reverseShift[1]), int(reverseShift[0])]
	
		return shift

	def calculate_shift_convolution(self, imageRef1:int, imageRef2:int) -> list:
		"""
		Input two images, either as PIL images or as indices.
		If indices are input, according images are fetched from an image list.
		Does a convolution between both images and finds the position of the maximum value of the convoluted image.
		Does a convolution between the first image and itself and finds the position of the maximum value.
		Finds the difference between these two maxima, resulting in a pixel shift.
		Applies mirroring or flipping to the result if prompted to.
		Returns the shift as [x, y]
		"""
		if type(imageRef1) and type(imageRef2) == int:
			index_given = True
			npimage1 = fman.read_file(filePath=self.directory + "/" + self.files[imageRef1], imageType="numpy")
			npimage2 = fman.read_file(filePath=self.directory + "/" + self.files[imageRef2], imageType="numpy")
		else:
			index_given = False
			npimage1 = np.asarray(imageRef1)
			npimage2 = np.asarray(imageRef2)

		relative_shift = fftconvolve(npimage1, npimage2[::-1, ::-1], mode='same')
		autocorr = fftconvolve(npimage1, npimage1[::-1, ::-1], mode='same')

		maxPeakShift = np.unravel_index(np.argmax(relative_shift), relative_shift.shape)
		maxPeakAutocorr = np.unravel_index(np.argmax(autocorr), autocorr.shape)

		shift = list(map(lambda i, j: i - j, maxPeakAutocorr, maxPeakShift))

		if index_given:
			if self.isMirrored:
				shift[1] *= -1
			if self.isFlipped:
				shift[0] *= -1

		shift = [shift[1], shift[0]]
		return shift

	def create_black_image(self, width=None, height=None):
		"""
		Calculates the size of the tile image.
		If width and height are None, calculate the size of the image according to hShift and vShift. Else, any black image can be created with an input at width and height in pixels.  
		Creates a 8-bit black PIL image of the size of the tile.  
		Returns a 8-bit black PIL image. 
		"""
		if width is None and height is None: 
			width = self.imageSize[0] + (abs(self.hShift[0]) * (self.tileD[0]-1)) + abs(self.vShift[0]) + 100
			height = self.imageSize[1] + (abs(self.vShift[1]) * (self.tileD[1]-1)) + abs(self.hShift[1]) + 100

		newImage = Image.new(mode="L", size=[width, height])
	
		return newImage

	def estimate_vertical_shift(self, index:int):
		if index == 0: 
			referenceIndex = 0
			movingIndex = self.tileD[0]
		else : 
			referenceIndex = index-self.tileD[0]
			movingIndex = index

		noise1 = np.asarray(Image.effect_noise((self.imageSize[0], self.imageSize[1]-300), 12))
		subNoise1 = self.subtract_value_on_all_pixels(value=100, image=noise1)
		noiseImage1 = Image.fromarray(subNoise1)

		noise2 = np.asarray(Image.effect_noise((self.imageSize[0], self.imageSize[1]-100), 50))
		#subNoise2 = self.subtract_value_on_all_pixels(value=, image=noise2)
		noiseImage2 = Image.fromarray(noise2)

		reference = fman.read_file(filePath=self.directory + "/" + self.files[referenceIndex], imageType="PIL", mirror=self.isMirrored)
		reference.paste(noiseImage1, (0, 0))
		path = "/Users/valeriepineaunoel/Desktop/noiseImage1-" + str(referenceIndex) + ".tiff"
		reference.save(fp=path)

		moving = fman.read_file(filePath=self.directory + "/" + self.files[movingIndex], imageType="PIL", mirror=self.isMirrored)
		moving.paste(noiseImage2, (0, 100))
		path = "/Users/valeriepineaunoel/Desktop/noiseImage2-" + str(movingIndex) + ".tiff"
		moving.save(fp=path)

		shift = self.calculate_shift_PCC(index1=reference, index2=moving)
		print(f"shift vertical : {shift}")

		return shift


	def merge_images_sidebyside(self, index1:int, index2:int):
	    """
	    Input the indexes of two images in a set. 
	    Merge two images into one, displayed side by side.
	    Returns the merged image object.
	    """
	    image1 = fman.read_file(filePath=self.directory + "/" + self.files[index1], imageType="PIL", mirror=self.isMirrored)
	    image2 = fman.read_file(filePath=self.directory + "/" + self.files[index2], imageType="PIL", mirror=self.isMirrored)
	
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
			Calculates the [x,y] coordinates (pixels) where the top-left pixel of the image has to be pasted. 
			Opens the image in PIL. 
			Pastes the image on the background tile image. 
		Returns the tile image with all the images pasted on it. 
		"""
		tile = self.create_black_image()

		firstImageCoordinates = self.calculate_coordinates_firstImage(image=tile)
		image = fman.read_file(filePath=self.directory + "/" + self.files[0], imageType="PIL", mirror=self.isMirrored)
		tile.paste(image, (firstImageCoordinates[0], firstImageCoordinates[1]))
		print(f"first image : {firstImageCoordinates}")

		position = [1,0] # [width,height]
		i = 1

		# keeps track of the coordinates of the first image of the row (vCoordinates) and of the previous image (hCoordinates)
		vCoordinates = [firstImageCoordinates[0], firstImageCoordinates[1]]
		hCoordinates = [firstImageCoordinates[0], firstImageCoordinates[1]]

		while position[1] < self.tileD[1]: # colonnes, y
			while position[0] < self.tileD[0]: # rangées, x
				# if first image of the row, use the image on top to calculate the shift
				if position[0] == 0:
					shift = self.estimate_vertical_shift(index=i)
					coordinates = [vCoordinates[0] + shift[0], vCoordinates[1] + shift[1]]
					hCoordinates = [coordinates[0], coordinates[1]]
					vCoordinates = [coordinates[0], coordinates[1]]
				# if not first image of the row, use the previous image to calcualte the shift
				else:
					shift = self.calculate_shift_PCC(index1=i-1, index2=i)
					print(f"shift last : {shift}")
					coordinates = [hCoordinates[0] + shift[0], hCoordinates[1] + shift[1]]
					hCoordinates = [coordinates[0], coordinates[1]]

				image = fman.read_file(filePath=self.directory + "/" + self.files[i], imageType="PIL", mirror=self.isMirrored)
				print(coordinates)
				tile.paste(image, (coordinates[0], coordinates[1]))
	
				position[0] += 1
				i += 1
			position[1] += 1
			position[0] = 0
	
		return tile
	
	















