from PIL import ImageOps, Image
import matplotlib.pyplot as plt
import numpy as np
from skimage.registration import phase_cross_correlation
import tifffile as tiff
import scipy.signal

import exceptions as exc
import filesManagement as fman
from imageTreatment import *
#from typing import *

class Stitching(ImageTreatment):
	def __init__(self, sourceDir:str, tileD:list, imageSize:list, isIntensityCorrection:bool = False, isMirrored:bool = False, isFlipped:bool = False):
		super().__init__(sourceDir=sourceDir)
		if isIntensityCorrection:
			# if true, the directory we are interested in is actually the one with intensity-corrected images. 
			self.directory, self.files = self.correct_intensity_envelop()
		else:
			self.directory = self.sourceDir

		self.tileD = tileD
		self.imageSize = imageSize
		self.isMirrored = isMirrored
		self.isFlipped = isFlipped

		# vertical (vShift) and horizontal (hShift) shifts between the first image and its neighbours. 
		self.hShift = self.calculate_shift_PCC(index1=0, index2=1)
		self.vShift = self.estimate_shift(index=0, stitchingSide="V")

	def calculate_coordinates_firstImage(self, image):
		"""
		The first image is the only one that needs to consider the vertical and the horizontal shifts to be positioned properly. After positioning this first image, all the other images can be positioned according to one single image. 
		Calculates the vertical (vShift) and horizontal (hShift) shifts with phase cross-correlation. 
		According to the shifts, calculates the coordinates of the top-left pixel of the first image.
		Verifies if the user wants to mirror the images. If so, the sign of the x coordinate changes. 
		Returns the coordinates in [x, y]. 
		"""

		# if an x value is negative, it means the neighbouring image goes to the left, so the first image must be pushed to the right. 
		if self.vShift[0] < 0 and self.hShift[0] < 0:
			x = image.size[0] - self.imageSize[0]
		elif self.hShift[0] < 0:
			x = (self.tileD[1] - 1) * abs(self.hShift[0])
		elif self.vShift[0] < 0:
			x = (self.tileD[1] - 1) * abs(self.vShift[0])
		else:
			x = 0

		# if an y value is negative, it means the neighbouring image goes upwards, so the first image must be positioned downwards. 
		if self.vShift[1] < 0 and self.hShift[1] < 0:
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

	def calculate_shift_PCC(self, imageRef1, imageRef2) -> list:
		"""
		Input the indexes of two images in a set.
		If index1 and index2 are int, opens the two images are numpy images. Else, it must be an image, pillow or numpy, that is input in the function. 
		Calculates the spatial shift between two images using the phase cross-correlation.
		Verifies if the user wants to mirror and/or flip the images. If so, the sign of the x and/or y coordinates change.
		Returns a list corresponding to the shift [width,height], where a positive height corresponds to a shift to the bottom and a positive weight corresponds to a shift to the right.
		"""
		if type(imageRef1) and type(imageRef2) == int:
			npimage1 = fman.read_file(filePath=self.directory + "/" + self.files[imageRef1], imageType="numpy")
			npimage2 = fman.read_file(filePath=self.directory + "/" + self.files[imageRef2], imageType="numpy")
		else: 
			npimage1 = np.asarray(imageRef1)
			npimage2 = np.asarray(imageRef2)

		reverseShift, error, disphase = phase_cross_correlation(reference_image=npimage1, moving_image=npimage2)
		
		# the sign of the x and/or y values of shift might need some change according to the flip or mirror state. 
		if self.isMirrored and type(imageRef1) == int:
			shift = [-int(reverseShift[1]), int(reverseShift[0])]
		elif self.isFlipped and type(imageRef1) == int:
			shift = [int(reverseShift[1]), -int(reverseShift[0])]
		elif self.isFlipped and self.isMirrored and type(imageRef1) == int:
			shift = [-int(reverseShift[1]), -int(reverseShift[0])]
		else:
			shift = [int(reverseShift[1]), int(reverseShift[0])]
	
		return shift

	def calculate_shift_convolution(self, index1:int, index2:int) -> list:
		"""
		NOT FINISHED!
		Input the indexes of two images in a set. 
		Generates an FFT convolution of two images. 
		Calculates the coordinates of the maximum peaks of the FFt convolution of the reference image with itself and of the reference image with the moving image. Ideally, this would be used to calculate the spatial shift between the two images. 
		Returns nothing for now, but would return the shift in [x,y].
		"""
		image1 = fman.read_file(filePath=self.directory + "/" + self.files[index1], imageType="numpy")
		image2 = fman.read_file(filePath=self.directory + "/" + self.files[index2], imageType="numpy")

		shift = scipy.signal.fftconvolve(image1, image2[::-1,::-1], mode='same')
		autocorr = scipy.signal.fftconvolve(image1, image1[::-1,::-1], mode='same')

		maxPeakShift = np.unravel_index(np.argmax(shift), shift.shape)
		maxPeakAutocorr = np.unravel_index(np.argmax(autocorr), autocorr.shape)

		print(f"MAX PEAKS : {maxPeakShift} and {maxPeakAutocorr}")

		return

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

	def estimate_shift(self, index:int, stitchingSide:str) -> list:
		"""
		If the shift is too high, using the whole image to estimate the shift (which whatever method) is biased. This function uses a cropped version of the reference and the moving images to help the shift estimation, if needed. 
		Input the index of the current image and the stitchingSide, which can be horizontal ("H") or vertical ("V").
		Defines the indexes of the reference and moving images. 
		Opens the images in pillow. 
		Calculates the region of both images to keep to estimate the shift. 
		Applies a low-pass filter on the crop images (to get rid of the noise). 
		Estimates the shift between the crop images. 
		Rescales the shift to make it correspond to the shift from the top-left coordinates. 
		Returns the shift [x,y]. 
		"""
		# Define the indexes. 
		if index == 0 and stitchingSide == "V": 
			referenceIndex = 0
			movingIndex = self.tileD[0]

		elif index != 0 and stitchingSide == "V":
			referenceIndex = index-self.tileD[0]
			movingIndex = index

		elif stitchingSide == "H": 
			referenceIndex = index
			movingIndex = index+1

		else:
			exc.define_variable(stitchingSide)
			exc.define_variable(index)

		reference = fman.read_file(filePath=self.directory + "/" + self.files[referenceIndex], imageType="PIL", mirror=self.isMirrored, flip=self.isFlipped)
		moving = fman.read_file(filePath=self.directory + "/" + self.files[movingIndex], imageType="PIL", mirror=self.isMirrored, flip=self.isFlipped)

		# Define the coordinates of the crop (which is the region you want to keep). 
		# I decided to use 250 rows of pixels for the vertical shift and 500 columns of pixels for the horizontal shift, but that can be changed by the user. 
		if stitchingSide == "V":
			rleft = 0
			rtop = self.imageSize[1] - 250
			mright = self.imageSize[0]
			mbottom = 250

		elif stitchingSide == "H":
			rleft = self.imageSize[0] - 500
			rtop = 0
			mright = 500
			mbottom = self.imageSize[1]

		else:
			exc.define_variable(stitchingSide)

		rright = self.imageSize[0]
		rbottom = self.imageSize[1]
		mleft = 0
		mtop = 0

		# Crop the reference and moving images and converts in np.asarray.
		cropReference = np.asarray(reference.crop((rleft, rtop, rright, rbottom)))
		cropMoving = np.asarray(moving.crop((mleft, mtop, mright, mbottom)))

		# Apply a low-pass filter on the cropped images. 
		lowCropReference = self.apply_low_pass_filter(image=cropReference)
		lowCropMoving = self.apply_low_pass_filter(image=cropMoving)

		# Estimates shift of low-passed cropped images.
		shift = self.calculate_shift_PCC(index1=lowCropReference, index2=lowCropMoving)

		# Rescales the shift to make it correspond to the top-left coordinate (makes up for the crop).
		if stitchingSide == "V":
			shift[1] = shift[1] + (self.imageSize[1]-250)

		if stitchingSide == "H":
			shift[0] = shift[0] + (self.imageSize[0]-500)

		return shift
	
	def stitching_scrapbooking_allImages(self):
		""" 
		Creates the background tile image of the right size. 
		Calculates the coordinates of the first image according to its two neighbouring images. 
		For all images of the list of files : 
			Calculates the [x,y] coordinates (pixels) where the top-left pixel of the image has to be pasted. 
			Opens the image in PIL. 
			Pastes the image on the background tile image. 
		Returns the tile image with all the images pasted on it. 
		"""
		tile = self.create_black_image()

		firstImageCoordinates = self.calculate_coordinates_firstImage(image=tile)
		image = fman.read_file(filePath=self.directory + "/" + self.files[0], imageType="PIL", mirror=self.isMirrored, flip=self.isFlipped)
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
					shift = self.estimate_shift(index=i, stitchingSide="V")
					coordinates = [vCoordinates[0] + shift[0], vCoordinates[1] + shift[1]]
					hCoordinates = [coordinates[0], coordinates[1]]
					vCoordinates = [coordinates[0], coordinates[1]]
				# if not first image of the row, use the previous image to calcualte the shift
				else:
					shift = self.calculate_shift_PCC(index1=i-1, index2=i)
					print(f"shift last : {shift}")
					coordinates = [hCoordinates[0] + shift[0], hCoordinates[1] + shift[1]]
					hCoordinates = [coordinates[0], coordinates[1]]

				image = fman.read_file(filePath=self.directory + "/" + self.files[i], imageType="PIL", mirror=self.isMirrored, flip=self.isFlipped)
				print(coordinates)
				tile.paste(image, (coordinates[0], coordinates[1]))
	
				position[0] += 1
				i += 1
			position[1] += 1
			position[0] = 0
	
		return tile
	
	















