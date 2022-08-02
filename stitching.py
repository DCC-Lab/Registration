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

		print(f"HSHIFT : {self.hShift}")
		print(f"VSHIFT : {self.vShift}")
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

	def estimate_shift(self, index:int, stitchingSide:str):
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

		reference = fman.read_file(filePath=self.directory + "/" + self.files[referenceIndex], imageType="PIL", mirror=self.isMirrored)
		moving = fman.read_file(filePath=self.directory + "/" + self.files[movingIndex], imageType="PIL", mirror=self.isMirrored)

		# crop defines the region you want to keep 
		if stitchingSide == "V":
			rleft = 0
			rtop = self.imageSize[1] - 250
			rright = self.imageSize[0]
			rbottom = self.imageSize[1]

			mleft = 0
			mtop = 0
			mright = self.imageSize[0]
			mbottom = 250

			print("VERTICAL")
		elif stitchingSide == "H":
			rleft = self.imageSize[0] - 500
			rtop = 0
			rright = self.imageSize[0]
			rbottom = self.imageSize[1]

			mleft = 0
			mtop = 0
			mright = 500
			mbottom = self.imageSize[1]

			print("HORIZONTAL")
		else:
			exc.define_variable(stitchingSide)

		cropReference = np.asarray(reference.crop((rleft, rtop, rright, rbottom)))
		cropMoving = np.asarray(moving.crop((mleft, mtop, mright, mbottom)))

		fftCropReference = np.fft.fftshift(np.fft.fft2(cropReference))
		fftCropMoving = np.fft.fftshift(np.fft.fft2(cropMoving))

		lowPassFilterReference = self.low_Pass_Filter(image=cropReference, sigmaFilter=1)
		lowPassFilterMoving = self.low_Pass_Filter(image=cropMoving, sigmaFilter=1)

		lowFFTCropReference = lowPassFilterReference * fftCropReference
		lowFFTCropMoving = lowPassFilterMoving * fftCropMoving

		lowCropReference = np.fft.ifft2(np.fft.ifftshift(lowFFTCropReference))
		lowCropMoving = np.fft.ifft2(np.fft.ifftshift(lowFFTCropMoving))

		shift = self.calculate_shift_PCC(index1=lowCropReference, index2=lowCropMoving)

		if stitchingSide == "V":
			shift[1] = shift[1] + (self.imageSize[1]-250)

		if stitchingSide == "H":
			shift[0] = shift[0] + (self.imageSize[0]-500)

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

				image = fman.read_file(filePath=self.directory + "/" + self.files[i], imageType="PIL", mirror=self.isMirrored)
				print(coordinates)
				tile.paste(image, (coordinates[0], coordinates[1]))
	
				position[0] += 1
				i += 1
			position[1] += 1
			position[0] = 0
	
		return tile
	
	















