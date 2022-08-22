from PIL import ImageOps, Image
import matplotlib.pyplot as plt
import numpy as np
from skimage.registration import phase_cross_correlation
import tifffile as tiff
from scipy.signal import fftconvolve
from tqdm import tqdm

import filesManagement as fman
from imageTreatment import *
#from typing import *

class Stitching(ImageTreatment):
	def __init__(self, sourceDir:str, tileD:list, imageSize:list, isIntensityCorrection:bool = True, shift = None, isMirrored:bool = False, isFlipped:bool = False):
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
		if shift is None:
			shiftMethod = input("Which stitching method do you want to use (Answer 1, 2, 3 or 4): \n 1. Estimate the shift with phase cross-correlation. \n 2. Estimate the shift with FFT convolution \n 3. Use the position in the file name. \n 4. Use a predefined shift.")
			if shiftMethod == "1":
				self.hShift = self.calculate_shift_PCC(imageRef1=0, imageRef2=1)
				self.vShift = self.estimate_shift(index=0, stitchingSide="V", shiftMethod="PCC")
				image, averageH, averageV = self.stitch_with_estimated_shift(shiftMethod="PCC")
				self.save_image(image)
			elif shiftMethod == "2":
				self.hShift = self.calculate_shift_convolution(imageRef1=0, imageRef2=1)
				self.vShift = self.estimate_shift(index=0, stitchingSide="V", shiftMethod="FFTConvolution")
				image, averageH, averageV = self.stitch_with_estimated_shift(shiftMethod="FFTConvolution")
				self.save_image(image)
			elif shiftMethod == "3":
				self.hShift = self.calculate_shift_from_file_name(imageRef1=0, imageRef2=1)
				#self.vShift = self.calculate_shift_from_file_name(imageRef1=0, imageRef2=self.tileD[0])
				self.vShift = [0,0]
				image = self.stitch_with_position_in_file_name()
				self.save_image(image)
			elif shiftMethod == "4":
				hx, hy = input("Please enter the x and y values of the HORIZONTAL shift (x, y):")
				vx, vy = input("Please enter the x and y values of the VERTICAL shift (x, y):")
				self.hShift = [hx, hy]
				self.vShift = [vx, vy]
				image = self.stitch_with_known_shifts()
				self.save_image(image)

		else:
			if shift == "PCC":
				self.hShift = self.calculate_shift_PCC(imageRef1=0, imageRef2=1)
				self.vShift = self.estimate_shift(index=0, stitchingSide="V", shiftMethod="PCC")
				image, averageH, averageV = self.stitch_with_estimated_shift(shiftMethod="PCC")
				self.save_image(image)
			elif shiftMethod == "FFTConvolution":
				self.hShift = self.calculate_shift_convolution(imageRef1=0, imageRef2=1)
				self.vShift = self.estimate_shift(index=0, stitchingSide="V", shiftMethod="FFTConvolution")
				image, averageH, averageV = self.stitch_with_estimated_shift(shiftMethod="FFTConvolution")
				self.save_image(image)
			elif shiftMethod == "FileName":
				self.hShift = self.calculate_shift_from_file_name(imageRef1=0, imageRef2=1)
				self.vShift = self.calculate_shift_from_file_name(imageRef1=0, imageRef2=self.tileD[0])
				image = self.stitch_with_position_in_file_name()
				self.save_image(image)
			elif shiftMethod == "Manual":
				hx, hy = input("Please enter the x and y values of the HORIZONTAL shift (x, y):")
				vx, vy = input("Please enter the x and y values of the VERTICAL shift (x, y):")
				self.hShift = [hx, hy]
				self.vShift = [vx, vy]
				image = self.stitch_with_known_shifts()
				self.save_image(image)


	def average_shifts(self, allShifts):
		"""
		Input a list of shifts [x,y] and returns the average of all the shifts element-wise. 
		"""
		length = len(allShifts)
		allShifts = np.dstack(tuple(allShifts))
		sumAllShifts = np.sum(allShifts, axis=2)
		average = sumAllShifts/length
		average = [int(average[0][0]), int(average[0][1])]

		return average

	def calculate_coordinates_firstImage(self, background):
		"""
		The first image is the only one that needs to consider the vertical and the horizontal shifts to be positioned properly. After positioning this first image, all the other images can be positioned according to one single image. 
		Calculates the vertical (vShift) and horizontal (hShift) shifts with phase cross-correlation. 
		According to the shifts, calculates the coordinates of the top-left pixel of the first image.
		Verifies if the user wants to mirror the images. If so, the sign of the x coordinate changes. 
		Returns the coordinates in [x, y]. 
		"""
		print(f"h and v : {self.hShift} and {self.vShift}")

		# if an x value is negative, it means the neighbouring image goes to the left, so the first image must be pushed to the right. 
		if self.vShift[0] < 0 and self.hShift[0] < 0:
			x = background.size[0] - self.imageSize[0]
		elif self.hShift[0] < 0:
			x = (self.tileD[1] - 1) * abs(self.hShift[0])
		elif self.vShift[0] < 0:
			x = (self.tileD[1] - 1) * abs(self.vShift[0])
		else:
			x = 0

		# if an y value is negative, it means the neighbouring image goes upwards, so the first image must be positioned downwards. 
		if self.vShift[1] < 0 and self.hShift[1] < 0:
			y = background.size[1] - self.imageSize[1]
		elif self.hShift[1] < 0:
			y = (self.tileD[0] - 1) * abs(self.hShift[1])
		elif self.vShift[1] < 0: 
			y = (self.tileD[0] - 1) * abs(self.vShift[1])
		else:
			y = 0

		if self.isMirrored == True:
			x = -x

		print(f"First image coordinates done.")

		return [int(x), int(y)]

	def calculate_shift_from_file_name(self, imageRef1, imageRef2):
		"""
		Finds the position of the image at indexes imageRef1 and imageRef2 from their file name. 
		Indeed, when images are acquired with Nirvana, the file name includes the position of the image. 
		By extracting this position in x and y, it would be possible to know the exact shift between two images. 
		Ex : 
			file name 1 is ADJCORR20220620-334B-1A-Cryoprotectant-Ave30-tiling-zoom2-004-Bliq VMS-Ch-2-Ch 2-VOI_1-Slice_180-X001-Y180-Z001-(X11465.31-Y5069.28-Z4374.00)-Time_238.163s-Frame_007260-Projection_Average.tif
			The index 112 is the parenthesis where the information related to the image position starts (X11465.31-Y5069.28-Z4374.00).
			file name 2 is ADJCORR20220620-334B-1A-Cryoprotectant-Ave30-tiling-zoom2-004-Bliq VMS-Ch-2-Ch 2-VOI_1-Slice_181-X001-Y181-Z001-(X11465.31-Y5094.52-Z4374.00)-Time_239.347s-Frame_007297-Projection_Average.tif
			The index 112 is the parenthesis where the information related to the image position starts (X11465.31-Y5094.52-Z4374.00).
			Use this information to know the shift between image 1 and 2 (0, -25.24). 
		We know that the X- and Y-axis are inversed, so X is actually Y and vice versa. 
		"""
		npimage1 = fman.read_file(filePath=self.directory + "/" + self.files[imageRef1], imageType="numpy")
		npimage2 = fman.read_file(filePath=self.directory + "/" + self.files[imageRef2], imageType="numpy")

		nameRef = self.files[imageRef1]
		nameMoving = self.files[imageRef2]

		stringIndexRef = nameRef.find("(X") # should be 112 with ADJ
		yRef = ""
		i = stringIndexRef + 2
		isAPosition = True
		while isAPosition:
			if nameRef[i].isnumeric() or nameRef[i] == ".":
				yRef += nameRef[i]
				i += 1
			else:
				isAPosition = False
		
		xRef = ""
		i += 2
		isAPosition = True
		while isAPosition:
			if nameRef[i].isnumeric() or nameRef[i] == ".":
				xRef += nameRef[i]
				i += 1
			else:
				isAPosition = False

		positionRef = [float(xRef), float(yRef)]

		stringIndexMoving = nameMoving.find("(X") # should be 112 with ADJ
		yMoving = ""
		i = stringIndexMoving + 2
		isAPosition = True
		while isAPosition:
			if nameMoving[i].isnumeric() or nameMoving[i] == ".":
				yMoving += nameMoving[i]
				i += 1
			else:
				isAPosition = False
		
		xMoving = ""
		i += 2
		isAPosition = True
		while isAPosition:
			if nameMoving[i].isnumeric() or nameMoving[i] == ".":
				xMoving += nameMoving[i]
				i += 1
			else:
				isAPosition = False

		positionMoving = [float(xMoving), float(yMoving)]

		if self.isMirrored and self.isFlipped:
			shift = [-int(positionRef[0]-positionMoving[0]), -int(positionRef[1]-positionMoving[1])]
		elif self.isMirrored:
			shift = [-int(positionRef[0]-positionMoving[0]), int(positionRef[1]-positionMoving[1])]
		else:
			shift = [int(positionRef[0]-positionMoving[0]), int(positionRef[1]-positionMoving[1])]

		if shift[0] > 30:
			raise Exception(f"x-shift weird here. {shift[0]}")
		elif shift[1] > 60:
			raise Exception(f"y-shift werid here. {shift[1]}")

		return shift


	def calculate_shift_PCC(self, imageRef1, imageRef2) -> list:
		"""
		Input the indexes of two images in a set.
		If indices are input, according images are fetched from an image list. Else, it must be an image, pillow or numpy, that is input in the function. 
		Calculates the spatial shift between two images using the phase cross-correlation.
		Applies mirroring or flipping to the result if prompted to.
		Returns a list corresponding to the shift [width,height], where a positive height corresponds to a shift to the bottom and a positive weight corresponds to a shift to the right.
		"""
		if type(imageRef1) and type(imageRef2) == int:
			indexGiven = True
			npimage1 = fman.read_file(filePath=self.directory + "/" + self.files[imageRef1], imageType="numpy")
			npimage2 = fman.read_file(filePath=self.directory + "/" + self.files[imageRef2], imageType="numpy")
		else: 
			indexGiven = False
			npimage1 = np.asarray(imageRef1)
			npimage2 = np.asarray(imageRef2)

		reverseShift, error, disphase = phase_cross_correlation(reference_image=npimage1, moving_image=npimage2)
		
		# the sign of the x and/or y values of shift might need some change according to the flip or mirror state. 
		if indexGiven:
			if self.isMirrored:
				reverseShift[1] *= -1
			if self.isFlipped:
				reverseShift[0] *= -1

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

		shift = list(map(lambda i, j: j - i, maxPeakAutocorr, maxPeakShift))

		if index_given:
			if self.isMirrored:
				shift[1] *= -1
			if self.isFlipped:
				shift[0] *= -1

		shift = [int(shift[1]), int(shift[0])]
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

	def estimate_shift(self, index:int, stitchingSide:str, shiftMethod:str="PCC") -> list:
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
			referenceIndex = index-1
			movingIndex = index

		else:
			raise Exception("Verify that all arguments are defined to estimate shift.")

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
			rleft = self.imageSize[0] - 900
			rtop = 0
			mright = 900
			mbottom = self.imageSize[1]

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
		if shiftMethod == "PCC":
			shift = self.calculate_shift_PCC(imageRef1=lowCropReference, imageRef2=lowCropMoving)
		elif shiftMethod == "FFTConvolution":
			shift = self.calculate_shift_convolution(imageRef1=sobelLowCropReference, imageRef2=sobelLowCropMoving)

		# Rescales the shift to make it correspond to the top-left coordinate (makes up for the crop).
		if stitchingSide == "V":
			shift[1] = shift[1] + (self.imageSize[1]-250)

		if stitchingSide == "H":
			shift[0] = shift[0] + (self.imageSize[0]-900)

		return shift

	def save_image(self, image):
		a = input("Want to save? (Answer y/n)")
		if a == "y":
			name = input("How do you want to name the stitched image?")
			path = "/Users/valeriepineaunoel/Desktop/" + name + ".tiff"
			image.save(fp=path, format="TIFF", subsampling=0, quality=100)
			print(fr"Saved as {path}")


	def stitch_with_position_in_file_name(self): 
		"""
		Stitch images according to a known shift between images.
		Creates the background tile image of the right size. 
		Pastes the first image at the top-left corner. 
		For all images of the list of files : 
			Calculates the [x,y] coordinates (pixels) according to a known shift where the top-left pixel of the image has to be pasted. 
			Opens the image in PIL. 
			Pastes the image on the background tile image. 
		Returns the tile image with all the images pasted on it.
		"""

		tile = self.create_black_image()

		i = 0

		print("Stitching images with positions in file name.")
		for y in tqdm(range(self.tileD[1])): # rangées, y
			for x in range(self.tileD[0]): # colonnes, x
				# if first image of the row, use the image on top to calculate the shift
				if x == 0:
					if y == 0:
						coordinates = [0,0]						
						vCoordinates = coordinates
						hCoordinates = coordinates
					else:
						shift = self.calculate_shift_from_file_name(imageRef1=i, imageRef2=i-self.tileD[0])
						coordinates = [vCoordinates[0] + shift[0], vCoordinates[1] + shift[1]]
						hCoordinates = coordinates
						vCoordinates = coordinates
				# if not first image of the row, use the previous image to calcualte the shift
				else:
					shift = self.calculate_shift_from_file_name(imageRef1=i-1, imageRef2=i)
					coordinates = [hCoordinates[0] + shift[0], hCoordinates[1] + shift[1]]
					hCoordinates = coordinates

				image = fman.read_file(filePath=self.directory + "/" + self.files[i], imageType="PIL", mirror=self.isMirrored, flip=self.isFlipped)
				coords = coordinates
				tile.paste(image, (coordinates[0], coordinates[1]))
	
				i += 1
	
		return tile


	def stitch_with_known_shifts(self):
		"""
		Stitch images according to a known shift between images.
		Creates the background tile image of the right size. 
		Pastes the first image at the top-left corner. 
		For all images of the list of files : 
			Calculates the [x,y] coordinates (pixels) according to a known shift where the top-left pixel of the image has to be pasted. 
			Opens the image in PIL. 
			Pastes the image on the background tile image. 
		Returns the tile image with all the images pasted on it.
		"""

		tile = self.create_black_image()

		i = 0

		print("Stitching images with known shift.")
		for y in tqdm(range(self.tileD[1])): # rangées, y
			for x in range(self.tileD[0]): # colonnes, x
				# if first image of the row, use the image on top to calculate the shift
				if x == 0:
					if y == 0:
						coordinates = self.calculate_coordinates_firstImage(background=tile)
						vCoordinates = coordinates
						hCoordinates = coordinates
					else:
						shift = self.vShift
						coordinates = [vCoordinates[0] + shift[0], vCoordinates[1] + shift[1]]
						hCoordinates = coordinates
						vCoordinates = coordinates
				# if not first image of the row, use the previous image to calcualte the shift
				else:
					shift = self.hShift
					coordinates = [hCoordinates[0] + shift[0], hCoordinates[1] + shift[1]]
					hCoordinates = coordinates

				image = fman.read_file(filePath=self.directory + "/" + self.files[i], imageType="PIL", mirror=self.isMirrored, flip=self.isFlipped)
				coords = coordinates
				tile.paste(image, (coordinates[0], coordinates[1]))
	
				i += 1
	
		return tile


	def stitch_with_estimated_shift(self, shiftMethod):
		""" 
		Stitch images by estimating the shift for each image. 
		Creates the background tile image of the right size. 
		Estimates the coordinates of the first image according to its two neighbouring images. 
		For all images of the list of files : 
			Calculates the [x,y] coordinates (pixels) according to the estimated shift done by the method of choice where the top-left pixel of the image has to be pasted. 
			Opens the image in PIL. 
			Pastes the image on the background tile image. 
		Calculates the average vertical and horizontal shifts. 
		Returns the tile image with all the images pasted on it, the mean horizontal shift and the mean vertical shift. 
		"""
		tile = self.create_black_image()

		i = 0
		allVerticalShifts = []
		allHorizontalShifts = []

		print("Stitching images with shift estimation.")
		for y in tqdm(range(self.tileD[1])): # rangées, y
			for x in range(self.tileD[0]): # colonnes, x
				# if first image of the row, use the image on top to calculate the shift
				if x == 0:
					if y == 0:
						coordinates = self.calculate_coordinates_firstImage(background=tile)
						vCoordinates = coordinates
						hCoordinates = coordinates
					else:
						shift = self.estimate_shift(index=i, stitchingSide="V", shiftMethod=shiftMethod)
						allVerticalShifts.append(shift)

						coordinates = [vCoordinates[0] + shift[0], vCoordinates[1] + shift[1]]
						hCoordinates = coordinates
						vCoordinates = coordinates
				# if not first image of the row, use the previous image to calcualte the shift
				else:
					shift = self.estimate_shift(index=i, stitchingSide="H", shiftMethod=shiftMethod)
					allHorizontalShifts.append(shift)

					coordinates = [hCoordinates[0] + shift[0], hCoordinates[1] + shift[1]]
					hCoordinates = coordinates

				image = fman.read_file(filePath=self.directory + "/" + self.files[i], imageType="PIL", mirror=self.isMirrored, flip=self.isFlipped)
				coords = coordinates
				tile.paste(image, (coordinates[0], coordinates[1]))
	
				i += 1

		averageHShifts = self.average_shifts(allHorizontalShifts)
		averageVShifts = self.average_shifts(allVerticalShifts)
	
		return tile, averageHShifts, averageVShifts
	
	















