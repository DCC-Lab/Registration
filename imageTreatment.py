import filesManagement as fman
import numpy as np
import scipy.ndimage as simg
import tifffile as tiff
import math as math

class ImageTreatment:
	def __init__(self, sourceDir:str):
		self.sourceDir = sourceDir

		self.files = fman.list_name_of_files(directory=sourceDir)

	def adjust_intensity(self, image, correction):
		""" 
		Takes an image and multiplies each pixels by the corresponding pixel in the correction image. 
		Returns the image corrected in intensity. 
		"""
		image = image.astype(np.float64)

		x = 0
		y = 0
		while x < image.shape[0]:
			while y < image.shape[1]: 
				image[x][y] = image[x][y] * correction[x][y]
				y += 1
			y = 0
			x += 1

		new8bitImage = np.uint8(self.rescale_image(image))

		# adjust the background for better contrast. Here the background is defined as any pixel under 20. 
		x = 0
		y = 0
		while x < new8bitImage.shape[0]:
			while y < new8bitImage.shape[1]:
				if new8bitImage[x][y] < 20:
					new8bitImage[x][y] = 0
				y += 1
			y = 0
			x += 1

		return new8bitImage

	def correct_intensity_envelop(self):
		"""
		Creates an intensity-correction image. 
		Adjusts the intensity of images, one at a time. 
		Saves the images in a new directory. 
		Returns the path of the new directory and the list of the corrected images. 
		"""
		aveImage = self.create_average_image()
		correctionImage = self.create_intensity_correction_image(image=aveImage)

		pathAfterCorrection = fman.create_new_directory(directory=self.sourceDir, newFileName="IntensityCorrection")
		
		for file in self.files:
			image = fman.read_file(filePath=self.sourceDir + "/" + file, imageType="numpy")
			correctedImage = self.adjust_intensity(image=image, correction=correctionImage)
			newFileName = pathAfterCorrection + "/" + "ADJ" + file 
			tiff.imwrite(newFileName, correctedImage)

		correctedFiles = fman.list_name_of_files(directory=pathAfterCorrection)

		return pathAfterCorrection, correctedFiles


	def create_average_image(self):
		"""
		With the image directory, averages them all to produce a final image for correction. 
		Returns the resultant average image. 
		"""
		pixels = self.sum_pixels()
		numberOfImages = len(self.files)

		x = 0
		y = 0
		while x < pixels.shape[0]:
			while y < pixels.shape[1]:
				pixels[x][y] = pixels[x][y]/numberOfImages 
				y += 1
			y = 0
			x += 1

		return pixels

	def create_intensity_correction_image(self, image):
		"""
		From an average image of all images in a set, generates the intensity correction image that should be applied on individual images. 
		Applies a gaussian blur on the correction image to smoothen the stuff. 
		Returns the intensity correction image with values between 0 and 1.
		"""
		inverseImage = self.inverse_pixels(image=image)
		rescaledImage = self.rescale_image(image=inverseImage)
		correction = simg.gaussian_filter(rescaledImage, sigma=10, mode="nearest")

		return correction

	def inverse_pixels(self, image):
		"""
		Takes the input image and inverses the pixels (0s become 1s and 1s become 0s). 
		"""
		inverseImage = np.zeros(shape=(image.shape[0], image.shape[1]))
		maxPixel = np.amax(image)
		x = 0
		y = 0
		while x < image.shape[0]:
			while y < image.shape[1]:
				inverseImage[x][y] = maxPixel - image[x][y]
				y += 1
			y = 0
			x +=1

		return inverseImage

	def low_Pass_Filter(self, image, sigmaFilter):
		"""
		Creates a low-pass filter in the Fourier space. 
		As an imput, 
			- the image must be a numpy array;
			- signamFilter defines the width of the filter.
		Returns the Fourier-image of the filter. This should be multiplied by the FFT image to filter, then iFFT the product to recover the image.
		"""
		x, y = np.meshgrid(np.linspace(-1,1,image.shape[1]), np.linspace(-1,1,image.shape[0]))
		d = np.sqrt(x*x+y*y)
		sigma = (sigmaFilter*0.18)/(2*math.sqrt(2*math.log(2)))
		mu = 0.0
		gauss = (1/(sigma*2*np.pi)) * np.exp(-((d-mu)**2/(2.0*sigma**2)))
		maxPixel = np.amax(gauss)
		gaussNorm = gauss/maxPixel

		return gaussNorm

	def normalize_image(self, image):
		""" 
		Creates an image with pixels in float 64 varying form 0 to 1 according to the intensity of each pixels from the input image.
		"""
		newImage = np.zeros(shape=(image.shape[0], image.shape[1]))
		maxPixel = np.amax(image)

		x = 0
		y = 0
		while x < image.shape[0]:
			while y < image.shape[1]:
				newImage[x][y] = image[x][y]/maxPixel
				y += 1
			y = 0
			x += 1

		return newImage

	def rescale_image(self, image, pixelRange=255):
		"""
		Normalizes the image with the function normalizeImage. 
		Rescales the image on the range defined with the variable pixelRange. 
		"""
		normImage = self.normalize_image(image=image)
		rescaledImage = np.zeros(shape=(image.shape[0], image.shape[1]))
		
		x = 0
		y = 0
		while x < image.shape[0]:
			while y < image.shape[1]:
				rescaledImage[x][y] = normImage[x][y] * pixelRange
				y += 1
			y = 0
			x += 1

		return rescaledImage

	def subtract_value_on_all_pixels(self, value, image):
		"""
		Takes a numpy image and a value to subtract in input. 
		For each element in the image that are higher than the value to be subtracted, the value is subtracted from the value. 
		Returns the subtracted image. 
		"""
		x = 0
		y = 0
		while y < image.shape[0]:
			while x < image.shape[1]:
				if image[y, x] > value:
					image[y, x] = image[y, x] - value
				x += 1
			x = 0
			y += 1

		return image

	def sum_pixels(self):
		"""
		Reads all images in directory and sums the value of each pixels. 
		Returns an image with the sum of all pixels. 
		"""
		firstPath = self.sourceDir + "/" + self.files[0]
		firstImage = fman.read_file(filePath=firstPath, imageType="numpy")

		pixels = np.zeros(shape=(firstImage.shape[0], firstImage.shape[1]), dtype=np.float64)

		i = 0
		for name in self.files:
			#filePath = self.sourceDir + "/" + name
			image = fman.read_file(filePath=self.sourceDir + "/" + name, imageType="numpy")
			x = 0
			y = 0
			while x < firstImage.shape[0]:
				while y < firstImage.shape[1]:
					pixels[x][y] = pixels[x][y] + image[x][y]
					y += 1
				y = 0
				x += 1
			i += 1

		return pixels




