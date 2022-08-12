import filesManagement as fman
import numpy as np
import scipy.ndimage as simg
import tifffile as tiff
import math as math
import matplotlib.pyplot as plt

class ImageTreatment:
	def __init__(self, sourceDir:str):
		self.sourceDir = sourceDir

		self.files = fman.list_name_of_files(directory=sourceDir)

	def adjust_intensity(self, image, correction, threshold:int=20):
		""" 
		Takes an image and multiplies each pixels by the corresponding pixel in the correction image. 
		Returns the image corrected in intensity. 
		"""
		image = image.astype(np.float64)
		image *= correction
		new8bitImage = np.uint8(self.rescale_image(image))

		# adjust the background for better contrast. Here the background is defined as any pixel under 20.
		thresholdIndices = new8bitImage < threshold
		new8bitImage[thresholdIndices] = 0 

		return new8bitImage

	def apply_low_pass_filter(self, image):
		"""
		Applies a low-pass filter on the image given in input. 
		Returns the filtered image. 
		"""
		# FFT of np.asarray
		fftimage = np.fft.fftshift(np.fft.fft2(image))

		# Create low-pass filter. 
		lowPassFilter = self.low_Pass_Filter(image=image, sigmaFilter=1)

		# Apply low-pass filters on respective cropped images. 
		lowFFTimage = lowPassFilter * fftimage

		# iFFT of low-passed cropped images. 
		lowImage = np.fft.ifft2(np.fft.ifftshift(lowFFTimage))

		return lowImage

	def apply_sobel_filter(self, image):
		newImage = simg.sobel(image, mode="nearest")
		#tiff.imshow(newImage)
		#plt.show()

		return newImage

	def correct_intensity_envelop(self):
		"""
		Creates an intensity-correction image. 
		Adjusts the intensity of images, one at a time. 
		Saves the images in a new directory. 
		Returns the path of the new directory and the list of the corrected images. 
		"""
		print(f"Start to correct intensity on raw images.")
		aveImage = self.create_average_image()
		print(f"Average image is done.")
		correctionImage = self.create_intensity_correction_image(image=aveImage)
		print(f"Correction image is done.")

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
		pixels = pixels/numberOfImages

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
		maxPixel = np.amax(image)
		image = maxPixel - image

		return image

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

	def merge_images_sidebyside(self, index1:int, index2:int):
	    """
	    Input the indexes of two images in a set. 
	    Merge two images into one, displayed side by side.
	    Returns the merged image object.
	    """
	    image1 = fman.read_file(filePath=self.directory + "/" + self.files[index1], imageType="PIL", mirror=self.isMirrored, flip=self.isFlipped)
	    image2 = fman.read_file(filePath=self.directory + "/" + self.files[index2], imageType="PIL", mirror=self.isMirrored, flip=self.isFlipped)
	
	    (width1, height1) = image1.size
	    (width2, height2) = image2.size
	
	    result_width = width1 + width2
	    result_height = max(height1, height2)
	
	    result = Image.new('RGB', (result_width, result_height))
	    result.paste(im=image1, box=(0, 0))
	    result.paste(im=image2, box=(width1, 0))
	    
	    return result

	def normalize_image(self, image):
		""" 
		Creates an image with pixels in float 64 varying form 0 to 1 according to the intensity of each pixels from the input image.
		"""
		maxPixel = np.amax(image)
		image /= maxPixel

		return image

	def rescale_image(self, image, pixelRange=255):
		"""
		Normalizes the image with the function normalizeImage. 
		Rescales the image on the range defined with the variable pixelRange. 
		"""
		normImage = self.normalize_image(image=image)
		normImage *= pixelRange

		return normImage

	def subtract_value_on_all_pixels(self, value, image):
		"""
		Takes a numpy image and a value to subtract in input. 
		For each element in the image that are higher than the value to be subtracted, the value is subtracted from the value. 
		Returns the subtracted image. 
		"""
		image -= value

		return image

	def sum_pixels(self):
		"""
		Reads all images in directory and sums the value of each pixels. 
		Returns an image with the sum of all pixels. 
		"""
		allImages = []
		for name in self.files:
			image = fman.read_file(filePath=self.sourceDir + "/" + name, imageType="numpy")
			allImages.append(image)	
		allImages = np.dstack(tuple(allImages))
		sumImage = np.sum(allImages, axis=2)

		return sumImage




