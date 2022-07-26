from PIL import ImageOps, Image
import matplotlib.pyplot as plt
from skimage.registration import phase_cross_correlation
import tifffile as tiff
import filesManagement as fman

class Stitching:
	def __init__(self, sourceDir:str, tileD:list, imageSize:list, vShift:list, hShift:list):
		self.sourceDir = sourceDir
		self.tileD = tileD
		self.vShift = vShift
		self.hShift = hShift
		self.imageSize = imageSize

		self.files = fman.listNameOfFiles(directory=sourceDir)


	def calculate_shift_PCC(self, index1:int, index2:int) -> list:
		"""
		Input the indexes of two images in a set.
		Calculates the spatial shift between two images using the phase cross-correlation.
		Returns a list corresponding to the shift [width,height], where a positive height 
		corresponds to a shift to the bottom and a positive weight corresponds to a shift to the 
		right.
		"""
		image1 = fman.read_file(filePath=self.sourceDir + "/" + self.files[index1], imageType="numpy")
		image2 = fman.read_file(filePath=self.sourceDir + "/" + self.files[index2], imageType="numpy")

		reverseShift, error, diffphase = phase_cross_correlation(image1, image2)
		#print(f'Shift, Error, diffphase : {shift, error, diffphase}')
		shift = [reverseShift[1], reverseShift[0]]
	
		return shift

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
	    image1 = fman.read_file(filePath=self.sourceDir + "/" + self.files[index1], imageType="PIL", mirror=True)
	    image2 = fman.read_file(filePath=self.sourceDir + "/" + self.files[index2], imageType="PIL", mirror=True)
	
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

		i = 0
		coordinates = [0,0] # [width,height]
	
		while coordinates[1] < self.tileD[1]: # colonnes, y
			row = 1
			coordinates[0] = 0
			while coordinates[0] < self.tileD[0]: # rangées, x
				x = (coordinates[1]*self.vShift[0]) + (coordinates[0]*self.hShift[0])
				y = ((self.tileD[0]-row)*abs(self.hShift[1])) + (coordinates[1]*self.vShift[1])
				print(f"coordinates [x,y] of image : {x} and {y}")
	
				image = fman.read_file(filePath=self.sourceDir + "/" + self.files[i], imageType="PIL", mirror=True)
				tile.paste(image, (x,y))
				coordinates[0] += 1
				row += 1
				i += 1
			coordinates[1] += 1
	
		return tile
	
	















