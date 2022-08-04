import unittest
import os

import stitching as stitch

pathTestImages = os.getcwd() + "/tests/testDataset"
tileDim = [2,3] #[x,y]
stitcherPCC = stitch.Stitching(sourceDir=pathTestImages, tileD=tileDim, imageSize=[1024,512], isIntensityCorrection=False, isMirrored=True, isFlipped=False)
stitcherConvolution = stitch.Stitching(sourceDir=pathTestImages, tileD=tileDim, imageSize=[1024,512], isIntensityCorrection=False, shiftEstimation="FFTConvolution", isMirrored=True, isFlipped=False)

class TestInputArguments(unittest.TestCase):
	# not sure how this is relevant because it might change everytime the user wants to use another method to estimate shift. 
	#def test_horizontal_and_vertical_shift_estimations_first_image(self):
	#	with self.subTest("Estimates horizontal and vertical shifts of first image."):
	#		assertEqual(self.stitcher.hShift, [130, 231])
	#		assertEqual(self.stitcher.vShift, [27, 301])

	def test_create_black_image_with_specified_size(self):
		with self.subTest("Create black image with specified size [px]."):
			image = stitcherPCC.create_black_image(width=1000, height=500)
			self.assertEqual(image.size, (1000,500))

	def test_create_black_image(self):
		with self.subTest("Create black image with no specified size [px]."):
			width = stitcherPCC.imageSize[0] + (abs(stitcherPCC.hShift[0]) * (stitcherPCC.tileD[0]-1)) + abs(stitcherPCC.vShift[0])
			height = stitcherPCC.imageSize[1] + (abs(stitcherPCC.vShift[1]) * (stitcherPCC.tileD[1]-1)) + abs(stitcherPCC.hShift[1])
			image = stitcherPCC.create_black_image()
			assert width <= image.size[0] 
			assert height <= image.size[1]


if __name__ == "__main__":
     unittest.main()
