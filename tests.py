import unittest
import os

import stitching as stitch
import filesManagement as fman

pathTestImages = os.getcwd() + "/tests/testDataset"
tileDim = [2,3] #[x,y]
# test images are already intensity adjusted. 
stitcherPCC = stitch.Stitching(sourceDir=pathTestImages, tileD=tileDim, imageSize=[1024,512], isIntensityCorrection=False, isMirrored=True, isFlipped=False)
stitcherConv = stitch.Stitching(sourceDir=pathTestImages, tileD=tileDim, imageSize=[1024,512], isIntensityCorrection=False, shiftEstimation="FFTConvolution", isMirrored=True, isFlipped=False)

class TestStitching(unittest.TestCase):
	def test_calculateShiftPCC_horizontalShift(self):
		with self.subTest("Horizontal shift estimation with PCC returns a list."):
			shift = stitcherPCC.calculate_shift_PCC(imageRef1=0, imageRef2=1)
			self.assertIsInstance(shift, list)

		with self.subTest("Horizontal shift estimation with PCC returns a list of two elements."):
			shift = stitcherPCC.calculate_shift_PCC(imageRef1=0, imageRef2=1)
			self.assertEqual(len(shift), 2)

		with self.subTest("Horizontal shift estimation with PCC returns a valid x-shift."):
			shift = stitcherPCC.calculate_shift_PCC(imageRef1=0, imageRef2=1)
			self.assertTrue(132 <= shift[0] <= 142)

		with self.subTest("Horizontal shift estimation with PCC returns a valid y-shift for horizontal stitch."):
			shift = stitcherPCC.calculate_shift_PCC(imageRef1=0, imageRef2=1)
			self.assertTrue(-18 <= shift[1] <= -8)

	def text_calculateShiftConvolution_horizontalShift(self):
		with self.subTest("Horizontal shift estimation with FFT convolution returns a list."):
			shift = stitcherConv.calculate_shift_convolution(imageRef1=0, imageRef2=1)
			self.assertIsInstance(shift, list)

		with self.subTest("Horizontal shift estimation with FFT convolution returns a list of two elements."):
			shift = stitcherConv.calculate_shift_convolution(imageRef1=0, imageRef2=1)
			self.assertEqual(len(shift), 2)

		with self.subTest("Horizontal shift estimation with FFT convolution returns valid x-shift."):
			shift = stitcherPCC.calculate_shift_convolution(imageRef1=0, imageRef2=1)
			self.assertTrue(132 <= shift[0] <= 142)

		with self.subTest("Horizontal shift estimation with FFT convolution returns valid y-shift."):
			shift = stitcherPCC.calculate_shift_convolution(imageRef1=0, imageRef2=1)
			self.assertTrue(-18 <= shift[1] <= -8)

	def test_createBlackImage(self):
		with self.subTest("Create black image with specified size (PCC for shift estimation)."):
			imagePCC = stitcherPCC.create_black_image(width=1000, height=500)
			self.assertEqual(imagePCC.size, (1000,500))

		with self.subTest("Create black image with specified size (FFT convolution for shift estimation)."):
			imageConv = stitcherConv.create_black_image(width=500, height=1000)
			self.assertEqual(imageConv.size, (500, 1000))

		with self.subTest("Create black image with no specified size (PCC for shift estimation)."):
			width = stitcherPCC.imageSize[0] + (abs(stitcherPCC.hShift[0]) * (stitcherPCC.tileD[0]-1)) + abs(stitcherPCC.vShift[0])
			height = stitcherPCC.imageSize[1] + (abs(stitcherPCC.vShift[1]) * (stitcherPCC.tileD[1]-1)) + abs(stitcherPCC.hShift[1])
			image = stitcherPCC.create_black_image()
			self.assertTrue(width <= image.size[0]) 
			self.assertTrue(height <= image.size[1])

		with self.subTest("Create black image with no specified size (FFT convolution for shift estimation)."):
			width = stitcherPCC.imageSize[0] + (abs(stitcherPCC.hShift[0]) * (stitcherPCC.tileD[0]-1)) + abs(stitcherPCC.vShift[0])
			height = stitcherPCC.imageSize[1] + (abs(stitcherPCC.vShift[1]) * (stitcherPCC.tileD[1]-1)) + abs(stitcherPCC.hShift[1])
			image = stitcherPCC.create_black_image()
			self.assertTrue(width <= image.size[0]) 
			self.assertTrue(height <= image.size[1])

	def test_estimateShift(self):
		with self.subTest("Vertical shift estimation with PCC returns a list."):
			shift = stitcherPCC.estimate_shift(index=0, stitchingSide="V", shiftMethod="PCC")
			self.assertIsInstance(shift, list)

		with self.subTest("Vertical shift estimation with FFT convolution returns a list."):
			shift = stitcherConv.estimate_shift(index=0, stitchingSide="V", shiftMethod="FFTConvolution")
			self.assertIsInstance(shift, list)

		with self.subTest("Vertical shift estimation with PCC returns a list of two elements."):
			shift = stitcherPCC.estimate_shift(index=0, stitchingSide="V", shiftMethod="PCC")
			self.assertEqual(len(shift), 2)

		with self.subTest("Vertical shift estimation with FFT convolution returns a list of two elements."):
			shift = stitcherConv.estimate_shift(index=0, stitchingSide="V", shiftMethod="FFTConvolution")
			self.assertEqual(len(shift), 2)

		with self.subTest("Vertical shift estimation with PCC returns a valid x-shift."):
			shift = stitcherPCC.estimate_shift(index=0, stitchingSide="V", shiftMethod="PCC")
			self.assertTrue(22 <= shift[0] <= 32)

		with self.subTest("Vertical shift estimation with FFT convolution returns a valid x-shift."):
			shift = stitcherConv.estimate_shift(index=0, stitchingSide="V", shiftMethod="FFTConvolution")
			self.assertTrue(22 <= shift[0] <= 32)

		with self.subTest("Vertical shift estimation with PCC returns a valid y-shift."):
			shift = stitcherPCC.estimate_shift(index=0, stitchingSide="V", shiftMethod="PCC")
			self.assertTrue(296 <= shift[1] <= 306)

		with self.subTest("Vertical shift estimation with FFT convolution returns a valid y-shift."):
			shift = stitcherConv.estimate_shift(index=0, stitchingSide="V", shiftMethod="FFTConvolution")
			self.assertTrue(296 <= shift[1] <= 306)

	def test_PCCAndFFTConvolutionShifts(self):
		with self.subTest("Vertical shift estimations with PCC and FFT convolution give approximately the same result for x-shift."):
			shiftPCC = stitcherPCC.estimate_shift(index=0, stitchingSide="V", shiftMethod="PCC")
			shiftConv = stitcherConv.estimate_shift(index=0, stitchingSide="V", shiftMethod="FFTConvolution")
			deltaX = abs(shiftPCC[0]-shiftConv[0])
			self.assertTrue(deltaX <= 10)

		with self.subTest("Vertical shift estimations with PCC and FFT convolution give approximately the same result for y-shift."):
			shiftPCC = stitcherPCC.estimate_shift(index=0, stitchingSide="V", shiftMethod="PCC")
			shiftConv = stitcherConv.estimate_shift(index=0, stitchingSide="V", shiftMethod="FFTConvolution")
			deltaY = abs(shiftPCC[1]-shiftConv[1])
			self.assertTrue(deltaY <= 10)

		with self.subTest("Horizontal shift estimations with PCC and FFT convolution give approximately the same result for x-shift."):
			shiftPCC = stitcherPCC.calculate_shift_PCC(imageRef1=0, imageRef2=1)
			shiftConv = stitcherConv.calculate_shift_convolution(imageRef1=0, imageRef2=1)
			deltaX = abs(shiftPCC[0]-shiftConv[0])
			self.assertTrue(deltaX <= 10)

		with self.subTest("Horizontal shift estimations with PCC and FFT convolution give approximately the same result for y-shift."):
			shiftPCC = stitcherPCC.calculate_shift_PCC(imageRef1=0, imageRef2=1)
			shiftConv = stitcherConv.calculate_shift_convolution(imageRef1=0, imageRef2=1)
			deltaY = abs(shiftPCC[1]-shiftConv[1])
			self.assertTrue(deltaY <= 10)

if __name__ == "__main__":
     unittest.main()
