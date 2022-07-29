from stitching import *

sourcePath = "/Users/valeriepineaunoel/Desktop/test20%"
tileDimensions = [2, 3]
# [x,y]. Positive y = going down. Positive x = going right. When comparing image2 with image1 : 
verticalShift = [47, 363]
horizontalShift = [163, -12]

stitch = Stitching(sourceDir=sourcePath, tileD=tileDimensions, imageSize=[1024,512], isIntensityCorrection=True, isMirrored=True)
stitch.calculate_shift_convolution(index1=0, index2=2)
stitchedImage = stitch.stitching_scrapbooking_allImages()
stitchedImage.save(fp="/Users/valeriepineaunoel/Desktop/stitchedImage.tiff")