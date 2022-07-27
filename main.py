from stitching import *

sourcePath = "/Users/valeriepineaunoel/Desktop/test20%"
tileDimensions = [2, 3]
# [x,y]. Positive y = going down. Positive x = going right. When comparing image2 with image1 : 
verticalShift = [47, 363]
horizontalShift = [163, -12]

print(sourcePath)

stitch = Stitching(sourceDir=sourcePath, tileD=tileDimensions, imageSize=[1024,512], vShift=verticalShift, hShift=horizontalShift)
stitchedImage = stitch.stitching_scrapbooking_allImages(correction=True)
stitchedImage.save(fp="/Users/valeriepineaunoel/Desktop/stitchedImage.tiff")
shift = stitch.calculate_shift_PCC(index1=0, index2=1, directory=sourcePath + "/IntensityCorrection")
print(f"SHIFT : {shift}")
#tiff.imwrite("/Users/valeriepineaunoel/Desktop/convolvedImage.tiff", shift)