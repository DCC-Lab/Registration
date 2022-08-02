from stitching import *

sourcePath = "/Users/valeriepineaunoel/Desktop/test20%/IntensityCorrection"
tileDimensions = [2, 3]
# [x,y]. Positive y = going down. Positive x = going right. When comparing image2 with image1 : 
verticalShift = [47, 363]
horizontalShift = [163, -12]

stitch = Stitching(sourceDir=sourcePath, tileD=tileDimensions, imageSize=[1024,512], isIntensityCorrection=False, isMirrored=True)
shift = stitch.estimate_shift(index=0, stitchingSide="V")
print(f"PCC : {shift}")

#stitchedImage = stitch.stitching_scrapbooking_allImages()
#stitchedImage.save(fp="/Users/valeriepineaunoel/Desktop/stitchedImage.tiff")