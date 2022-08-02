from stitching import *

sourcePath = "/Users/valeriepineaunoel/Desktop/test20%/IntensityCorrection"
tileDimensions = [2, 3]

stitch = Stitching(sourceDir=sourcePath, tileD=tileDimensions, imageSize=[1024,512], isIntensityCorrection=False, isMirrored=False)
stitchedImage = stitch.stitching_scrapbooking_allImages()
stitchedImage.save(fp="/Users/valeriepineaunoel/Desktop/stitchedImage.tiff")