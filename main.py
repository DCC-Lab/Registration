from stitching import *
from imageTreatment import *
import tifffile as tiff

sourcePath = "/Users/valeriepineaunoel/Desktop/20220621-20%-334A-1/LineCorrection"
tileDimensions = [11, 3]

stitch = Stitching(sourceDir=sourcePath, tileD=tileDimensions, imageSize=[1024,512], isIntensityCorrection=True, shiftEstimation="PCC", isMirrored=True, isFlipped=False)
#stitchedImage = stitch.stitching_with_known_shifts()
#plt.imshow(stitchedImage)
#plt.show()
#a = input("Want to save? (Answer y/n)")
#if a == "y":
#	path = "/Users/valeriepineaunoel/Desktop/stitchedImage.tiff"
#	stitchedImage.save(fp=path)
#	print(r"Saved as /Users/valeriepineaunoel/Desktop/stitchedImage.tiff")