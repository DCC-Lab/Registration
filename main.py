from stitching import *
from imageTreatment import *
import tifffile as tiff

sourcePath = "/Users/valeriepineaunoel/Library/Mobile Documents/com~apple~CloudDocs/Documents/Stitching-Scrapbooking/tests/testDataset"
tileDimensions = [2, 3]

stitch = Stitching(sourceDir=sourcePath, tileD=tileDimensions, imageSize=[1024,512], isIntensityCorrection=False, shiftEstimation="PCC", isMirrored=True, isFlipped=False)
stitchedImage = stitch.stitching_scrapbooking_allImages()
plt.imshow(stitchedImage)
plt.show()
a = input("Want to save? (Answer y/n)")
if a == "y":
	path = "/Users/valeriepineaunoel/Desktop/stitchedImage.tiff"
	stitchedImage.save(fp=path)
	print(r"Saved as /Users/valeriepineaunoel/Desktop/stitchedImage.tiff")