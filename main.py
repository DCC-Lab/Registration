from stitching import *
from imageTreatment import *
import tifffile as tiff

sourcePath = "/Users/valeriepineaunoel/Desktop/20220621-334A-1-Cryoprotectant-Ave30-tiling-zoom2/LineCorrection"
tileDimensions = [229, 81]

stitch = Stitching(sourceDir=sourcePath, tileD=tileDimensions, imageSize=[1024,512], isIntensityCorrection=True, isMirrored=True, isFlipped=False)