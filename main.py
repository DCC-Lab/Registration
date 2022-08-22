from stitching import *
from imageTreatment import *
import tifffile as tiff

sourcePath = "/Users/valeriepineaunoel/Desktop/20220620-334B-1A-Cryoprotectant-Ave30-tiling-zoom2-003/LineCorrection/IntensityCorrection"
tileDimensions = [176, 1]

stitch = Stitching(sourceDir=sourcePath, tileD=tileDimensions, imageSize=[1024,512], isIntensityCorrection=False, isMirrored=True, isFlipped=False)