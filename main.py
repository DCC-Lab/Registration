from stitching import *
from imageTreatment import *
import tifffile as tiff

sourcePath = "/Users/valeriepineaunoel/Desktop/20220822-305B-1/20220822-305B-1-Zoom2-Ave30-10%-12/LineCorrection"
tileDimensions = [44, 10]

stitch = Stitching(sourceDir=sourcePath, tileD=tileDimensions, imageSize=[1024,512], isIntensityCorrection=True, isMirrored=True, isFlipped=False)