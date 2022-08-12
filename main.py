from stitching import *
from imageTreatment import *
import tifffile as tiff

sourcePath = "/Users/valeriepineaunoel/Desktop/20220621-6umPolystyreneBeads-overlap20%-zoom2/LineCorrection"
tileDimensions = [6, 4]

stitch = Stitching(sourceDir=sourcePath, tileD=tileDimensions, imageSize=[1024,512], isIntensityCorrection=True, isMirrored=True, isFlipped=False)