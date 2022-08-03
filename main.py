from stitching import *
from imageTreatment import *
import tifffile as tiff

sourcePath = "/Users/valeriepineaunoel/Desktop/20220621-6umPolystyreneBeads-overlap10%-zoom2/LineCorrection"
tileDimensions = [6, 4]

#img = ImageTreatment(sourceDir=sourcePath)
#
#average = img.create_average_image()
#tiff.imwrite("/Users/valeriepineaunoel/Desktop/average.tiff", average)
#
#inv = img.inverse_pixels(image=average)
#tiff.imwrite("/Users/valeriepineaunoel/Desktop/inv.tiff", inv)
#res = img.rescale_image(image=inv)
#tiff.imwrite("/Users/valeriepineaunoel/Desktop/rescale.tiff", res)
#corr = img.create_intensity_correction_image(image=average)
#tiff.imwrite("/Users/valeriepineaunoel/Desktop/corr.tiff", corr)
#
#enveloppe = img.correct_intensity_envelop()

stitch = Stitching(sourceDir=sourcePath, tileD=tileDimensions, imageSize=[1024,512], isIntensityCorrection=True, isMirrored=True, isFlipped=False)
stitchedImage = stitch.stitching_scrapbooking_allImages()
plt.imshow(stitchedImage)
plt.show()
a = input("Want to save? (Answer y/n)")
if a == "y":
	path = "/Users/valeriepineaunoel/Desktop/stitchedImage.tiff"
	stitchedImage.save(fp=path)
	print(r"Saved as /Users/valeriepineaunoel/Desktop/stitchedImage.tiff")