from functions import *

sourcePath = "/Users/valeriepineaunoel/Desktop/test"
files = listNameOfFiles(directory=sourcePath)
tileDimensions = [2, 2]
# [x,y]. Positive y = going down. Positive x = going right. When comparing image2 with image1 : 
verticalShift = [-39, 384]
horizontalShift = [-169, -12]

# initialize images in numpy arrays
imagenp0 = read_file(filePath=sourcePath + "/" + files[0], imageType="numpy")
imagenp1 = read_file(filePath=sourcePath + "/" + files[1], imageType="numpy")
#image_offset = read_file(filePath=sourcePath + "/" + files[2], imageType="numpy")
imagenp2 = read_file(filePath=sourcePath + "/" + files[2], imageType="numpy")
imagenp3 = read_file(filePath=sourcePath + "/" + files[3], imageType="numpy")

# initialise images as pillow images
imagePIL0 = read_file(filePath=sourcePath + "/" + files[0], imageType="PIL")
imagePIL1 = read_file(filePath=sourcePath + "/" + files[1], imageType="PIL")
imagePIL2 = read_file(filePath=sourcePath + "/" + files[2], imageType="PIL")
imagePIL3 = read_file(filePath=sourcePath + "/" + files[3], imageType="PIL")

backgroundImage = create_tile_image(tile=tileDimensions, shift=manualShift)
# width, height = backgroundImage.size
backgroundImage.save(fp="/Users/valeriepineaunoel/Desktop/backgroundImage.tif")