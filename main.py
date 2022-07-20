from functions import *

sourcePath = "/Users/valeriepineaunoel/Desktop/test"
files = listNameOfFiles(directory=sourcePath)
tileDimensions = [2, 2]
manualShift = [384, 39] # [y,x]. Positive y = going down. Positive x = going right. 

# initialize images in numpy arrays
imagenp0 = read_file(filePath=sourcePath + "/" + files[0], imageType="numpy")
imagenp1 = read_file(filePath=sourcePath + "/" + files[1], imageType="numpy")
imagenp2 = read_file(filePath=sourcePath + "/" + files[2], imageType="numpy")
imagenp3 = read_file(filePath=sourcePath + "/" + files[3], imageType="numpy")

# initialise images as pillow images
imagePIL0 = read_file(filePath=sourcePath + "/" + files[0], imageType="PIL")
imagePIL1 = read_file(filePath=sourcePath + "/" + files[1], imageType="PIL")
imagePIL2 = read_file(filePath=sourcePath + "/" + files[2], imageType="PIL")
imagePIL3 = read_file(filePath=sourcePath + "/" + files[3], imageType="PIL")

# calculate shift with cross-phase correlation does not give the right answer. 
#coordinatesOfShift = calculate_shift(image1, image2)

#stitchedImage.save(fp="stitchedImage.tif")