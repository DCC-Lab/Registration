from functions import *

sourcePath = "/Users/valeriepineaunoel/Desktop/test"
files = listNameOfFiles(directory=sourcePath)
tileDimensions = [2, 2]

pathFile1 = sourcePath + "/" + files[1]
pathFile2 = sourcePath + "/" + files[3]
print(files[1])
print(files[3])
image1 = read_file(pathFile1)
image2 = read_file(pathFile2)

allo = calculate_shift(image1, image2)

#stitchedImage.save(fp="stitchedImage.tif")