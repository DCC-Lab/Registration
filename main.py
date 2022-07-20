from functions import *

sourcePath = "/Users/valeriepineaunoel/Desktop/test"
files = listNameOfFiles(directory=sourcePath)
tileDimensions = [2, 2]

pathFile1 = sourcePath + "/" + files[0]
pathFile2 = sourcePath + "/" + files[1]

stitchedImage = merge_images(file1=pathFile1, file2=pathFile2)
stitchedImage.save(fp="stitchedImage.tif")