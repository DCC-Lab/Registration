import os
import fnmatch

def listNameOfFiles(directory: str, extension="tif") -> list:
	"""
	Fetch files name. 
	Does not consider .DS_Store files. 
	Returns a list of the names of the files. 
	"""
	foundFiles = []
	for file in os.listdir(directory):
		if fnmatch.fnmatch(file, f'*.{extension}'):
			if file == ".DS_Store":
				pass
			else:
				foundFiles.append(file)

	foundFiles.sort()
	return foundFiles

def read_file(file_path):
	"""
	Reads the .tif file and convert them in a np.array. 
	Returns the file as a np.array. 
	"""
	image_array = tiff.imread(file_path)
	return image_array