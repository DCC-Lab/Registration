# Stitching-Scrapbooking
This Python code uses straighforward and well-known image-treatment methods --such as correlations, Fourier transforms and arithmetic calculations-- to stitch images that contain overlapped regions. Its development is inspired by the principle of scrapbooking (hence its name). In summary, the stitching algorithm:
  - creates a black image of the right size, which is used as a background image to paste individual images on it;
  - estimates the lateral shift of each image according to the position of its neighbouring images;
  - pastes each image at the right lateral coordinates on the background image.

