# DFT-drift-correction
Script for microscope drift correction by cross-correlating the brightfield images of cells.

*Author: Aleks Ponjavic, University of Leeds.*

*Based on algorithm from Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, "Efficient subpixel image registration algorithms," Opt. Lett. 33, 156-158 (2008).*

*Function for reading TIFF files >4gb from: Tristan Ursell (2020). imread_big -- read in TIFF stacks larger than 4GB (https://www.mathworks.com/matlabcentral/fileexchange/61376-imread_big-read-in-tiff-stacks-larger-than-4gb), MATLAB Central File Exchange. Retrieved October 9, 2020.*

Disclaimer
I have made every effort to evaluate the proper working of this code under many different conditions. However, it is the responsibility of the user to ensure that this registration code is adequate and working correcntly for their application.
Feel free to e-mail me with questions or comments at A.Ponjavic@leeds.ac.uk.

To run this algorithm, the main efficient_subpixel_registration.m script, as well as the auxiliary functions dftregistration.m and imread_big.m need to be in the working directory. 

To perform drift correction, open the efficient_subpixel_registration.m file, edit the necessary parameters, and run.

## Inputs - files
- Main folder with data.
- Reference brightfield images.
  -  Since reading large files slows the scripts down a lot, there is an option to divide the single brightfield series into several files. The user needs to specify the last frame number in each file.
- Localisations file for correction

## Inputs - parameters
 - Pixel size in nm.
 - Upsampling factor: image cross-correlation will be with accuracy 1/upFactor pixels.
 - Are the images flipped?
 - Occasional correction with averaging: used to remove the need for continuous BF imaging and also to improve cross-correlation precision.
   - Imaging might be done in two modes: taking a burst of images once per cycle (i.e. WL_image_WL_image_WL_image_WL), or at the beginning and end of each cycle (i.e. WL_image_WL; WL_image_WL; WL_image_WL). NB! Even if images are taken once per cycle, there needs to be a WL image both at the very beginning and end of the series!
    - Two parameters that need to be specified in occasional mode are the number of super-res images taken in one cycle, and the number of BF images taken in a burst for averaging purposes.

## Output
- Output is written in the same folder as input, with an appended user-specified string.
- If desired, whitelight images can be plotted while the script runs.
- If input is a single bead, SD before and after correction will be calculated nd displayed for comparison.
