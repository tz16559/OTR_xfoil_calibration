
The analysis chain starts with the png images from the Linux OTR DAQ, these are PNGs and are the wrong size, we first convert to PGM and crop to the correct size

>python convert_png_to_pgm_and_crop.py path_to_input_files path_to_output

This will traverse all png files in path_to_input_files and recreate this structure in path_to_output (which does not need to exist before running this)

you may now wish to stack multiple images together to improve the signal to noise ratio and subtract the pedestal. (left as an excercise for the reader)

Next, these images need to be converted to a root histogram, edit the convert_to_root() function in convert_to_root.C to point to the location of the pgm files, then:

>root -l
>.L convert_to_root.C
>convert_to_root()


This will create a root file containing the images in TH2 histograms

The spot finding algorithm can now be used, edit ./calibration_analysis/src/holeFinder.cpp to look for this root file 

>cd ./calibration_analysis/
>make
>./bin/holeFinder files_to_analyse

where files_to_analyse is a list of images that you want to analyse e.g. V_00000001 V_00000002

this will output a series of plots showing spot positions, differences between the two images etc.

