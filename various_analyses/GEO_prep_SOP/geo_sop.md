## This SOP documents procedures to prepare the txt data files for microarray data submission to GEO

###1. cd to the "Standard Analysis" folder

###2. copy the microarray data files to a working directory:

`find . -name "*Data.xls" -exec cp {} GEO_working_dir/ \;` 

###3. generate fileList.txt to a keep a list for all data files:

`ls *xls > fileList.txt`

###4. run the R script to convert the data files into txt files

`Rscript makeGEOrawTxt.R`


