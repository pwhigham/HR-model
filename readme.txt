This repository contains the basic code required to perform the Hill-Robertson experiments.

The files are:

fixedspace.cpp - the C code to produce the HR model.  Set to N=64.  Needs to be compiled for use on a 
computer, however, fixed64.exe will run on windows machines.  To run for other population sizes the 
file must be edited and N changed (and then recompiled).  This was done to try and make the software
run as fast as possible -- and it wasn't too inconvenient to have to build several different executables for different size networks (only done for the initial analysis of varying N).

The other C programs implement the supplementary material used for analyzing the steady-state behaviour of the HR model. 

********************************************************
NOTE: The output folder has been zipped to reduce size.
********************************************************

networks.R - R code to build the network matrices.  Assumes that you will place each network in a folder, and these are named 1.txt, 2.txt, etc.  fixedspace.cpp assumes this naming convention for the files.

Fig3.bat - example windows batch files to run panmictic/ring/lattice/sf1/sf2/star models (N=64) varying Nbeta and q.

Fig5.bat - example windows batch files to run panmictic/ring/lattice/sf1/sf2/star models (N=64) varying Nbeta and Nc. See Fig5.bat for comments regarding how to run this batch file.

Fig7.bat - example windows batch file to run data for Figure 7 haplotype data, varying Nibeta and Nc.

The R file:  paperFigures.R has a large number of example functions that produce the figures (and additional material) for the paper.  Note that the default paths used in these functions assume that the folder containing the data is Spacefixed/output/*.  


** Note **  The main aspect to validate is the C program fixedspace.cpp which runs the HR model.

  