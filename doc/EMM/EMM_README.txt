Enhanced Magnetic Model Software and support documents
=======================================================


Sublibrary Files
================
GeomagnetismLibrary.c	                Geomagnetism library, C functions
GeomagnetismHeader.h			Geomagnetism library, C header file 


Main Programs
===============
emm_sph_point.c			Command prompt version for single point computation
emm_sph_grid.c			Grid, profile and time series computation, C main function
emm_sph_file.c			C program which takes a coordinate file as input


Data Files
===============
EMM20XX.COF				EMM Static Coefficients file (2000 - 2017)
EMM20XXSV.COF				EMM secular variation coefficient file (2000 - 2017)


Supporting Documents
=====================
Geomagnetism_Library_Software_Manual.pdf		Description of subroutines used in EMM software
EMM2017TestValues.pdf   	Test values for  EMM spherical harmonic calculation.


Excecutables (Windows and Linux)
============
emm_point.exe					Command Prompt program for single point
emm_grid.exe					Grid, time series or profile
emm_file.exe					File processing

Test Files
===============
EMM2017_TEST_VALUES.pdf				Test values for EMM2017
sample_coords.txt				Sample input file for program emm_sph_file.exe 
sample_output_file.txt				Sample output file for program emm_sph_file.exe


Compiling with gcc
===================
gcc inputfile [dependencies] -lm -o outputfile
For example, the emm_sph_file.c can be compiled as
gcc emm_sph_file.c GeomagnetismLibrary.c -lm  -o emm_sph_file.exe



Model Software Support
======================

*  National Centers for Environmental Information
*  NOAA E/NE42
*  325 Broadway
*  Boulder, CO 80303 USA
*  Attn: Manoj Nair or Adam Woods
*  Phone:  (303) 497-4642 or -4460
*  Email:  Geomag.Models@noaa.gov
For more details about the World Magnetic Model visit 
http://www.ngdc.noaa.gov/geomag/EMM/
