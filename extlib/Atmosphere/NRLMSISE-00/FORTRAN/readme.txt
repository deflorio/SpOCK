NRLMSISE-00 Model	2001

Authors:  M. Picone, A.E. Hedin, D. Drob
Naval Research Laboratory<br>

Parameter: Neutral densities and temperature from ground to
thermosphere

Brief Description: 
The NRLMSIS-00 empirical atmosphere model was developed by Mike Picone, Alan Hedin, 
and Doug Drob based on the MSISE90 model. The main differences to MSISE90 are noted
in the comments at the top of the computer code. They involve (1) the extensive use
of drag and accelerometer data on total mass density,  (2) the addition of a
component to the total mass density that accounts for possibly significant 
contributions of O+ and hot oxygen at altitudes above 500 km, and (3) the inclusion
of the SMM UV occultation data on [O2]. 
The MSISE90 model describes the neutral temperature and densities in Earth's 
atmosphere from ground to thermospheric heights. Below 72.5 km the model is
primarily based on the MAP Handbook (Labitzke et al., 1985) tabulation of
zonal average temperature and pressure by Barnett and Corney, which was also
used for the CIRA-86. Below 20 km these data were supplemented with averages
from the National Meteorological Center (NMC). In addition, pitot tube,
falling sphere, and grenade sounder rocket measurements from 1947 to 1972 were 
taken into consideration. Above 72.5 km MSISE-90 is essentially a revised
MSIS-86 model taking into account data derived from space shuttle flights and
newer incoherent scatter results. For someone interested only in the
thermosphere (above 120 km), the author recommends the MSIS-86 model. MSISE
is also not the model of preference for specialized tropospheric work. It is
rather for studies that reach across several atmospheric boundaries.

Availability: 
(1) FORTRAN source code is available in this directory:
nrlmsise00_driver.for  Example driver program  
                       (f77 -o nrlmsis nrlmsise00_driver.for nrlmsise00_sub.for)
nrlmsise00_output.for  Output from driver program 
nrlmsise00_sub.for     NRLMSISE90 FORTRAN subroutines

(2) FORTRAN source code is also available from http://uap-www.nrl.navy.mil/models_web/msis/msis_home.htm. 

(3) A C version of the code was provided by D. Brodowski and is available in directory 
/nrlmsis00_c_version, or can be downloaded from D. Brodowski's site at 
http://www.brodo.de/english/pub/nrlmsise/

National Space Science Data Center / Dieter Bilitza / bilitza@gsfc.ansa.gov /Nov 5,01

CORRECTIONS:
July 14, 2003 - Beta version was replaced with offical version. Very minor changes 
	primarily concerning the O2 and He profiles in the lower thermosphere and upper 
	mesosphere.

The following modificatiions were made by Lutz Rastaetter at the CCMC:

March 23, 2011 - declare data types: CHARACTER*4 NAME,ISTIME,ISDATE 
                                 and INTEGER IMR 

Sept. 16, 2014 - array dimensions were increased: P(1) -> P(150), AP(1) -> AP(7)

July 9, 2019   - closed model coverage gap at altitude ZN3(1) (32.5 km) 
