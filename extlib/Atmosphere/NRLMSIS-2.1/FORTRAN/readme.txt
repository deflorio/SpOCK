#######################################################################
 MSIS® (NRL-SOF-014-1) SOFTWARE
 NRLMSIS® empirical atmospheric model software. Use is governed by the
 Open Source Academic Research License Agreement contained in the file
 nrlmsis2.1_license.txt, which is part of this software package. BY
 USING OR MODIFYING THIS SOFTWARE, YOU ARE AGREEING TO THE TERMS AND
 CONDITIONS OF THE LICENSE.  
#######################################################################

NRLMSIS 2.1 Whole-Atmosphere Empirical Model of Temperature and Neutral Species
  Densities

VERSION HISTORY
  08 MAR 19 Version 1.97 (Beta version)
  26 MAY 20 Version 2.0 (Release version)
  04 APR 22 Version 2.1 (Release version with NO density)

AUTHORS
  John Emmert (john.emmert@nrl.navy.mil)
  Douglas Drob (douglas.drob@nrl.navy.mil)
  McArthur Jones Jr. (mcarthur.jones@nrl.navy.mil)
  
REFERENCE FOR NRLMSIS 2.0
  Emmert, J. T., Drob, D. P., Picone, J. M., Siskind, D. E., Jones, M. Jr., 
  Mlynczak, M. G., et al. (2021). NRLMSIS 2.0: A whole-atmosphere empirical model
  of temperature and neutral species densities. Earth and Space Science, 8, 
  e2020EA001321. https://doi.org/10.1029/2020EA001321

PACKAGE CONTENTS
  readme.txt                This file
  nrlmsis2.1_license.txt    Open Source Academic Research License Agreement
  msis2.1_test.F90          Test program
  msis_init.F90             Subroutines to initialize the model, set switches and
                              options, and load parameter file
  msis_gtd8d.F90            Subroutine to evaluate the model using the legacy
                              interface
  msis_calc.F90             Subroutine for evaluating the model using the new
                              interface
  msis_constants.F90        Module containing model constants
  msis_gfn.F90              Subroutines to calculate horizontal expansion
                              functions
  msis_tfn.F90              Subroutines to calculate the vertical temperature
                              profile
  msis_dfn.F90              Subroutines to calculate vertical density profiles
  msis_utils.F90            Subroutines to convert between geodetic height and
                              geopotential height, and other support subroutines
  msis21.parm               Binary data file containing model parameters
  msis2.1_test_in.txt       ASCII file containing input for test program.
  msis2.1_test_ref_dp.txt   ASCII file containing expected output of test program
                             (double-precision internally)
 
RELEASE NOTES: MODEL FORMULATION
  Minor changes to the NRLMSIS 2.0 formulation include:
  - Addition of new terms to support fitting of NO densities.
  - Reorganization of support subroutines (alt2gph, gph2alt, bspline, dilog)
    into msis_utils module.
 
RELEASE NOTES: PARAMETER ESTIMATION
  - NO density parameters were tuned to six NO data sets (ENVISAT/MIPAS, SNOE,
    ACE/FTS, AIM/SOFIE, UARS/HALOE, and Odin/SMR).
  - Temperature and all other species densities are the same as in NRLMSIS 2.0.

COMPILING THE MODEL CODE
  The model package was tested on Windows, Linux, and Mac systems using the
  following Fortran compilers and compile statements:
    gfortran 4.8.5, 7.5.0, 9.3.0
      gfortran -O3 -cpp -o msis2.1_test.exe msis_constants.F90 msis_utils.F90
      msis_init.F90 msis_gfn.F90 msis_tfn.F90 msis_dfn.F90 msis_calc.F90
      msis_gtd8d.F90 msis2.1_test.F90
        NOTES:
        - The following optimization flags may improve performance:
            -march=native -ffast-math
    Intel 2017.2.163, 18.0.1.156, 2021.1
      ifort -O2 -fpp -o msis2.1_test.exe msis_constants.F90 msis_utils.F90
      msis_init.F90 msis_gfn.F90 msis_tfn.F90 msis_dfn.F90 msis_calc.F90
      msis_gtd8d.F90 msis2.1_test.F90
        NOTES:
        - The following optimization flags may improve performance:
            Windows:     -Qipo -QxHost
            Linux/macOS: -ipo -xHost
  For double precision, add the flag -DDBLE. Double precision is not necessary
  for most applications, but for testing purposes it ensures that the test
  output exactly matches the expected output in msis2.1_test_ref_dp.txt,
  regardless of the compiler or compiler settings.

INITIALIZING AND RUNNING THE MODEL
  - The model must be initialized using the MSISINIT subroutine, which sets
    switches and options and loads the model parameter values from a file.
  - The switch_legacy optional argument to MSISINIT performs the same function
    as TSELEC(SW) in NRLSMSISE-00, except that switches 15-25 are not used in
    NRLMSIS 2.1. The change in the switch-setting call is illustrated as
    follows, where SW is the 25-element array of switches:
        NRLMSISE-00: CALL TSELEC(SW)
        NRLMSIS 2.1: call msisinit(switch_legacy=SW)
  - The MSISCALC subroutine checks for initialization and does a default
    initialization if necessary. This self-initialization will be removed in
    future versions.
  - The model can be called using either the legacy interface (subroutine
    GTD8D) or the new interface (subroutine MSISCALC).
  - Details of the input and output arguments of MSISINIT, GTD8D, and MSISCALC
    are provided in the headers of the respective source code files.
  
ACKNOWLEDGEMENTS
  This work was supported by the Office of Naval Research and NASA.