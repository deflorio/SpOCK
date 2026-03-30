//==========================================================================
/*
*    Copyright 2020 Sergio De Florio
*    All rigths reserved
*
*    This file is part of SpOCK
* 
*    SpOCK is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation version 3
* 
*    SpOCK is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
*    GNU General Public License for more details.
* 
*    You should have received a copy of the GNU General Public License
*    along with SpOCK. If not, see <https://www.gnu.org/licenses/>.
*
*/
//==========================================================================

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <math.h>
//#include <windows.h>

#include <typeinfo>
#include <thread>
#include <string>

#include <Constants.h>
#include <VarTypes.h>
#include <Transformations.h>
#include <IO_utils.h>
#include <Interpolation.h>

// External libraries: Eigen
#include <Eigen/Core>
#include <boost/math/interpolators/cubic_hermite.hpp>
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>

using namespace std;
using namespace constants;

using namespace boost;
using namespace boost::math::interpolators;

int main(int argc, char *argv[])
    {
    if(argc != 5)
      {
      // Tell the user how to run the program
      cerr << "Usage: " << argv[0] << " Path to files       OrbitEphemeris.csv       Type of interpolation (LN = linear, CH = cubic Hermite, LG = Lagrange, FILL = connect ephemerides point with constant values)       Number of interpolation points (even number)       Desired interpolation step [s]" << endl;
      return 1;
      }
    
    ////////////////////////////// Variables //////////////////////////////
    string OrbEph_file(argv[1]);
    string IntplType(argv[2]);
    int intpl_size = atoi(argv[3]);
    if ( !(intpl_size % 2 == 0) )
        {
        cerr << "\nNumber of interpolation points has to be even" << endl;
        return 1;
        }
    if( IntplType.compare("LN") == 0 || IntplType.compare("FILL") == 0 ) intpl_size = 2;
    string IntplStep_str = argv[4];
    double IntplStep = atof(argv[4]);
    
    if( IntplType.compare("LN") != 0 && IntplType.compare("CH") != 0 && IntplType.compare("LG") != 0 && IntplType.compare("FILL") != 0 )
        {
        cerr << "\nPlease insert correct interpolation type (LN, CH, LG or FILL)" << endl;
        return 1;
        }
        
    if(IntplStep < 0.1)
        {
        cerr << "\nInterpolation step value has to be >= 0.1" << endl;
        return 1;
        }
    
    string InplOrbEph_csvfilename;
    InplOrbEph_csvfilename = OrbEph_file.substr(0, OrbEph_file.size() - 4) + "_" + IntplType + "_" + IntplStep_str + "s" + "_step.csv";
    
    ofstream InplOrbEph_file;
    InplOrbEph_file.open(InplOrbEph_csvfilename);
    
    
    
    MatrixXd loaded_OrbEph, timeposvelECI, timeposvelECEF;
    Vector6d ECIstate, ECEFstate;
    //Vector6d RTNstate1, RTNstate2, RTN_diff;
    Mat3x3d T_RTN2ECI;
    int OrbEph_rows;
    double timeGPS, intpl_initime, intpl_endtime;
    VectorXd intpl_initstate, time_interpol, t_IntplInterval, xECI_IntplIn, yECI_IntplIn, zECI_IntplIn, vxECI_IntplIn, vyECI_IntplIn, vzECI_IntplIn, xECEF_IntplIn, yECEF_IntplIn, zECEF_IntplIn, vxECEF_IntplIn, vyECEF_IntplIn, vzECEF_IntplIn;
    VectorNi<2> IntplInt;
    
    bool header = false;
    bool valid = true;
    const int Eph_cols = 14;
    
    timeposvelECI = MatrixXd::Zero(intpl_size, 7); // GPSsecs, x, y, z, vx, vy, vz
    timeposvelECEF = MatrixXd::Zero(intpl_size, 7); // GPSsecs, x, y, z, vx, vy, vz
    
    ////////////////////////////// Load orbit ephemeris 1 //////////////////////////////
    loaded_OrbEph = read_csvfile(OrbEph_file.c_str(), Eph_cols, header);
    OrbEph_rows = loaded_OrbEph.rows();
    
    // Compute sampling time of OrbEph
    VectorXd OrbEph_time_1, OrbEph_time_2, timediffs;
    
    OrbEph_time_1 = VectorXd::Zero(OrbEph_rows - 1);
    OrbEph_time_2 = VectorXd::Zero(OrbEph_rows - 1);
    timediffs = VectorXd::Zero(OrbEph_rows - 1);
    
    OrbEph_time_1 = loaded_OrbEph.col(1).segment(0,OrbEph_rows-1);
    OrbEph_time_2 = loaded_OrbEph.col(1).segment(1,OrbEph_rows-1);
    
    timediffs = OrbEph_time_2 - OrbEph_time_1;
    
    double OrbEph_step = timediffs.mean();
    
    // Exit if there is no common epochs between the time intervals of the two ephemerides
    if(OrbEph_step == IntplStep)
        {
        cerr << "\nNo need of interpolation, chosen interpolation step is equal to input ephemerides step" << endl;
        return 1;
        }
    
    ////////////////////////////// Create splines for interpolation of ECI orbit state of ephemeris 1 //////////////////////////////
    int intpl_idx = intpl_size/2 - 1;
    int end_intpl_idx = OrbEph_rows - intpl_size/2 - 1;
    
    intpl_initime = loaded_OrbEph(intpl_idx,1);
    intpl_endtime = loaded_OrbEph(end_intpl_idx,1);
    intpl_initstate = loaded_OrbEph.row(intpl_idx).segment(2,12);
    timeGPS = intpl_initime;
    
    time_interpol = loaded_OrbEph.col(1);
    
    ////////////////////////////// Interpolation of ECI orbit state in ephemeris 1 on times of ephemeris 2 and computation of orbital elements //////////////////////////////
    double secs = 0.0;
    
    InplOrbEph_file << fixed << secs << "," << timeGPS << "," << loaded_OrbEph(intpl_idx,2) << "," << loaded_OrbEph(intpl_idx,3) << "," << loaded_OrbEph(intpl_idx,4) << "," << loaded_OrbEph(intpl_idx,5) << "," << loaded_OrbEph(intpl_idx,6) << "," << loaded_OrbEph(intpl_idx,7) << "," << loaded_OrbEph(intpl_idx,8) << "," << loaded_OrbEph(intpl_idx,9) << "," << loaded_OrbEph(intpl_idx,10) << "," << loaded_OrbEph(intpl_idx,11) << "," << loaded_OrbEph(intpl_idx,12) << "," << loaded_OrbEph(intpl_idx,13) << endl;
    
    timeGPS += IntplStep;
    secs += IntplStep;
    
    while(timeGPS <= intpl_endtime)
        {
        IntplInt = IntplInterval(time_interpol, timeGPS, intpl_size, valid);
        intpl_idx = IntplInt(0);
        
        if( IntplInt(0) == IntplInt(1) ) // No need of interpolation
            {
            InplOrbEph_file << fixed << secs << "," << timeGPS << "," << loaded_OrbEph(intpl_idx,2) << "," << loaded_OrbEph(intpl_idx,3) << "," << loaded_OrbEph(intpl_idx,4) << "," << loaded_OrbEph(intpl_idx,5) << "," << loaded_OrbEph(intpl_idx,6) << "," << loaded_OrbEph(intpl_idx,7) << "," << loaded_OrbEph(intpl_idx,8) << "," << loaded_OrbEph(intpl_idx,9) << "," << loaded_OrbEph(intpl_idx,10) << "," << loaded_OrbEph(intpl_idx,11) << "," << loaded_OrbEph(intpl_idx,12) << "," << loaded_OrbEph(intpl_idx,13) << endl;
            
            timeGPS += IntplStep;
            secs += IntplStep;
            }
        else
            {
            t_IntplInterval = time_interpol.segment(intpl_idx,intpl_size);
            //cout << fixed << t_IntplInterval.transpose() << endl;
                
            // ECI
            timeposvelECI.col(0) = t_IntplInterval;
            timeposvelECI.block(0,1,intpl_size,6) = loaded_OrbEph.block(intpl_idx,2,intpl_size,6);
            
            xECI_IntplIn = timeposvelECI.col(1);
            yECI_IntplIn = timeposvelECI.col(2);
            zECI_IntplIn = timeposvelECI.col(3);
            vxECI_IntplIn = timeposvelECI.col(4);
            vyECI_IntplIn = timeposvelECI.col(5);
            vzECI_IntplIn = timeposvelECI.col(6);
            
            // ECEF
            timeposvelECEF.col(0) = t_IntplInterval;
            timeposvelECEF.block(0,1,intpl_size,6) = loaded_OrbEph.block(intpl_idx,8,intpl_size,6);
            
            xECEF_IntplIn = timeposvelECEF.col(1);
            yECEF_IntplIn = timeposvelECEF.col(2);
            zECEF_IntplIn = timeposvelECEF.col(3);
            vxECEF_IntplIn = timeposvelECEF.col(4);
            vyECEF_IntplIn = timeposvelECEF.col(5);
            vzECEF_IntplIn = timeposvelECEF.col(6);
                    
            // Interpolation
            while( timeGPS < t_IntplInterval(intpl_size/2) )
                {
                //cout << fixed << timeGPS << endl;
                
                ECIstate = Vector6d::Zero();
                ECEFstate = Vector6d::Zero();
                
                if( IntplType.compare("CH") == 0 )
                    {
                    CH_Interpolation(timeGPS, timeposvelECI, ECIstate);
                    CH_Interpolation(timeGPS, timeposvelECEF, ECEFstate);
                    }
                if( IntplType.compare("LN") == 0 )
                    {
                    ECIstate(0) = LERP(t_IntplInterval, xECI_IntplIn, timeGPS, valid);
                    ECIstate(1) = LERP(t_IntplInterval, yECI_IntplIn, timeGPS, valid);
                    ECIstate(2) = LERP(t_IntplInterval, zECI_IntplIn, timeGPS, valid);
                    ECIstate(3) = LERP(t_IntplInterval, vxECI_IntplIn, timeGPS, valid);
                    ECIstate(4) = LERP(t_IntplInterval, vyECI_IntplIn, timeGPS, valid);
                    ECIstate(5) = LERP(t_IntplInterval, vzECI_IntplIn, timeGPS, valid);
                    
                    ECEFstate(0) = LERP(t_IntplInterval, xECEF_IntplIn, timeGPS, valid);
                    ECEFstate(1) = LERP(t_IntplInterval, yECEF_IntplIn, timeGPS, valid);
                    ECEFstate(2) = LERP(t_IntplInterval, zECEF_IntplIn, timeGPS, valid);
                    ECEFstate(3) = LERP(t_IntplInterval, vxECEF_IntplIn, timeGPS, valid);
                    ECEFstate(4) = LERP(t_IntplInterval, vyECEF_IntplIn, timeGPS, valid);
                    ECEFstate(5) = LERP(t_IntplInterval, vzECEF_IntplIn, timeGPS, valid);
                    }
                if( IntplType.compare("LG") == 0 )
                   {
                    ECIstate(0) = LG_Interpolation(t_IntplInterval, xECI_IntplIn, timeGPS, valid);
                    ECIstate(1) = LG_Interpolation(t_IntplInterval, yECI_IntplIn, timeGPS, valid);
                    ECIstate(2) = LG_Interpolation(t_IntplInterval, zECI_IntplIn, timeGPS, valid);
                    ECIstate(3) = LG_Interpolation(t_IntplInterval, vxECI_IntplIn, timeGPS, valid);
                    ECIstate(4) = LG_Interpolation(t_IntplInterval, vyECI_IntplIn, timeGPS, valid);
                    ECIstate(5) = LG_Interpolation(t_IntplInterval, vzECI_IntplIn, timeGPS, valid);
                    
                    ECEFstate(0) = LG_Interpolation(t_IntplInterval, xECEF_IntplIn, timeGPS, valid);
                    ECEFstate(1) = LG_Interpolation(t_IntplInterval, yECEF_IntplIn, timeGPS, valid);
                    ECEFstate(2) = LG_Interpolation(t_IntplInterval, zECEF_IntplIn, timeGPS, valid);
                    ECEFstate(3) = LG_Interpolation(t_IntplInterval, vxECEF_IntplIn, timeGPS, valid);
                    ECEFstate(4) = LG_Interpolation(t_IntplInterval, vyECEF_IntplIn, timeGPS, valid);
                    ECEFstate(5) = LG_Interpolation(t_IntplInterval, vzECEF_IntplIn, timeGPS, valid);
                   }
                if( IntplType.compare("FILL") == 0 )
                    {
                    ECIstate = loaded_OrbEph.row(intpl_idx).segment(2,6);
                    ECEFstate = loaded_OrbEph.row(intpl_idx).segment(8,6);    
                    }
                
                InplOrbEph_file << fixed << secs << "," << timeGPS << "," << ECIstate(0) << "," << ECIstate(1) << "," << ECIstate(2) << "," << ECIstate(3) << "," << ECIstate(4) << "," << ECIstate(5) << "," << ECEFstate(0) << "," << ECEFstate(1) << "," << ECEFstate(2) << "," << ECEFstate(3) << "," << ECEFstate(4) << "," << ECEFstate(5) << endl;
                
                timeGPS = round( (timeGPS + IntplStep)*1e2 )/1e2;
                secs += IntplStep;
                }
            }
        }
        
    // Close files        
    InplOrbEph_file.close();
    
    return(0);
    } // End of main()
