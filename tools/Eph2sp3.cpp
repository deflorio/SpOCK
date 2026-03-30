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

// External libraries: Eigen
#include <Eigen/Core>

extern "C"
        {
        #include "../extlib/cspice/include/SpiceUsr.h"
        }

using namespace std;
using namespace constants;

int main(int argc, char *argv[])
    {
    if(argc < 3)
      {
      // Tell the user how to run the program
      cerr << "Usage: " << argv[0] << "       Path_to_file_Orbit_Ephemeris.csv       Position/velocity flag (P for position only and V for position and velocity in sp3 file)\nOptional arguments: Satellite name (3 chars)       Data used (5 chars)       Coordinate system (5 chars)       Orbit type (3 chars)       Agency (4 chars)" << endl;
      return 1;
      }
      
    // Load SPICE Kernels
    string eop = "../data/cspice/earth_latest_high_prec.bpc";
    string leapsecond = "../data/cspice/pck00010.tpc";
    string pck_data = "../data/cspice/naif0012.tls";
    
    furnsh_c(eop.c_str( ));
    furnsh_c(leapsecond.c_str( ));
    furnsh_c(pck_data.c_str( ));
    
    ////////////////////////////// Variables //////////////////////////////
    string OrbEph_file(argv[1]);
    string posvel(argv[2]);
    
    string SatName = "L01";
    if(argc >= 4)
        {
        SatName = argv[3];
        if(SatName.length() != 3)
          {
          cerr << "Satellite name should be 3 chars long" << endl;
          return 1;
          }
        }
    
    string DataType = "ORBIT";
    if(argc >= 5)
        {
        DataType = argv[4];
        if(DataType.length() != 5)
          {
          cerr << "Data type should be 5 chars long" << endl;
          return 1;
          }
        }
        
    string CoordSys = " ITRF";
    if(argc >= 6)
        {
        CoordSys = argv[5];
        if(CoordSys.length() != 5)
          {
          cerr << "Orbit type should be 5 chars long" << endl;
          return 1;
          }
        }
    
    string OrbType = "EXT";
    if(argc >= 7)
        {
        OrbType = argv[6];
        if(OrbType.length() != 3)
          {
          cerr << "Orbit type should be 3 chars long" << endl;
          return 1;
          }
        }
        
    string AgencyName = " SDF";
    if(argc >= 8)
        {
        AgencyName = argv[7];
        if(AgencyName.length() != 4)
          {
          cerr << "Agency name should be 4 chars long" << endl;
          return 1;
          }
        }
        
    //string IntplType(argv[2]);
    //int intpl_size = atoi(argv[3]);
    
    string sp3_cvsfilename = OrbEph_file.substr(0, OrbEph_file.size() - 4) + ".sp3";
    
    ofstream sp3_file;
    sp3_file.open(sp3_cvsfilename);
    
    MatrixXd loaded_OrbEph;
    VectorXd OrbEph_epochs;
    int year, month, day, hour, min, week, MJD0;
    double GPS_UTC, sec, weeksecs;
    Vec3d posECEF, velECEF;
    int OrbEph_rows;
    
    bool header = false;
    const int Eph_cols = 14;
    
    ////////////////////////////// Load orbit ephemeris //////////////////////////////
    loaded_OrbEph = read_csvfile(OrbEph_file.c_str(), Eph_cols, header);
    OrbEph_rows = loaded_OrbEph.rows();
    
    OrbEph_epochs = loaded_OrbEph.col(1);
    
    //GPStime Epoch;
    //Epoch = GPStime( (unsigned int)OrbEph_epochs(0), 0 );
    
    Vector6d Date;
    VectorNd<2> GPSws;
    double fracOfDay, days;
    
    // Compute sampling time of OrbEph
    VectorXd OrbEph_time_1, OrbEph_time_2, timediffs;
    
    OrbEph_time_1 = VectorXd::Zero(OrbEph_rows - 1);
    OrbEph_time_2 = VectorXd::Zero(OrbEph_rows - 1);
    timediffs = VectorXd::Zero(OrbEph_rows - 1);
    
    OrbEph_time_1 = loaded_OrbEph.col(1).segment(0,OrbEph_rows-1);
    OrbEph_time_2 = loaded_OrbEph.col(1).segment(1,OrbEph_rows-1);
    
    timediffs = OrbEph_time_2 - OrbEph_time_1;
    
    double OrbEph_step = timediffs.mean();
    
    ////////////////////////////// Write sp3 file header //////////////////////////////
    //Epoch.Date(year,month,day);
    GPS_UTC = leapsec(OrbEph_epochs(0)) - timescales::TAI_GPS;
    Date = GPS2UTCdate(OrbEph_epochs(0) + GPS_UTC, "ISOC"); // (GPS - UTC) is added to have the date in GPS time and not UTC
    
    year = (int)Date(0);
    month = (int)Date(1);
    day = (int)Date(2);
    hour = (int)Date(3);
    min = (int)Date(4);
    sec = round(Date(5));
    
    GPSws = GPS2GPSws(OrbEph_epochs(0));
    week = (int)GPSws(0);
    weeksecs = GPSws(1);
    
    MJD0 = int( GPSws2MJD(GPSws(0), GPSws(1)) );
    
    days = OrbEph_epochs(0)/86400.0;
    fracOfDay = days - floor(days);
  
    // First line 
    sp3_file << fixed << "#c" << posvel << setw(4) << year  << " " << setw(2) << month << " " << setw(2) << day << " " << setfill('0') << setprecision(8) << setw(2)  << hour << " " << setw(2)  << min << " " << setw(11) << sec << " " << setfill(' ') << setw(7) << OrbEph_rows << " " << DataType << " " << CoordSys << " " << OrbType << " " << AgencyName << endl;
  
    // 2nd header line
    sp3_file << "##" << setw(5) << week << setprecision(8) << setw(16) << weeksecs << setprecision(8) << setw(15) << OrbEph_step << setw(6) << MJD0 << setprecision(13) << setw(16) << fracOfDay << endl;
  
    // Header lines 3 to 7 (List of satellites)
    sp3_file << "+ " << setw(4) << 1 << "   " << SatName << "  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0" << endl;
    sp3_file << "+          0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0" << endl;
    sp3_file << "+          0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0" << endl;
    sp3_file << "+          0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0" << endl;
    sp3_file << "+          0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0" << endl;
  
    // Header lines 8 to 12
    sp3_file  << "++         0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0" << endl;
    sp3_file  << "++         0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0" << endl;
    sp3_file  << "++         0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0" << endl;
    sp3_file  << "++         0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0" << endl;
    sp3_file  << "++         0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0" << endl;
    
    // Header lines 13 to 18
    sp3_file  << "%c L  cc GPS ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc" << endl;
    sp3_file  << "%c cc cc ccc ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc" << endl;
    sp3_file  << "%f  0.0000000  0.000000000  0.00000000000  0.000000000000000" << endl;
    sp3_file  << "%f  0.0000000  0.000000000  0.00000000000  0.000000000000000" << endl;
    sp3_file  << "%i    0    0    0    0      0      0      0      0         0" << endl;
    sp3_file  << "%i    0    0    0    0      0      0      0      0         0" << endl;
  
    // Header comment lines 19 to 22
    sp3_file << "/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC" << endl;
    sp3_file << "/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC" << endl;
    sp3_file << "/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC" << endl;
    sp3_file << "/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC" << endl;
    
    ////////////////////////////// Write sp3 file data section //////////////////////////////
    for(int i = 0; i < OrbEph_rows; i++)
      {
      int width;
      //Epoch.Date(year,month,day);
      GPS_UTC = leapsec(OrbEph_epochs(0)) - timescales::TAI_GPS;
      Date = GPS2UTCdate(OrbEph_epochs(i) + GPS_UTC, "ISOC");
      
      year = (int)Date(0);
      month = (int)Date(1);
      day = (int)Date(2);
      hour = (int)Date(3);
      min = (int)Date(4);
      sec = round(Date(5));
      
      posECEF = loaded_OrbEph.row(i).segment(8,3);
      velECEF = loaded_OrbEph.row(i).segment(11,3);
    
      // Epoch
      sp3_file << "*  " << setfill('0') << setw(4) << year  << " " << setw(2) << month << " " << setw(2) << day << " " << setw(2) << hour << " " << setw(2) << min << " " << setprecision(8) << setw(11) << sec << setfill(' ') << endl;  
        
      if (posECEF.norm() < 10000.0e3 ) width = 7;
      else width = 6;
      
      // State
      sp3_file << "P" + SatName << setprecision(width) << setw(14) << posECEF(0)/1000.0 << setprecision(width) << setw(14) << posECEF(1)/1000.0 << setprecision(width) << setw(14) << posECEF(2)/1000.0 << setprecision(6) << setw(14) << 999999.999999 << endl;
      if( posvel.compare("V") == 0 ) sp3_file << "V" + SatName << setprecision(width) << setw(14) << velECEF(0)*10 << setprecision(width) << setw(14) << velECEF(1)*10 << setprecision(width) << setw(14) << velECEF(2)*10 << setprecision(6) << setw(14) << 999999.999999 << endl;
      
      //Epoch += OrbEph_step;
      }
    
    // Close file        
    sp3_file.close();
    
    return(0);
    } // End of main()
