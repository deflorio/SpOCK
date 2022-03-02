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

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <IO_utils.h>
#include <Transformations.h>
#include <Eigen/Core>
#include <Eigen/Geometry> 

#ifdef _WIN32
#include <io.h>
#endif

#include <sgp4ext.h>
#include <sgp4unit.h>
#include <sgp4io.h>

using namespace std;
using namespace Eigen;

//int main()
int main(int argc, char *argv[])
	{
	/////////////// Simmulation parameters XML file ///////////////
	if(argc < 2)
		{
	// Tell the user how to run the program
	cerr << "Usage: " << argv[0] << " path_to/simulation/parameters/file.xml\nEx: ./bin/OrbitPropagator /home/username/path1/path2/input/simparam.xml" << endl;
	return 1;
	}
	
	string XML_simparam_file(argv[1]);
		
	//	if(argc < 5)
	//      {
	//	  // Tell the user how to run the program
	//	  cerr << "Usage: " << argv[0] << " path_to_TLEs/TLE_file.tle" << endl;
	//	  return 1;
	//	  }
	//	  
	//	double step, duration;
	//	
	//	step = atof(argv[1]); // [s]
	//	duration = atof(argv[2]); // [d]
	//	
	//	string TLE_file(argv[3]);
	//	string ephemname(argv[4]);
	
	// Input files path
	string TLE_file, Data_path, eop, pck_data, leapsecond;
	// Output files path
	string ephemname;
	// Simulation step
	int SIM_STEP;
	// Simulation duration
	int SIM_DURATION;
	
	SGP4_XML_parser(XML_simparam_file, TLE_file, Data_path, eop, pck_data, leapsecond, ephemname, SIM_STEP, SIM_DURATION);
	
	size_t charnum_TLEname = TLE_file.size() + 10;
	
	char* fgetsout;
	char str[2];
	char infilename[charnum_TLEname];
	double ro[3];
	double vo[3];
	
	Vector6d UTCdate;
	Vector6d ECI_state;
	Vector6d ECEF_state;
	
	string sec_str;
	//string UTC_epoch;
	int propsecs = 0;
	int init_GPSsecs = 0;
	int GPSsecs = 0;
	char typerun, typeinput, opsmode;
	gravconsttype whichconst;
	//int whichcon;
	ofstream orbstate_file;
	
	//FILE *infile, *outfile, *outfilee;
	FILE *infile;
	
	// Load SPICE Kernels            
	#ifdef USE_SPICE
	
			eop = Data_path + "/cspice/" + eop;
			leapsecond = Data_path + "/cspice/" + leapsecond;
			pck_data = Data_path + "/cspice/" + pck_data;
	
			furnsh_c(eop.c_str( ));
			furnsh_c(leapsecond.c_str( ));
			furnsh_c(pck_data.c_str( ));
				
	#endif

	// ----------------------------  Locals  -------------------------------
  double sec,  jd, startmfe, stopmfe, deltamin;
	double tsince = 0.0;
  double tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2;
	int  year; int mon; int day; int hr; int min;
	char longstr1[130];
	typedef char str3[4];
	str3 monstr[13];
	char longstr2[130];
	elsetrec satrec;
	// ------------------------  Implementation   --------------------------
	strcpy(monstr[1], "Jan");
	strcpy(monstr[2], "Feb");
	strcpy(monstr[3], "Mar");
	strcpy(monstr[4], "Apr");
	strcpy(monstr[5], "May");
	strcpy(monstr[6], "Jun");
	strcpy(monstr[7], "Jul");
	strcpy(monstr[8], "Aug");
	strcpy(monstr[9], "Sep");
	strcpy(monstr[10], "Oct");
	strcpy(monstr[11], "Nov");
	strcpy(monstr[12], "Dec");

	startmfe = 0.0;
	stopmfe = SIM_DURATION/60.0; // [min] 10.0*1440.0; // [min]
	deltamin = SIM_STEP/60.0; // [min] 1.0/6.0; // 10 s

	printf("\n%s\n",SGP4Version );

			//opsmode = 'a' best understanding of how afspc code works
			//opsmode = 'i' imporved sgp4 resulting in smoother behavior
			//printf("input operation mode a, i \n\n");
			//opsmode = getchar();
			//fflush(stdin);
	
	opsmode = 'i';

			//typerun = 'c' compare 1 year of full satcat data
			//typerun = 'v' verification run, requires modified elm file with
			//              start, stop, and delta times
			//typerun = 'm' maunual operation- either mfe, epoch, or dayof yr also
			//printf("input type of run c, v, m \n\n");
			//typerun = getchar();
			//fflush(stdin);
	
	typerun = 'm';

			//typeinput = 'm' input start stop mfe
			//typeinput = 'e' input start stop ymd hms
			//typeinput = 'd' input start stop yr dayofyr
			//if ((typerun != 'v') && (typerun != 'c'))
			//  {
			//    printf("input mfe, epoch (YMDHMS), or dayofyr approach, m,e,d \n\n");
			//    typeinput = getchar();
			//  }
			//  else
			//    typeinput = 'e';
		
	typeinput = 'm';

			//printf("input which constants 721 72 84 \n");
			//scanf( "%i",&whichcon );
			//if (whichcon == 721) whichconst = wgs72old;
			//if (whichcon == 72) whichconst = wgs72;
			//if (whichcon == 84) whichconst = wgs84;
	
	whichconst = wgs84;

	getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );

	// ---------------- Setup files for operation ------------------
	strcpy(infilename, TLE_file.c_str());
	
	infile = fopen(infilename, "r");
	if(infile == NULL)
		{
		printf("Failed to open file: %s\n", infilename);
		return 1;
		}
	
	while (feof(infile) == 0)
		{
		do
			{
			fgetsout = fgets( longstr1,130,infile);
			if(fgetsout == NULL){};
				strncpy(str, &longstr1[0], 1);
				str[1] = '\0';
			} while ((strcmp(str, "#")==0)&&(feof(infile) == 0));

		if (feof(infile) == 0)
			{
			fgetsout = fgets( longstr2,130,infile);
			if(fgetsout == NULL){};
			// Convert the char string to sgp4 elements
			// Includes initialization of sgp4
			twoline2rv(longstr1, longstr2, typerun, typeinput, opsmode, whichconst, startmfe, stopmfe, deltamin, satrec);
			//ephemname = "SGP4_Ephemerides/SGP4_ephemeris_ECI-ECEF_" + to_string(satrec.satnum) + ".csv";
			orbstate_file.open(ephemname);

			cout << "TLE NORAD ID " + to_string(satrec.satnum) + "\n";
			// call the propagator to get the initial state vector value
			sgp4 (whichconst, satrec,  0.0, ro,  vo);
			
			ECI_state(0) = ro[0]*1E3; // [m]
			ECI_state(1) = ro[1]*1E3; // [m]
			ECI_state(2) = ro[2]*1E3; // [m]
			ECI_state(3) = vo[0]*1E3; // [m]
			ECI_state(4) = vo[1]*1E3; // [m]
			ECI_state(5) = vo[2]*1E3; // [m]
			
			jd = satrec.jdsatepoch + tsince/1440.0;
			invjday(jd,year,mon,day,hr,min,sec);
			
			UTCdate << year, mon, day, hr, min, sec;
			GPSsecs = UTCdate2GPSsecs(UTCdate);
			init_GPSsecs = GPSsecs;
			
			propsecs = GPSsecs - init_GPSsecs;
			
			ECEF_state = ECI2ECEF(GPSsecs,ECI_state);
			
			//stringstream sec_stream;
			//sec_stream << setfill('0') << setw(2) << round(sec);
			//sec_str = sec_stream.str();
			//UTC_epoch = to_string(year) + "-" + to_string(mon) + "-" + to_string(day) + " " + to_string(hr) + ":" + to_string(min) + ":" + sec_str;
			
			orbstate_file << propsecs << "," << fixed << GPSsecs << "," << ECI_state(0) << "," << ECI_state(1) << "," << ECI_state(2) << "," << ECI_state(3) << "," << ECI_state(4) << "," << ECI_state(5) << "," << ECEF_state(0) << "," << ECEF_state(1) << "," << ECEF_state(2) << "," << ECEF_state(3) << "," << ECEF_state(4) << "," << ECEF_state(5) << endl;
		
			tsince = startmfe;
			// check so the first value isn't written twice
			if ( fabs(tsince) > 1.0e-8 )
					tsince = tsince - deltamin;
		
			// ----------------- Loop to perform the propagation ----------------
			while ((tsince < stopmfe) && (satrec.error == 0))
				{
				tsince = tsince + deltamin;
		
				if(tsince > stopmfe) tsince = stopmfe;
		
				sgp4(whichconst, satrec,  tsince, ro,  vo);
		
				if (satrec.error > 0) printf("# *** error: t:= %f *** code = %3d\n", satrec.t, satrec.error);
		
				if (satrec.error == 0)
					{
					jd = satrec.jdsatepoch + tsince/1440.0;
					invjday( jd, year,mon,day,hr,min, sec );
			
					ECI_state(0) = ro[0]*1E3; // [m]
					ECI_state(1) = ro[1]*1E3; // [m]
					ECI_state(2) = ro[2]*1E3; // [m]
					ECI_state(3) = vo[0]*1E3; // [m]
					ECI_state(4) = vo[1]*1E3; // [m]
					ECI_state(5) = vo[2]*1E3; // [m]
			
					UTCdate << year, mon, day, hr, min, sec;
					GPSsecs = UTCdate2GPSsecs(UTCdate);
					propsecs = GPSsecs - init_GPSsecs;
					
					ECEF_state = ECI2ECEF(GPSsecs,ECI_state);
					
					//stringstream sec_stream;
					//sec_stream << setfill('0') << setw(2) << round(sec);
					//sec_str = sec_stream.str();
					//UTC_epoch = to_string(year) + "-" + to_string(mon) + "-" + to_string(day) + " " + to_string(hr) + ":" + to_string(min) + ":" + sec_str;
					
					orbstate_file << propsecs << "," << fixed << GPSsecs << "," << ECI_state(0) << "," << ECI_state(1) << "," << ECI_state(2) << "," << ECI_state(3) << "," << ECI_state(4) << "," << ECI_state(5) << "," << ECEF_state(0) << "," << ECEF_state(1) << "," << ECEF_state(2) << "," << ECEF_state(3) << "," << ECEF_state(4) << "," << ECEF_state(5) << endl;
					} // if satrec.error == 0
				} // while propagating the orbit
		
			//fprintf(outfilee," END Ephemeris \n");
			//fclose (outfilee);
			
			orbstate_file.close();
			} // if not eof
    } // while through the input file

  return 0;
	}  // end testcpp
