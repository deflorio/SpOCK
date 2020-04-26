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

#include <chrono>
#include <thread>

#include <IO_utils.h>
#include <VarTypes.h>
#include <Constants.h>
#include <Transformations.h>
#include <Solarsys.h>
#include <Events.h>
// External libraries: Eigen
#include <Eigen/Core>

#ifdef USE_SPICE

extern "C"
      {
      #include <SpiceUsr.h>
      }
      
#endif

  
      
using namespace std;
using namespace SC;
using namespace solarsystem;
using namespace constants;
using namespace mathconst;
using namespace math;

using namespace Eigen;

////////////////////////// Initialization ///////////////////////////

int main(int argc, char *argv[])
    {
    chrono::time_point<chrono::high_resolution_clock> clockstart, clockend;
    
    clockstart = chrono::high_resolution_clock::now();
    
    /////////////// Simmulation parameters XML file ///////////////
    if(argc < 2)
      {
	  // Tell the user how to run the program
	  cerr << "Usage: " << argv[0] << " path_to/simulation/parameters/file.xml\nEx: ./bin/EventsComputation /home/username/path1/path2/input/simparam.xml" << endl;
	  return 1;
	  }
	  
	string XML_events_file(argv[1]);
    // Step of input ephemerides
    int simstep;
    // Time span considered for the events computation
    int duration;
    // Field of view half angle perpendicolar to the direction of motion
    double FOV_cross;
    // Field of view half angle along the direction of motion
    double FOV_along;
    // Spacecraft number to start with
    int SC_start;
    // Spacecraft number to end with
    int SC_end;
    // Orbital plane number to start with
    int PL_start;
    // Orbital plane number to end with
    int PL_end;
    // Compute contacts of payload with targets
    bool TGs_on;
    // Compute contacts of spacecraft with ground stations
    bool GSs_on;
    // Use targets grid for computation of contacts of payload with targets
    bool TGs_grid_on;
    // Compute eclipse times
    bool Eclipse_on;
    // Minimum and maximum longitute and latitude of targets grid
    Vec4d TG_grid_limits;
    // Step of the targets grid
    double gridstep;
	// Targets list
	ground::TG TGs_list[1000];
    // Ground stations list
	ground::GS GSs_list[1000];
	// Input files path
    string ephem_file, Orbit_ephemeris_path, Orbit_ephemeris_rootname, Data_path, planetephemeris, eop, pck_data, leapsecond;
    // Output files path
    string TG_filename, GS_filename, Eclipse_filename;
    
    XML_parser_events(XML_events_file, simstep, duration, FOV_cross, FOV_along, SC_start, SC_end, PL_start, PL_end, TGs_on, GSs_on, TGs_grid_on, Eclipse_on, TG_grid_limits, gridstep, TGs_list, GSs_list, Orbit_ephemeris_path, Orbit_ephemeris_rootname, Data_path, planetephemeris, eop, pck_data, leapsecond, TG_filename, GS_filename, Eclipse_filename);
    //const string ReadXML_TXT_file = "../input/readXMLevents.txt";
    size_t last_slash = XML_events_file.find_last_of("/");
    string ReadXML_TXT_file_name = XML_events_file.substr(last_slash+1);
    size_t lastspoint = ReadXML_TXT_file_name.find_last_of(".");
    ReadXML_TXT_file_name = ReadXML_TXT_file_name.substr(0,lastspoint);
    const string ReadXML_TXT_file = XML_events_file.substr(0,last_slash) + "/Read_" + ReadXML_TXT_file_name + ".txt"; 
    // Put read XML in a text file (for check purposes)
    ReadXMLeventstoTXT(ReadXML_TXT_file, simstep, duration, FOV_cross, FOV_along, SC_start, SC_end, PL_start, PL_end, TGs_on, GSs_on, TGs_grid_on, Eclipse_on, TG_grid_limits, gridstep, TGs_list, GSs_list, Orbit_ephemeris_path, Orbit_ephemeris_rootname, Data_path, planetephemeris, eop, pck_data, leapsecond, TG_filename, GS_filename, Eclipse_filename);
	
    string TG_output_path, GS_output_path, Umbras_output_path;
    
    last_slash = TG_filename.rfind("/");
    TG_output_path = TG_filename;
    TG_output_path.erase(last_slash, TG_filename.npos);
    
    last_slash = GS_filename.rfind("/");
	GS_output_path = GS_filename;
    GS_output_path.erase(last_slash, GS_filename.npos);
	
	last_slash = Eclipse_filename.rfind("/");
	Umbras_output_path = Eclipse_filename;
    Umbras_output_path.erase(last_slash, Eclipse_filename.npos);
    
    double grid_minlon, grid_maxlon, grid_minlat, grid_maxlat;
    
    grid_minlon = TG_grid_limits(0);
    grid_maxlon = TG_grid_limits(1);
    grid_minlat = TG_grid_limits(2);
    grid_maxlat = TG_grid_limits(3);
    
    // SPICE Kernels            
    #ifdef USE_SPICE
    
        planetephemeris = Data_path + "/cspice/" + planetephemeris;
        eop = Data_path + "/cspice/" + eop;
        leapsecond = Data_path + "/cspice/" + leapsecond;
        pck_data = Data_path + "/cspice/" + pck_data;
    
        furnsh_c(planetephemeris.c_str( ));
        furnsh_c(eop.c_str( ));
        furnsh_c(leapsecond.c_str( ));
        furnsh_c(pck_data.c_str( ));
          
    #endif
    
    /////////////////////////////////////////////////////////////////////////
    ////////////////////////// ANALYSIS PARAMETERS //////////////////////////
    /////////////////////////////////////////////////////////////////////////
    
    ///////////////////
    // Targets grid //
    //////////////////
	
	VectorXd TGs_grid_lons, TGs_grid_lats;
	
	if(TGs_grid_on)
	  {
	  if(grid_maxlon < grid_minlon)
		{
		// Tell the user how to run the program
		std::cerr << "Maximum targets grid longitude has to be larger than minimum grid longitude" << std::endl;
		return 1;
		}
	  
	  if(grid_maxlat < grid_minlat)
		{
		// Tell the user how to run the program
		std::cerr << "Maximum targets grid latitude has to be larger than minimum grid latitude" << std::endl;
		return 1;
		}
	  
	  if( (fabs(grid_maxlon-grid_minlon)/gridstep < 1) && (grid_minlon != grid_maxlon) )
		{
		// Tell the user how to run the program
		std::cerr << "The grid step has to be smaller than the difference between the maximum and minimum longitudes of the targets grid" << std::endl;
		return 1;
		}
	  
	  if( (fabs(grid_maxlat-grid_minlat)/gridstep < 1) && (grid_minlat != grid_maxlat) )
		{
		// Tell the user how to run the program
		std::cerr << "The grid step has to be smaller than the difference between the maximum and minimum latitudes of the targets grid" << std::endl;
		return 1;
		}
	  
	  int grid_rows = fabs(grid_maxlat-grid_minlat)/gridstep + 1;
	  int grid_cols = fabs(grid_maxlon-grid_minlon)/gridstep + 1;
    
	  TGs_grid_lons.resize(grid_cols);
	  TGs_grid_lats.resize(grid_rows);
    
	  double lon = grid_minlon;
	  double lat = grid_minlat;
    
	  for(int i = 0 ; i < grid_rows; i++)
        {
		TGs_grid_lats(i) = lat;  
		lat = lat + gridstep;
		}
		
	  for(int i = 0 ; i < grid_cols; i++)
        {
		TGs_grid_lons(i) = lon;  
		lon = lon + gridstep;
		}
	  }
		
    ////////////////////////////////////////////////////////////////
    /////////////// FILE WITH ALL CONTACTS INFORMATION /////////////
    ////////////////////////////////////////////////////////////////    
    
    ofstream AllContacts_file;
    if(TGs_on || TGs_grid_on)
        {
        AllContacts_file.open(TG_filename);
        
        AllContacts_file << "TG,Epoch UTC in,GPS time [s] in,Duration [s],El. in,El. out,Epoch UTC Max. El.,Max. El. [deg],Az. in [deg],Az. out [deg],Az. Max El. [deg],lon [deg],lat [deg],Orbital plane,Spacecraft" << endl;
        }
	  
	ofstream GS_AllContacts_file;
    if(GSs_on)
        {
        GS_AllContacts_file.open(GS_filename);
        
        GS_AllContacts_file << "GS,AOS UTC,LOS UTC,AOS [GPS secs],Duration [m],Epoch UTC Max. El.,Max Elevation [deg],Az. AOS [deg],Az. LOS [deg],Az. Max El. [deg],lon [deg],lat [deg],Orbital plane,Spacecraft" << endl;
        }
      
    SOLSYS Suncompute;
    ofstream AllUmbras_file;
	
    if(Eclipse_on)
        {
        AllUmbras_file.open(Eclipse_filename);
        
        AllUmbras_file << "Penumbra start UTC,Umbra start UTC,Umbra end UTC,Penumbra end UTC,Umbra duration [m],Penumbra duration [m],Orbital plane,Spacecraft" << endl;
        }
        
    if(TGs_on)  cout << "Compute TG contacts (list)" << endl;
    if(TGs_grid_on) cout << "Compute TG contacts (grid)" << endl;
    if(GSs_on)  cout << "Compute GS contacts" << endl;
    if(Eclipse_on)  cout << "Compute umbra/penumbra times" << endl;
    
    //////////////////////////////////////////////////////////////////////
    ////////////////////////// COMPUTE CONTACTS //////////////////////////
    //////////////////////////////////////////////////////////////////////
    
    MatrixXd loaded_ephem;
    VectorNd<15> ephem_row = VectorNd<15>::Zero();
    Vector6d ECIstate = Vector6d::Zero();
    Vector6d orbel = Vector6d::Zero();
    Vector6d mean_orbel = Vector6d::Zero();
    double GPStime;
    double ephem_duration;
    
    int matrows = 0;
    
    //#pragma omp parallel for 
    for(int p_ind = PL_start ; p_ind <= PL_end; p_ind++)
            {
            //#pragma omp for
            for(int s_ind = SC_start ; s_ind <= SC_end; s_ind++)
                {
                // Load spacecraft ephemeris file
                bool valid;
                Vec4d a_ex_ey_i;
                ephem_file = Orbit_ephemeris_path + "/" + "S" + to_string(s_ind) + "-P" + to_string(p_ind) + "_" + Orbit_ephemeris_rootname;
				if(PL_start == 1 || PL_end == 1 || SC_start == 1 || SC_end == 1) ephem_file = Orbit_ephemeris_path + "/" + Orbit_ephemeris_rootname;
                //cout << ephem_file << endl;
                loaded_ephem = read_csvfile(ephem_file.c_str(),15);
                matrows = loaded_ephem.rows();
                //cout << matrows << endl;
                int indECI = round(matrows/2); // A state placed in the middle of ephemeris file
                ephem_row = loaded_ephem.row(indECI);
                
                ECIstate = ephem_row.segment(2,6);
                orbel = rv2oe(ECIstate, valid); // Osculating elements
                mean_orbel = osc2mean(orbel, valid);
                a_ex_ey_i = mean_orbel.segment(0,4);
                
                //Matrix<double,Dynamic,7> orbstateECEF;
                MatrixXd orbstateECEF = loaded_ephem.block(0,7,matrows,7); // Column 0 is ECI v_Z to be replaced by GPS time
                MatrixXd time_posECI = loaded_ephem.block(0,1,matrows,4); // Column 0 is GPS secs
                
                ephem_duration = loaded_ephem(matrows-1,1) - loaded_ephem(0,1);
                if(duration > ephem_duration)
                    {
                    cerr << "\nThe duration of the computation (<simduration>) cannot be longer than the duration of the input orbit ephemeris" << endl;
                    exit(EXIT_FAILURE);
                    }
                
                for(int i = 0 ; i < matrows; i++)
                    {
                    ephem_row = loaded_ephem.row(i);
                    
                    GPStime = ephem_row(1);
                    orbstateECEF(i,0) = GPStime;
                    }
                    
                if(TGs_on || TGs_grid_on) TGsContacts(s_ind, p_ind, a_ex_ey_i, orbstateECEF, TGs_grid_lons, TGs_grid_lats, TGs_list, TG_output_path, FOV_cross, FOV_along, duration, simstep, TGs_on, TGs_grid_on, AllContacts_file);
                  
                if(GSs_on) GSsContacts(s_ind, p_ind, a_ex_ey_i, orbstateECEF, GSs_list, GS_output_path, duration, simstep, GS_AllContacts_file);
				
				if(Eclipse_on) Umbras(s_ind, p_ind, Suncompute, time_posECI, Umbras_output_path, duration, simstep, AllUmbras_file);
                }
            }
            
    //////////////////////////////////////////////////////////////////////
    /////////////// WRITE FILE WITH ALL CONTACTS INFORMATION /////////////
    //////////////////////////////////////////////////////////////////////        
    
    if(TGs_on || TGs_grid_on) AllContacts_file.close();
	if(GSs_on) GS_AllContacts_file.close();
	
	if(Eclipse_on) AllUmbras_file.close();
    
    clockend = chrono::high_resolution_clock::now();
        
    chrono::duration<double,milli> elapsed_millisecs = clockend - clockstart;
    cout << "Elapsed seconds: " << elapsed_millisecs.count()/1000.0 << endl;
    
    return(0);
    }
