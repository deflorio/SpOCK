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

#ifndef IO_UTILS_H_
#define IO_UTILS_H_

#include <string>

#include <VarTypes.h>

#include <boost/tokenizer.hpp>

#ifdef USE_MATLAB

	#include "/opt/MATLAB/R2017b/extern/include/mat.h"
	
#endif
  
using namespace std;
using namespace SC;

// Load and read a csv file
Eigen::MatrixXd read_csvfile(const char* filename, int cols, bool header = false);
// Determine if a string is a number
bool isNumber(const string& str);
// Load and read a geopotential model file
Matrix3D read_gfc(const char* filename, int maxdeg, double epoch, double& mu, double& Re);
// Load and read a space weather indices file
Eigen::MatrixXf read_SpaceWeather(const char* filename, double start_epoch, int sim_duration);

// Write data to a Matlab mat file
//void write_matfile(const Eigen::VectorXd& ephem, int cols, const char* matvarname, const char* filename);
//#ifdef USE_MATLAB
//
//int Vec2matfile(const char *filename, const Eigen::VectorXd& vec_in);
//
//#endif

// Load and read simulation parameters XML file and validate it against an XML schema file
int XML_parser(const string XML_simparam_file, string& Orbit_ephemeris, string& Attitude_ephemeris, string& TLE_file, string& Data_path, string& planetephemeris, string& eop, string& pck_data, string& leapsecond, string& magneticfield, string& gravityfield, string& atmosphere, string& sunmoon, string& orbfile_name, string& attfile_name, string& sensors_filename, string& csv_torques_name, string& csv_accelerations_name, int& SIM_STEP, int& SIM_DURATION, Vector6d& init_orbtime, Vector6d& init_orbstate, double& phi, double& theta, double& psi, double& om_x, double& om_y, double& om_z, bool& initstate_in_RTN, bool& realtime, double& realtime_wait, bool& ggrad_on, bool& mag_on, bool& drag_on, bool& srp_on, int& nMAX, bool& sunmoon_on, string& Drag_Model, string& SRP_Model, string& AttitudeType, bool& attctrl_on, string& AttCtrlType, bool& orbctrl_on, string& OrbCtrlType, double& SC_mass, Mat3x3d& MoI, Vec3d& CoG, double& SC_Cd, double& SC_Cr, double& SC_Area_D, double& SC_Area_R, Vec3d& Mdip, Face& F_Xplus, Face& F_Xminus, Face& F_Yplus, Face& F_Yminus, Face& F_Zplus, Face& F_Zminus, SYS_params& Sensor_prm_SUN, SYS_params& Sensor_prm_EARTH, SYS_params& Sensor_prm_CSS1, SYS_params& Sensor_prm_CSS2, SYS_params& Sensor_prm_CSS3, SYS_params& Sensor_prm_CSS4, SYS_params& Sensor_prm_CSS5, SYS_params& Sensor_prm_CSS6, SYS_params& Sensor_prm_MAG, SYS_params& Sensor_prm_MAGstowed, SYS_params& Sensor_prm_RS, SYS_params& Sensor_prm_MAGTRQ, SYS_params& Sensor_prm_WHEEL1, SYS_params& Sensor_prm_WHEEL2, SYS_params& Sensor_prm_WHEEL3, SYS_params& Solarpan1_prm, SYS_params& Solarpan2_prm, SYS_params& Solarpan3_prm, SYS_params& OrbitPropulsion1_prm, SYS_params& OrbitPropulsion2_prm, vector<maneuver>& all_maneuvers);

// Load and read simulation parameters XML file and validate it against an XML schema file
int SGP4_XML_parser(const string XML_simparam_file, string& TLE_file, string& Data_path, string& eop, string& pck_data, string& leapsecond, string& orbfile_name, int& SIM_STEP, int& SIM_DURATION);

// Write simulation parameters read by parser XML_parser in a text file for verification purposes
void ReadXMLtoTXT(const string txt_file, string Orbit_ephemeris, string Attitude_ephemeris, string TLE_file, string Data_path, string planetephemeris, string eop, string pck_data, string leapsecond, string magneticfield, string gravityfield, string atmosphere, string sunmoon, string orbfile_name, string attfile_name, string sensors_filename, string csv_torques_name, string csv_accelerations_name, int SIM_STEP, int SIM_DURATION, Vector6d init_orbtime, Vector6d init_orbstate, double phi, double theta, double psi, double om_x, double om_y, double om_z, bool initstate_in_RTN, bool realtime, double realtime_wait, bool ggrad_on, bool mag_on, bool drag_on, bool srp_on, int nMAX, bool sunmoon_on, string Drag_Model, string SRP_Model, string AttitudeType, bool attctrl_on, string AttCtrlType, bool orbctrl_on, string OrbCtrlType, double SC_mass, Mat3x3d MoI, Vec3d CoG, double SC_Cd, double SC_Cr, double SC_Area_D, double SC_Area_R, Vec3d Mdip, Face F_Xplus, Face F_Xminus, Face F_Yplus, Face F_Yminus, Face F_Zplus, Face F_Zminus, SYS_params Sensor_prm_SUN, SYS_params Sensor_prm_EARTH, SYS_params Sensor_prm_CSS1, SYS_params Sensor_prm_CSS2, SYS_params Sensor_prm_CSS3, SYS_params Sensor_prm_CSS4, SYS_params Sensor_prm_CSS5, SYS_params Sensor_prm_CSS6, SYS_params Sensor_prm_MAG, SYS_params Sensor_prm_MAGstowed, SYS_params Sensor_prm_RS, SYS_params Sensor_prm_MAGTRQ, SYS_params Sensor_prm_WHEEL1, SYS_params Sensor_prm_WHEEL2, SYS_params Sensor_prm_WHEEL3, SYS_params Solarpan1_prm, SYS_params Solarpan2_prm, SYS_params Solarpan3_prm, SYS_params OrbitPropulsion1_prm, SYS_params OrbitPropulsion2_prm, vector<maneuver> all_maneuvers);

// Load and read XML file for events computation and validate it against an XML schema file
int XML_parser_events(const string XML_events_file, int& simstep, int& duration, double& FOV_cross, double& FOV_along, int& SC_start, int& SC_end, int& PL_start, int& PL_end, bool& TGs_on, bool& GSs_on, bool& TGs_grid_on, bool& Eclipse_on, Vec4d& TG_grid_limits, double& gridstep, ground::TG* TGs_list, ground::GS* GSs_list, string& Orbit_ephemeris_path, string& Orbit_ephemeris_rootname, string& Data_path, string& planetephemeris, string& eop, string& pck_data, string& leapsecond, string& TG_filename, string& GS_filename, string& Eclipse_filename);
      
// Write events computation parameters read by parser XML_parser_events in a text file for verification purposes
void ReadXMLeventstoTXT(const string txt_file, int simstep, int duration, double FOV_cross, double FOV_along, int SC_start, int SC_end, int PL_start, int PL_end, bool TGs_on, bool GSs_on, bool TGs_grid_on, bool Eclipse_on, Vec4d TG_grid_limits, double gridstep, ground::TG* TGs_list, ground::GS* GSs_list, string Orbit_ephemeris_path, string Orbit_ephemeris_rootname, string Data_path, string planetephemeris, string eop, string pck_data, string leapsecond, string TG_filename, string GS_filename, string Eclipse_filename);

// Run-start display message
void RunStartMessage(Vector6d init_orbtime, Vector6d init_orbstate, int SIM_DURATION, bool* T_model, int nMAX, string Drag_Model, string SRP_Model, string magneticfield, string gravityfield, string atmosphere, string sunmoon, string proptype);

// Display in terminal the simulation execution status bar
void RunStatusBar(double t, int simduration, int barwidth);


#endif // IO_UTILS_H_