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

#include <typeinfo>

#include <chrono>
#include <thread>

#include <Orbit.h>
#include <IO_utils.h>
#include <VarTypes.h>
#include <Constants.h>
#include <Transformations.h>
#include <HIL_interface.h>

#include <Orbpropulsion.h>
#include <EarthSunSensor.h>
#include <Magneto.h>
#include <WheelMR.h>
#include <SolarPanel.h>

// External libraries: Eigen
#include <Eigen/Core>
#include "boost/multi_array.hpp"

#ifdef USE_SPICE

extern "C"
        {
        #include "extlib/cspice/include/SpiceUsr.h"
        }
      
#endif

using namespace propagator;
using namespace orbit;

using namespace std;
using namespace SC;
using namespace constants;
using namespace mathconst;
using namespace math;

using namespace earthsun;
using namespace magneto;
using namespace mrwheel;
using namespace solarpan;
using namespace orbpropulsion;

using namespace boost;


int main(int argc, char *argv[])
    {
    chrono::time_point<chrono::high_resolution_clock> clockstart, clockend;
    
    clockstart = chrono::high_resolution_clock::now();
    
    /////////////// Simmulation parameters XML file ///////////////
    
    if(argc < 2)
      {
	  // Tell the user how to run the program
	  cerr << "Usage: " << argv[0] << " path_to/simulation/parameters/file.xml\nEx: ./bin/OrbitPropagator /home/username/path1/path2/input/simparam.xml" << endl;
	  return 1;
	  }
	  
	string XML_simparam_file(argv[1]);
    
    /////////////////////////////////////////////////////////////
    //////////// CREATE ORBIT PROPAGATOR CLASS //////////////////
    /////////////////////////////////////////////////////////////
    
    ORB SC_Orbit;
    
    ///////////////////////////////////////////////////////////////
    //////////// SIMULATION PARAMETERS VARIABLES //////////////////
    ///////////////////////////////////////////////////////////////
    
    ///////////////////////// Files paths /////////////////////////
    
    // Input files path
    string Orbit_ephemeris, Attitude_ephemeris, TLE_file, Data_path, planetephemeris, eop, pck_data, leapsecond, magneticfield, gravityfield, atmosphere, sunmoon;
    // Output files path
    string orbfile_name, attfile_name, sensors_filename, csv_torques_name, csv_accelerations_name;
    
    ///////////////////////// Simulation parameters /////////////////////////
    
    // Simulation step
    int SIM_STEP;
    // Simulation duration
    int SIM_DURATION;
    // Initial orbit UTC date and time
    Vector6d init_orbtime = Vector6d::Zero();
    // Initial orbit state
    Vector6d init_orbstate = Vector6d::Zero();
    // Initial attitude state
    double phi, theta, psi, om_x, om_y, om_z;
    // Initial state in RTN frame on/off
    bool initstate_in_RTN;
    // Real time simulation on/off
    bool realtime;
    // Step execution waiting time for hardware-in-the-loop simulations
    double realtime_wait;
    // Gravity gradient torque on/off
    bool ggrad_on;
    // Magnetic torque on/off
    bool mag_on;
    // Atmospheric drag torque on/off
    bool drag_on;
    // Solar radiation pressure torque on/off
    bool srp_on;
    // Maximum order and degree of gravitational field model used for the orbit propagation
    int nMAX;
    // Third body perturbation on/off
    bool sunmoon_on;
    // Atmospheric drag model used
    string Drag_Model;
    // Solar radiation pressure model used
    string SRP_Model;
    // Attitude type during orbit propagation
    string AttitudeType;
    // Attitude control on/off
    bool attctrl_on;
    // Attitude control type
    string AttCtrlType;
    // Orbit control on/off
    bool orbctrl_on;
    // Orbit control type
    string OrbCtrlType;
    
    ///////////////////////// Spacecraft properties /////////////////////////
    
    // Spacecraft mass
    double SC_mass;
    // Spacecraft center of mass position vector in body-fixed coordinates
    static Vec3d CoG;
    // Moments of inertia matrix. Moment of inertia taken at the center of mass and aligned with the body-fixed frame [kg*m^2]
    static Mat3x3d MoI = Mat3x3d::Zero();
    // Drag coefficient
    static double SC_Cd;
    // SRP coefficient
    static double SC_Cr;
    // Drag area to be used with atmospheric drag simple model
    static double SC_Area_D;
    // Radiation area to be used with solar radiation pressure simple model
    static double SC_Area_R;
    // Spacecraft magnetic dipole moment vector in body-fixed coordinates
    static Vec3d Mdip;
    // Spacecraft surfaces;
    Face F_Xplus, F_Xminus, F_Yplus, F_Yminus, F_Zplus, F_Zminus;
    F_Xplus.n = spacecraft::Normals.at("+X");
    F_Xminus.n = spacecraft::Normals.at("-X");
    F_Yplus.n = spacecraft::Normals.at("+Y");
    F_Yminus.n = spacecraft::Normals.at("-Y");
    F_Zplus.n = spacecraft::Normals.at("+Z");
    F_Zminus.n = spacecraft::Normals.at("-Z");
    
    ///////////////// ADCS sensors and actuators /////////////////
    
    // Sensors and actuators classes
    SYS_params Sensor_prm_SUN, Sensor_prm_EARTH, Sensor_prm_CSS1, Sensor_prm_CSS2, Sensor_prm_CSS3, Sensor_prm_CSS4, Sensor_prm_CSS5, Sensor_prm_CSS6, Sensor_prm_MAG, Sensor_prm_MAGstowed, Sensor_prm_RS, Sensor_prm_MAGTRQ, Sensor_prm_WHEEL1, Sensor_prm_WHEEL2, Sensor_prm_WHEEL3, Solarpan1_prm, Solarpan2_prm, Solarpan3_prm, OrbitPropulsion1_prm, OrbitPropulsion2_prm;
    ///////////////// Commanded attitude maneuvers /////////////////
    vector<maneuver> all_maneuvers; // Struct maneuver defined in VarTypes.h
    
    //////////////////////////////////////////////////////////////
    ////////// PARSING OF XML SIMULATION PARAMETERS FILE /////////
    //////////////////////////////////////////////////////////////
    
    XML_parser(XML_simparam_file, Orbit_ephemeris, Attitude_ephemeris, TLE_file, Data_path, planetephemeris, eop, pck_data, leapsecond, magneticfield, gravityfield, atmosphere, sunmoon, orbfile_name, attfile_name, sensors_filename, csv_torques_name, csv_accelerations_name, SIM_STEP, SIM_DURATION, init_orbtime, init_orbstate, phi, theta, psi, om_x, om_y, om_z, initstate_in_RTN, realtime, realtime_wait, ggrad_on, mag_on, drag_on, srp_on, nMAX, sunmoon_on, Drag_Model, SRP_Model, AttitudeType, attctrl_on, AttCtrlType, orbctrl_on, OrbCtrlType, SC_mass, MoI, CoG, SC_Cd, SC_Cr, SC_Area_D, SC_Area_R, Mdip, F_Xplus, F_Xminus, F_Yplus, F_Yminus, F_Zplus, F_Zminus, Sensor_prm_SUN, Sensor_prm_EARTH, Sensor_prm_CSS1, Sensor_prm_CSS2, Sensor_prm_CSS3, Sensor_prm_CSS4, Sensor_prm_CSS5, Sensor_prm_CSS6, Sensor_prm_MAG, Sensor_prm_MAGstowed, Sensor_prm_RS, Sensor_prm_MAGTRQ, Sensor_prm_WHEEL1, Sensor_prm_WHEEL2, Sensor_prm_WHEEL3, Solarpan1_prm, Solarpan2_prm, Solarpan3_prm, OrbitPropulsion1_prm, OrbitPropulsion2_prm, all_maneuvers);
    //cout << "Sono qui" << endl;
    size_t lastslash = XML_simparam_file.find_last_of("/");
    string ReadXML_TXT_file_name = XML_simparam_file.substr(lastslash+1);
    size_t lastspoint = ReadXML_TXT_file_name.find_last_of(".");
    ReadXML_TXT_file_name = ReadXML_TXT_file_name.substr(0,lastspoint);
    const string ReadXML_TXT_file = XML_simparam_file.substr(0,lastslash) + "/Read_" + ReadXML_TXT_file_name + ".txt"; 
    //const string ReadXML_TXT_file = "input/readXML.txt";
    // Put read XML in a text file (for check purposes)
    cout << ReadXML_TXT_file << endl;
    
    ReadXMLtoTXT(ReadXML_TXT_file, Orbit_ephemeris, Attitude_ephemeris, TLE_file, Data_path, planetephemeris, eop, pck_data, leapsecond, magneticfield, gravityfield, atmosphere, sunmoon, orbfile_name, attfile_name, sensors_filename, csv_torques_name, csv_accelerations_name, SIM_STEP, SIM_DURATION, init_orbtime, init_orbstate, phi, theta, psi, om_x, om_y, om_z, initstate_in_RTN, realtime, realtime_wait, ggrad_on, mag_on, drag_on, srp_on, nMAX, sunmoon_on, Drag_Model, SRP_Model, AttitudeType, attctrl_on, AttCtrlType, orbctrl_on, OrbCtrlType, SC_mass, MoI, CoG, SC_Cd, SC_Cr, SC_Area_D, SC_Area_R, Mdip, F_Xplus, F_Xminus, F_Yplus, F_Yminus, F_Zplus, F_Zminus, Sensor_prm_SUN, Sensor_prm_EARTH, Sensor_prm_CSS1, Sensor_prm_CSS2, Sensor_prm_CSS3, Sensor_prm_CSS4, Sensor_prm_CSS5, Sensor_prm_CSS6, Sensor_prm_MAG, Sensor_prm_MAGstowed, Sensor_prm_RS, Sensor_prm_MAGTRQ, Sensor_prm_WHEEL1, Sensor_prm_WHEEL2, Sensor_prm_WHEEL3, Solarpan1_prm, Solarpan2_prm, Solarpan3_prm, OrbitPropulsion1_prm, OrbitPropulsion2_prm, all_maneuvers);
    
    //////////////////////////////////////////////////////////////
    /////////////// PROCESSING OF PARSED VARIABLES ///////////////
    //////////////////////////////////////////////////////////////
    
    //////////////////// Orbit and environmental models /////////////////////
    
    // Load SPICE Kernels            
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
    
    EnvModels envmodels_paths;
    
    envmodels_paths.datapath = Data_path;
    envmodels_paths.sunmoon = sunmoon;
    envmodels_paths.magneticfield = magneticfield;
    envmodels_paths.gravityfield = gravityfield;
    envmodels_paths.atmosphere = atmosphere;
    
    ///////////////// Spacecraft parameters /////////////////
    
    SC_params SC_prms;
    
    SC_prms.SC_mass = SC_mass;
    SC_prms.MoI = MoI;
    SC_prms.Segment["+X"] = F_Xplus;
    SC_prms.Segment["-X"] = F_Xminus;
    SC_prms.Segment["+Y"] = F_Yplus;
    SC_prms.Segment["-Y"] = F_Yminus;
    SC_prms.Segment["+Z"] = F_Zplus;
    SC_prms.Segment["-Z"] = F_Zminus;
    SC_prms.Mdip = Mdip;
    SC_prms.CD = SC_Cd;
    SC_prms.C_SRP = SC_Cr;
    SC_prms.Area_D = SC_Area_D;
    SC_prms.Area_R = SC_Area_R;
    
    ///////////////////////// Attitude management //////////////////////
	
	Vec4d attstate = Vec4d::Zero();
    Mat3x3d ECItoBody, T_ECI2RTN, T_RTN2ECI, T_RTN2Body;
    
    VectorNd<7> attitudeRTN_state_vec = VectorNd<7>::Zero();
    Eigen::MatrixXd loaded_ephem;
    VectorNd<8> ephem_row = VectorNd<8>::Zero();
    int ind = 0;
    
    // Conversion to radians and normalization to 2*pi
    phi = mod(phi*DEG2RAD,PI2);
    theta = mod(theta*DEG2RAD,PI2);
    psi = mod(psi*DEG2RAD,PI2);
    
    if( AttitudeType.compare("RTN_fixed") == 0 )
      {
      T_ECI2RTN = ECI2RTN_Matrix(init_orbstate);
      T_RTN2ECI = T_ECI2RTN.transpose();
      T_RTN2Body = RotationMatrix321(phi, theta, psi);
      
      ECItoBody = T_RTN2Body*T_ECI2RTN;
  
      attstate = RotationMatrix2Quaternion(ECItoBody);
      }
    else if( AttitudeType.compare("Ephemeris") == 0 )
      {
      // Load orbit ephemerides
      try{ loaded_ephem = read_csvfile(Attitude_ephemeris.c_str(),8); }
      catch(const string errmsg)
            {
            cerr << "Attitude ephemerides: " + errmsg << endl;
            exit(EXIT_FAILURE);
            }
      }
	else
      {
      ECItoBody = RotationMatrix321(phi, theta, psi);
      attstate = RotationMatrix2Quaternion(ECItoBody);
      }
		
	//////////////////// Environment models //////////////////////////
    bool T_model[5] = {ggrad_on, mag_on, drag_on, srp_on, sunmoon_on};
    
    //////////////////// Run-start display message //////////////////////////
    RunStartMessage(init_orbtime, init_orbstate, SIM_DURATION, T_model, nMAX, Drag_Model, SRP_Model, magneticfield, gravityfield, atmosphere, sunmoon, "ORB");
    
    ///////////////// Orbit propulsion systems objects /////////////////
    ORBPROPULSION OrbitPropulsion1(OrbitPropulsion1_prm);
    OrbitPropulsion1.Init();
    OrbitPropulsion1.thrust2dv(SC_prms.SC_mass);
    
    ORBPROPULSION OrbitPropulsion2(OrbitPropulsion2_prm);
    OrbitPropulsion2.Init();
    OrbitPropulsion2.thrust2dv(SC_prms.SC_mass);
    
    ///////////////// Commanded orbital maneuvers /////////////////
    
    // Sort maneuvers by type and put them in vectors
    vector<maneuver> impulsive_maneuvers1; // Impulsive maneuvers
    vector<maneuver> impulsive_maneuvers2; // Impulsive maneuvers
    vector<maneuver> continuous_maneuvers1; // Continuous maneuvers executed by prop1
    vector<maneuver> continuous_maneuvers2; // Continuous maneuvers executed by prop2
    
    vector<maneuver> all_enabled_orbman;
    
    unsigned int man_ind;
    bool man_on;
    
    string man_name;
    
    for(man_ind = 0 ; man_ind < all_maneuvers.size(); man_ind++)
        {
        man_name = all_maneuvers[man_ind].name;
        man_on = all_maneuvers[man_ind].maneuver_on;
        
        if( (man_name.find("Prop") != string::npos) && man_on ) all_enabled_orbman.push_back(all_maneuvers[man_ind]);
        
        if( (man_name.find("ImpulsiveManeuver") != string::npos) && (man_name.find("Prop1") != string::npos) && man_on ) impulsive_maneuvers1.push_back(all_maneuvers[man_ind]);
        else if( (man_name.find("ImpulsiveManeuver") != string::npos) && (man_name.find("Prop2") != string::npos) && man_on ) impulsive_maneuvers2.push_back(all_maneuvers[man_ind]);
        else if( (man_name.find("ContinuousManeuver") != string::npos) && (man_name.find("Prop1") != string::npos) && man_on ) continuous_maneuvers1.push_back(all_maneuvers[man_ind]);
        else if( (man_name.find("ContinuousManeuver") != string::npos) && (man_name.find("Prop2") != string::npos) && man_on ) continuous_maneuvers2.push_back(all_maneuvers[man_ind]);
        }
        
    if( !continuous_maneuvers1.empty() ) OrbitPropulsion1.impman2contman(continuous_maneuvers1, SIM_STEP);
    if( !continuous_maneuvers2.empty() ) OrbitPropulsion2.impman2contman(continuous_maneuvers2, SIM_STEP);
    
    vector<maneuver>::iterator impman_ind1, impman_ind2, contman_ind1, contman_ind2; // To be used in propagation loop
    
    impman_ind1 = impulsive_maneuvers1.begin();
    impman_ind2 = impulsive_maneuvers2.begin();
    contman_ind1 = continuous_maneuvers1.begin();
    contman_ind2 = continuous_maneuvers2.begin();
    
    double m_init_time; // To be used in propagation loop
    
    Vec3d dv1_imp_CMD = Vec3d::Zero(); // To be used in propagation loop
    Vec3d dv1_cont_CMD = Vec3d::Zero();
    Vec3d dv2_imp_CMD = Vec3d::Zero();
    Vec3d dv2_cont_CMD = Vec3d::Zero();
    Vec3d dv_CMD = Vec3d::Zero();
    Vec3d dv_CTRL = Vec3d::Zero();
    
    string man_sys; // To be used in propagation loop
    
    double ini_GPSorbtime, GPStime;
    ini_GPSorbtime = UTCdate2GPSsecs(init_orbtime);
    
    Vector6d orbit_state_vec_ECEF = Vector6d::Zero();
    orbit_state_vec_ECEF = ECI2ECEF(ini_GPSorbtime, init_orbstate);
	
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////// ORBIT INITIALIZATION /////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    SPACEENV::n_max = nMAX;
    
    SC_Orbit.Setup(SC_prms,envmodels_paths);
    SC_Orbit.Init(ini_GPSorbtime, attstate, init_orbstate);
    SC_Orbit.drag_on = T_model[2];
    SC_Orbit.srp_on = T_model[3];
    SC_Orbit.sunmoon_on = T_model[4];
    SC_Orbit.Drag_Model = Drag_Model;
    SC_Orbit.SRP_Model = SRP_Model;
    SC_Orbit.simdur = SIM_DURATION;
    SC_Orbit.ForceModelsSetup();
    
    cout << "Start\n" << endl;
			
    //////////////////////////// Vector for csv files /////////////////////////
    VectorNd<14> orbit_state_vec;
    
    orbit_state_vec(0) = 0.0;
    orbit_state_vec(1) = ini_GPSorbtime;
    //orbit_state_vec(2) = 0.0;
    orbit_state_vec.segment(2,6) = init_orbstate;
    orbit_state_vec.segment(8,6) = orbit_state_vec_ECEF;
    
    // Orbit state
    ofstream orbstate_file;
    orbstate_file.open(orbfile_name);
    
    orbstate_file << fixed << orbit_state_vec(0) << "," << orbit_state_vec(1) << "," << orbit_state_vec(2) << "," << orbit_state_vec(3) << "," << orbit_state_vec(4) << "," << orbit_state_vec(5) << "," << orbit_state_vec(6) << "," << orbit_state_vec(7) << "," << orbit_state_vec(8) << "," << orbit_state_vec(9) << "," << orbit_state_vec(10) << "," << orbit_state_vec(11) << "," << orbit_state_vec(12) << "," << orbit_state_vec(13) << endl;
    
    // Accelerations (spacecraft body-fixed frame)
    ofstream accelerations_file;
    accelerations_file.open(csv_accelerations_name);
    
    accelerations_file << "GPS Time [s],GravR [m/s²],GravT [m/s²],GravN [m/s²],SunMoonR [m/s²],SunMoonT [m/s²],SunMoonN [m/s²],SRP_R [m/s²],SRP_T [m/s²],SRP_N [m/s²],DragR [m/s²],DragT [m/s²],DragN [m/s²],AccR [m/s²],AccT [m/s²],AccN [m/s²],dvR [m/s],dvT [m/s],dvN [m/s]" << endl;
    
    VectorNd<15> Accenv = VectorNd<15>::Zero();
    VectorNd<15> AccenvRTN = VectorNd<15>::Zero();
    Vec3d Acc_act = Vec3d::Zero();
	
    // Attitude state
    ofstream attstate_file;
    attstate_file.open(attfile_name);
	
    Vector6d orbprop_state = Vector6d::Zero();
    
    const int output_rows = SIM_DURATION/SIM_STEP;
    int step = 0;
    
    MatrixXd orbstate_to_file = MatrixXd::Zero(output_rows+1,14);
    orbstate_to_file.row(step) = orbit_state_vec;
    step++;
    MatrixXd accelerations_to_file = MatrixXd::Zero(output_rows+1,19);
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////// ORBIT PROPAGATION LOOP ////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    man_ind = 0;
    
    //double part_dur = SIM_DURATION/10.0;
    //int sim_done = 10;
    //bool barinit = false;
    int barwidth = 100;
    //int barpos = 0;
    
    for( double t = 0.0 ; t < SIM_DURATION ; t += SIM_STEP )
        {
        //forstart = chrono::high_resolution_clock::now();
        
        GPStime = ini_GPSorbtime + t + SIM_STEP;
        //////////////// Commanded orbit maneuvers /////////////////
        dv1_imp_CMD << 0, 0, 0;
        dv1_cont_CMD << 0, 0, 0;
        dv2_imp_CMD << 0, 0, 0;
        dv2_cont_CMD << 0, 0, 0;
        dv_CTRL << 0, 0, 0;
        
        ////////////////////////////////////////////// ORBIT MANEUVERS //////////////////////////////////////////////
        if( !all_enabled_orbman.empty() )
            {
            // Maneuver execution message    
            if( t >= all_enabled_orbman[man_ind].init_time )
                {
                cout << all_enabled_orbman[man_ind].name << ", dv = " << all_enabled_orbman[man_ind].ManVec(0) <<  " , " <<  all_enabled_orbman[man_ind].ManVec(1) <<  " , " << all_enabled_orbman[man_ind].ManVec(2) << " [m/s]\n" << endl;
                
                man_ind++;
                }
        
            // Impulsive maneuvers of propulsion system 1
            if( !impulsive_maneuvers1.empty() && impman_ind1 != impulsive_maneuvers1.end() )
              {
              m_init_time = impman_ind1->init_time;
              
              if( t >= m_init_time )
                {
                dv1_imp_CMD = impman_ind1->ManVec;
                man_name = impman_ind1->name;
                man_sys = man_name.substr(man_name.find("_") + 1);
                
                OrbitPropulsion1.dv_CMD = dv1_imp_CMD;
                OrbitPropulsion1.man_sys = man_sys;
                
                //cout << "Orbit propulsion system 1 impulsive maneuver in " << man_sys << "frame: dv = " << dv1_imp_CMD(0) <<  " , " <<  dv1_imp_CMD(1) <<  " , " << dv1_imp_CMD(2) << " [m/s]\n" << endl;
                
                dv1_imp_CMD = OrbitPropulsion1.Output(GPStime, attstate, orbprop_state);
                
                impman_ind1++;
                }
              }
            // Continuous maneuvers of propulsion system 1  
            if( !continuous_maneuvers1.empty() && contman_ind1 != continuous_maneuvers1.end() )
              {
              m_init_time = contman_ind1->init_time;
              
              if( t >= m_init_time )
                {
                dv1_cont_CMD = contman_ind1->ManVec;
                man_name = contman_ind1->name;
                man_sys = man_name.substr(man_name.find("_") + 1);
                
                OrbitPropulsion1.dv_CMD = dv1_cont_CMD;
                OrbitPropulsion1.man_sys = man_sys;
                
                //cout << "Orbit propulsion system 1 continuous maneuver in " << man_sys << "frame: dv = " << dv1_cont_CMD(0) <<  " , " <<  dv1_cont_CMD(1) <<  " , " << dv1_cont_CMD(2) << " [m/s]\n" << endl;
                
                dv1_cont_CMD = OrbitPropulsion1.Output(GPStime, attstate, orbprop_state);
                
                contman_ind1++;
                }
              }
              
            // Impulsive maneuvers of propulsion system 2
            if( !impulsive_maneuvers2.empty() && impman_ind2 != impulsive_maneuvers2.end() )
              {
              m_init_time = impman_ind2->init_time;
              
              if( t >= m_init_time )
                {
                dv2_imp_CMD = impman_ind2->ManVec;
                man_name = impman_ind2->name;
                man_sys = man_name.substr(man_name.find("_") + 1);
                
                OrbitPropulsion2.dv_CMD = dv2_imp_CMD;
                OrbitPropulsion2.man_sys = man_sys;
                
                //cout << "Orbit propulsion system 2 impulsive maneuver in " << man_sys << "frame: dv = " << dv2_imp_CMD(0) <<  " , " <<  dv2_imp_CMD(1) <<  " , " << dv2_imp_CMD(2) << " [m/s]\n" << endl;
                
                dv2_imp_CMD = OrbitPropulsion2.Output(GPStime, attstate, orbprop_state);
                
                impman_ind2++;
                }
              }
            // Continuous maneuvers of propulsion system 2  
            if( !continuous_maneuvers2.empty() && contman_ind2 != continuous_maneuvers2.end() )
              {
              m_init_time = contman_ind2->init_time;
              
              if( t >= m_init_time )
                {
                dv2_cont_CMD = contman_ind2->ManVec;
                man_name = contman_ind2->name;
                man_sys = man_name.substr(man_name.find("_") + 2);
                
                OrbitPropulsion2.dv_CMD = dv2_cont_CMD;
                OrbitPropulsion2.man_sys = man_sys;
                
                //cout << "Orbit propulsion system 2 continuous maneuver in " << man_sys << "frame: dv = " << dv2_cont_CMD(0) <<  " , " <<  dv2_cont_CMD(1) <<  " , " << dv2_cont_CMD(2) << " [m/s]\n" << endl;
                
                dv2_cont_CMD = OrbitPropulsion2.Output(GPStime, attstate, orbprop_state);
                
                contman_ind2++;
                }
              }
          }
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////  
        /////////////////////////////////// PUT YOUR ORBIT CONTROL FUNCTION HERE ///////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////  
        
        ////////////////////////////////////////// UNCOMMENT THE FOLLOWING /////////////////////////////////////////
        
        //double control_vec[4]; // control_vec[0] = GPSsecs time of maneuver issued by orbit controller control_vec[1], control_vec[1], control_vec[1] = components in RTN frame of dv vector issued by orbit controller
        //
        //if(orbctrl_on)
        //    {
        //    if( OrbCtrlType.compare("ControllerName") == 0 )
        //        {
        //        ////////////// User own code for orbit control //////////////
        //        control_vec = YourOrbitControlFunction();
        //        ////////////////////////////////////////////////////////////
        //        
        //        // Put maneuver in ctrl_maneuver vector for execution by propulsion system
        //        if( control_vec[0] != 0.0 && ( fabs(control_vec[1]) + fabs(control_vec[2]) + fabs(control_vec[3]) ) != 0.0 && GPStime >= control_vec[0] && ctrl_maneuver.size() == 1 ) // Maneuver time is different than 0 and there is at least one component of the maneuver vector which different than 0 
        //            {
        //            ctrl_maneuver[0].name = "ContinuousManeuver" + CtrProp_name + "_RTN";
        //            ctrl_maneuver[0].maneuver_on = true;
        //            ctrl_maneuver[0].init_time = control_vec[0] - ini_GPSorbtime; // Conversion in simulation time
        //            for(int k = 0; k < 3; k ++) ctrl_maneuver[0].ManVec(k) = control_vec[k+1];
        //            
        //            if( CtrProp_name.compare("Prop1") == 0 ) OrbitPropulsion1.impman2contman(ctrl_maneuver, SIM_STEP);
        //            else if( CtrProp_name.compare("Prop2") == 0 ) OrbitPropulsion2.impman2contman(ctrl_maneuver, SIM_STEP);
        //            
        //            ctrl_contman_ind = ctrl_maneuver.begin();
        //            }
        //        }
        //    
        //    if( ( ctrl_maneuver[0].init_time != 0 ) && ( ctrl_maneuver[0].ManVec.norm() != 0 ) && ( ctrl_contman_ind != ctrl_maneuver.end() ) )
        //        {
        //        m_init_time = ctrl_contman_ind->init_time;
        //        
        //        if( t >= m_init_time )
        //            {
        //            dv_CTRL = ctrl_contman_ind->ManVec;
        //            man_name = ctrl_contman_ind->name;
        //            man_sys = "RTN";
        //            
        //            if( CtrProp_name.compare("Prop1") == 0 )
        //                {
        //                OrbitPropulsion1.dv_CMD = dv_CTRL;
        //                OrbitPropulsion1.man_sys = man_sys;
        //                dv_CTRL = OrbitPropulsion1.Output(GPStime, attstate, orbprop_state);
        //                }
        //            else if( CtrProp_name.compare("Prop2") == 0 )
        //                {
        //                OrbitPropulsion2.dv_CMD = dv_CTRL;
        //                OrbitPropulsion2.man_sys = man_sys;
        //                dv_CTRL = OrbitPropulsion2.Output(GPStime, attstate, orbprop_state);
        //                }
        //         
        //            ctrl_contman_ind++;
        //            
        //            if(ctrl_contman_ind == ctrl_maneuver.end())
        //                {
        //                // Resize and reinitialize vector ctrl_maneuver
        //                ctrl_maneuver.resize(1);
        //                ctrl_maneuver[0].ManVec = Vec3d::Zero();
        //                ctrl_maneuver[0].name = "";
        //                ctrl_maneuver[0].init_time = 0.0;
        //                ctrl_maneuver[0].duration = 0.0;
        //                ctrl_maneuver[0].maneuver_on = false;
        //                // Reset control maneuver counter
        //                ctrl_contman_ind = ctrl_maneuver.begin();
        //                }
        //            }
        //        }
        //    }
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////  
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////  
          
        dv_CMD = dv1_imp_CMD + dv1_cont_CMD + dv2_imp_CMD + dv2_cont_CMD + dv_CTRL;
    
		//////////////// Do propagation step /////////////////
		
        SC_Orbit.state = attstate;
        SC_Orbit.Maneuver(dv_CMD);
    
        SC_Orbit.Integrate(t,SIM_STEP);
        Accenv = SC_Orbit.Acceleration_env;
    
        orbprop_state = SC_Orbit.orbstate;
        
        #ifdef VERBOSE
        
            cout << "Simulation time: " << t << "\n" << endl;
            
            cout << orbprop_state(0) << "   " << orbprop_state(1) << "   " << orbprop_state(2) << "   " << orbprop_state(3) << "   " << orbprop_state(4) << "   " << orbprop_state(5) << endl;
              
        #endif
        
        ///////////////////////////////////////////////////////
        
        if( AttitudeType.compare("Ephemeris") == 0 )
          {    
          ephem_row = loaded_ephem.row(ind);
          attstate = ephem_row.segment(1,4);
          
          ind++;
          }
        
        T_ECI2RTN = ECI2RTN_Matrix(orbprop_state);
        
        //////////////// Display propagation progress /////////////////
        RunStatusBar(t, SIM_DURATION, barwidth);
		
        //////////////////////////// Orbit ephemeris //////////////////////////////
        
        orbit_state_vec(0) = t + SIM_STEP;
        orbit_state_vec(1) = GPStime;
        orbit_state_vec.segment(2,6) = orbprop_state;
        
        //////////////////////////// Accelerations /////////////////////////
        
        AccenvRTN.segment(0,3) = T_ECI2RTN*Accenv.segment(0,3);
        AccenvRTN.segment(3,3) = T_ECI2RTN*Accenv.segment(3,3);
        AccenvRTN.segment(6,3) = T_ECI2RTN*Accenv.segment(6,3);
        AccenvRTN.segment(9,3) = T_ECI2RTN*Accenv.segment(9,3);
        AccenvRTN.segment(12,3) = T_ECI2RTN*Accenv.segment(12,3);
        
        if(realtime) // Write output file inside propagation loop
            {
            orbit_state_vec_ECEF = ECI2ECEF(GPStime,orbprop_state);
            orbit_state_vec.segment(8,6) = orbit_state_vec_ECEF;
        
            // Orbit state
            orbstate_file << fixed << orbit_state_vec(0) << "," << orbit_state_vec(1) << "," << orbit_state_vec(2) << "," << orbit_state_vec(3) << "," << orbit_state_vec(4) << "," << orbit_state_vec(5) << "," << orbit_state_vec(6) << "," << orbit_state_vec(7) << "," << orbit_state_vec(8) << "," << orbit_state_vec(9) << "," << orbit_state_vec(10) << "," << orbit_state_vec(11) << "," << orbit_state_vec(12) << "," << orbit_state_vec(13) << endl;
            // Accelerations
            accelerations_file << setprecision(20) << ini_GPSorbtime + t << "," << AccenvRTN(0) << "," << AccenvRTN(1) << "," << AccenvRTN(2) << "," << AccenvRTN(3) << "," << AccenvRTN(4) << "," << AccenvRTN(5) << "," << AccenvRTN(6) << "," << AccenvRTN(7) << "," << AccenvRTN(8) << "," << AccenvRTN(9) << "," << AccenvRTN(10) << "," << AccenvRTN(11) << "," << AccenvRTN(12) << "," << AccenvRTN(13) << "," << AccenvRTN(14) << "," << Acc_act(0) << "," << Acc_act(1) << "," << Acc_act(2) << endl;
            }
        else
            {
            // Orbit state
            orbstate_to_file.row(step) = orbit_state_vec;
            // Accelerations
            accelerations_to_file(step,0) = ini_GPSorbtime + t;
            accelerations_to_file.block<1,15>(step,1) = AccenvRTN;
            accelerations_to_file.block<1,3>(step,16) = Acc_act;
            }
        
        /////////////////////// Attitude state output ///////////////////////
        
        if( AttitudeType.compare("RTN_fixed") == 0 ) T_RTN2Body = RotationMatrix321(phi, theta, psi);
        if( AttitudeType.compare("Fixed") == 0 )
          {
          T_RTN2ECI = T_ECI2RTN.transpose();    
          ECItoBody = Quaternion2RotationMatrix(attstate);    
          T_RTN2Body = ECItoBody*T_RTN2ECI;
          }
        
        // Write attitude to attitude output file in case option "Ephemeris" is not selected
        if( AttitudeType.compare("Ephemeris") == 1 )
          {
          Vec3d euler_ang = EulerAngles321(T_RTN2Body);
          
          Vector6d attstateRTN;
          attstateRTN << mod(euler_ang(0),PI2)*RAD2DEG, mod(euler_ang(1),PI2)*RAD2DEG, mod(euler_ang(2),PI2)*RAD2DEG, 0.0, 0.0, 0.0;
          
          attitudeRTN_state_vec(0) = GPStime;
          attitudeRTN_state_vec.segment(1,6) = attstateRTN;
          
          if(realtime)
            {
            attstate_file << fixed << attitudeRTN_state_vec(0) << "," << attitudeRTN_state_vec(1) << "," << attitudeRTN_state_vec(2) << "," << attitudeRTN_state_vec(3) << "," << attitudeRTN_state_vec(4) << "," << attitudeRTN_state_vec(5) << "," << attitudeRTN_state_vec(6) << endl;
            }
          else
            {
            // TO BE IMPLEMENTED
            }
          }
		//////////////// Get new attitude in case it is RTN-fixed /////////////////
		if( AttitudeType.compare("RTN_fixed") == 0 )
          {
          T_RTN2Body = RotationMatrix321(phi, theta, psi);
          ECItoBody = T_RTN2Body*T_ECI2RTN;
      
          attstate = RotationMatrix2Quaternion(ECItoBody);
          }
          
        ++step;
        }
    
    if(!realtime)
        {
        #pragma omp parallel sections
            {
            #pragma omp section
                {
                for(int i = 1; i < step; i++)
                    {
                    GPStime = orbstate_to_file(i,1);
                    orbprop_state = orbstate_to_file.row(i).segment(2,6);
                    
                    orbit_state_vec_ECEF = ECI2ECEF(GPStime,orbprop_state);
                    
                    orbstate_to_file.row(i).segment(8,6) = orbit_state_vec_ECEF;
                    // Orbit state
                    orbstate_file << fixed << orbstate_to_file(i,0) << "," << orbstate_to_file(i,1) << "," << orbstate_to_file(i,2) << "," << orbstate_to_file(i,3) << "," << orbstate_to_file(i,4) << "," << orbstate_to_file(i,5) << "," << orbstate_to_file(i,6) << "," << orbstate_to_file(i,7) << "," << orbstate_to_file(i,8) << "," << orbstate_to_file(i,9) << "," << orbstate_to_file(i,10) << "," << orbstate_to_file(i,11) << "," << orbstate_to_file(i,12) << "," << orbstate_to_file(i,13) << endl;
                    }
                }
            
            #pragma omp section
                {
                for(int i = 1; i < step; i++)
                    {
                    // Aceelerations
                    accelerations_file << accelerations_to_file(i,0) << "," << accelerations_to_file(i,1) << "," << accelerations_to_file(i,2) << "," << accelerations_to_file(i,3) << "," << accelerations_to_file(i,4) << "," << accelerations_to_file(i,5) << "," << accelerations_to_file(i,6) << "," << accelerations_to_file(i,7) << "," << accelerations_to_file(i,8) << "," << accelerations_to_file(i,9) << "," << accelerations_to_file(i,10) << "," << accelerations_to_file(i,11) << "," << accelerations_to_file(i,12) << "," << accelerations_to_file(i,13) << "," << accelerations_to_file(i,14) << "," << accelerations_to_file(i,15) << "," << accelerations_to_file(i,16) << "," << accelerations_to_file(i,17) << "," << accelerations_to_file(i,18) << endl;
                    }
                }
            }
        }
    
    orbstate_file.close();
    
    accelerations_file.close();
	
	attstate_file.close();
    
    clockend = chrono::high_resolution_clock::now();
        
    chrono::duration<double,milli> elapsed_millisecs = clockend - clockstart;
    cout << "Elapsed seconds: " << elapsed_millisecs.count()/1000.0 << endl;
    
  return(0);  
  
  
  } // End of main()
  
  