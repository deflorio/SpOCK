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

#include <Attitude.h>
#include <AttitudeQ.h>
#include <IO_utils.h>
#include <VarTypes.h>
#include <Constants.h>
#include <Transformations.h>
#include <HIL_interface.h>

#include <EarthSunSensor.h>
#include <Magneto.h>
#include <WheelMR.h>
#include <SolarPanel.h>
#include <Solarsys.h>

// External libraries: Eigen
#include <Eigen/Core>
#include "boost/multi_array.hpp"

#ifdef USE_SPICE

extern "C"
        {
        #include "extlib/cspice/include/SpiceUsr.h"
        }
      
#endif

using namespace std;
using namespace SC;
using namespace constants;
using namespace mathconst;
using namespace math;

using namespace earthsun;
using namespace magneto;
using namespace mrwheel;
using namespace solarpan;

using namespace boost;


int main(int argc, char *argv[])
    {
    chrono::time_point<chrono::high_resolution_clock> clockstart, clockend;
    
    clockstart = chrono::high_resolution_clock::now();
    
    /////////////// Simmulation parameters XML file ///////////////
	
	if(argc < 2)
      {
	  // Tell the user how to run the program
	  cerr << "Usage: " << argv[0] << " path_to/simulation/parameters/file.xml\nEx: ./bin/AttitudePropagator /home/username/path1/path2/input/simparam.xml" << endl;
	  return 1;
	  }
	  
	string XML_simparam_file(argv[1]);
    
    /////////////////////////////////////////////////////////////
    //////////// CREATE ATTITUDE PROPAGATOR CLASS ///////////////
    /////////////////////////////////////////////////////////////
    
    #ifdef USE_QUATERNION

        using namespace attitudeq;
        ATTQ SC_Attitude;
        //PROP* SC_Attitude_ptr = &SC_Attitude;
    
    #else
    
        using namespace attitude;
        ATT SC_Attitude;
    
    #endif
    
    ///////////////////////////////////////////////////////////////
    //////////// SIMULATION PARAMETERS VARIABLES //////////////////
    ///////////////////////////////////////////////////////////////
    
    ///////////////////////// Files paths /////////////////////////
    
    // Input files path
    string Orbit_ephemeris, Attitude_ephemeris, Data_path, planetephemeris, eop, pck_data, leapsecond, magneticfield, gravityfield, atmosphere;
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
    
    XML_parser(XML_simparam_file, Orbit_ephemeris, Attitude_ephemeris, Data_path, planetephemeris, eop, pck_data, leapsecond, magneticfield, gravityfield, atmosphere, orbfile_name, attfile_name, sensors_filename, csv_torques_name, csv_accelerations_name, SIM_STEP, SIM_DURATION, init_orbtime, init_orbstate, phi, theta, psi, om_x, om_y, om_z, initstate_in_RTN, realtime, realtime_wait, ggrad_on, mag_on, drag_on, srp_on, nMAX, sunmoon_on, Drag_Model, SRP_Model, AttitudeType, attctrl_on, AttCtrlType, orbctrl_on, OrbCtrlType, SC_mass, MoI, CoG, SC_Cd, SC_Cr, SC_Area_D, SC_Area_R, Mdip, F_Xplus, F_Xminus, F_Yplus, F_Yminus, F_Zplus, F_Zminus, Sensor_prm_SUN, Sensor_prm_EARTH, Sensor_prm_CSS1, Sensor_prm_CSS2, Sensor_prm_CSS3, Sensor_prm_CSS4, Sensor_prm_CSS5, Sensor_prm_CSS6, Sensor_prm_MAG, Sensor_prm_MAGstowed, Sensor_prm_RS, Sensor_prm_MAGTRQ, Sensor_prm_WHEEL1, Sensor_prm_WHEEL2, Sensor_prm_WHEEL3, Solarpan1_prm, Solarpan2_prm, Solarpan3_prm, OrbitPropulsion1_prm, OrbitPropulsion2_prm, all_maneuvers);
    //cout << "Sono qui" << endl;
    size_t lastslash = XML_simparam_file.find_last_of("/");
    string ReadXML_TXT_file_name = XML_simparam_file.substr(lastslash+1);
    size_t lastspoint = ReadXML_TXT_file_name.find_last_of(".");
    ReadXML_TXT_file_name = ReadXML_TXT_file_name.substr(0,lastspoint);
    const string ReadXML_TXT_file = XML_simparam_file.substr(0,lastslash) + "/Read_" + ReadXML_TXT_file_name + ".txt"; 
    //const string ReadXML_TXT_file = "input/readXML.txt";
    // Put read XML in a text file (for check purposes)
    ReadXMLtoTXT(ReadXML_TXT_file, Orbit_ephemeris, Attitude_ephemeris, Data_path, planetephemeris, eop, pck_data, leapsecond, magneticfield, gravityfield, atmosphere, orbfile_name, attfile_name, sensors_filename, csv_torques_name, csv_accelerations_name, SIM_STEP, SIM_DURATION, init_orbtime, init_orbstate, phi, theta, psi, om_x, om_y, om_z, initstate_in_RTN, realtime, realtime_wait, ggrad_on, mag_on, drag_on, srp_on, nMAX, sunmoon_on, Drag_Model, SRP_Model, AttitudeType, attctrl_on, AttCtrlType, orbctrl_on, OrbCtrlType, SC_mass, MoI, CoG, SC_Cd, SC_Cr, SC_Area_D, SC_Area_R, Mdip, F_Xplus, F_Xminus, F_Yplus, F_Yminus, F_Zplus, F_Zminus, Sensor_prm_SUN, Sensor_prm_EARTH, Sensor_prm_CSS1, Sensor_prm_CSS2, Sensor_prm_CSS3, Sensor_prm_CSS4, Sensor_prm_CSS5, Sensor_prm_CSS6, Sensor_prm_MAG, Sensor_prm_MAGstowed, Sensor_prm_RS, Sensor_prm_MAGTRQ, Sensor_prm_WHEEL1, Sensor_prm_WHEEL2, Sensor_prm_WHEEL3, Solarpan1_prm, Solarpan2_prm, Solarpan3_prm, OrbitPropulsion1_prm, OrbitPropulsion2_prm, all_maneuvers);

    // Load orbit ephemerides
    Eigen::MatrixXd loaded_ephem;
    try{ loaded_ephem = read_csvfile(Orbit_ephemeris.c_str(),15); }
    catch(const string errmsg)
          {
          cerr << "Orbit ephemerides: " + errmsg << endl;
          exit(EXIT_FAILURE);
          }
    
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
    envmodels_paths.planetephem = planetephemeris;
    envmodels_paths.magneticfield = magneticfield;
    envmodels_paths.gravityfield = gravityfield;
    envmodels_paths.atmosphere = atmosphere;
    
    ////////////////////////////// Outputs //////////////////////////////////
    
    // Files
    //string mat_attfile_name = "output/" + attfile_name + ".mat";            
    //string csv_attfile_name = "output/" + attfile_name + ".csv";
    
    ///////////////// ADCS sensors and actuators objects /////////////////
    
    EARTHSUNSENS SunCamera(Sensor_prm_SUN);
    
    EARTHSUNSENS EarthCamera(Sensor_prm_EARTH);
    
    EARTHSUNSENS CSS1(Sensor_prm_CSS1);
    CSS1.Init();
    
    EARTHSUNSENS CSS2(Sensor_prm_CSS2);
    CSS2.Init();
    
    EARTHSUNSENS CSS3(Sensor_prm_CSS3);
    CSS3.Init();
    
    EARTHSUNSENS CSS4(Sensor_prm_CSS4);
    CSS4.Init();
    
    EARTHSUNSENS CSS5(Sensor_prm_CSS5);
    CSS5.Init();
    
    EARTHSUNSENS CSS6(Sensor_prm_CSS6);
    CSS6.Init();
    
    Sensor_prm_MAG.SpaceEnv.magneticfield = envmodels_paths.magneticfield;
    MAGNETO Magnetometer1(Sensor_prm_MAG);
    Magnetometer1.Init();
    
    Sensor_prm_MAGstowed.SpaceEnv.magneticfield = envmodels_paths.magneticfield;
    MAGNETO Magnetometer2(Sensor_prm_MAGstowed);
    Magnetometer2.Init();
    
    Mat3x3d Mag_SC2SYS;
    if(Magnetometer1.On) Mag_SC2SYS = Sensor_prm_MAG.SC2SYS;
    if(Magnetometer2.On) Mag_SC2SYS = Sensor_prm_MAGstowed.SC2SYS;
    
    MRWHEEL RateSensor(Sensor_prm_RS);
    
    Sensor_prm_MAGTRQ.SpaceEnv.magneticfield = envmodels_paths.magneticfield;
    MAGNETO Magnetorquer(Sensor_prm_MAGTRQ);
    Magnetorquer.Init();
    
    MRWHEEL Wheel1(Sensor_prm_WHEEL1);
    Wheel1.Init();
    
    MRWHEEL Wheel2(Sensor_prm_WHEEL2);
    Wheel2.Init();
    
    MRWHEEL Wheel3(Sensor_prm_WHEEL3);
    Wheel3.Init();
    
    ///////////////// Solar panels objects /////////////////
    
    SOLARPAN SolarPanel1(Solarpan1_prm);
    SolarPanel1.Init();
    
    SOLARPAN SolarPanel2(Solarpan2_prm);
    SolarPanel2.Init();
    
    SOLARPAN SolarPanel3(Solarpan3_prm);
    SolarPanel3.Init();
    
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
    
    ///////////////////////// Initial state parameters //////////////////////
    
    // Orbit initial state vector
    double ini_GPStime, GPStime, omega_orb;
    Vec3d posECI = Vec3d::Zero();
    Vec3d velECI = Vec3d::Zero();
    Vec3d posECEF = Vec3d::Zero();
    Vector6d current_orbstate = Vector6d::Zero();
    Vector6d current_orbstateECEF = Vector6d::Zero();
    VectorNd<15> ephem_row = VectorNd<15>::Zero();
    
    ephem_row = loaded_ephem.row(0);
    ini_GPStime = ephem_row(1);
    
    current_orbstate = ephem_row.segment(2,6);
    posECI = current_orbstate.segment(0,3);
    velECI = current_orbstate.segment(3,3);
    
    current_orbstateECEF = ephem_row.segment(8,6);
    posECEF = current_orbstateECEF.segment(0,3);
    
    omega_orb = velECI.norm()/posECI.norm(); //sqrt( astro::GM_EARTH/(orbrho*orbrho*orbrho) );
    omega_orb = omega_orb*RAD2DEG;
    
    // Attitude initial state vector
    Vector6d init_attstateRTN;
    VectorXd init_attstate;
    
    // Conversion to radians and normalization to 2*pi
    phi = mod(phi*DEG2RAD,PI2);
    theta = mod(theta*DEG2RAD,PI2);
    psi = mod(psi*DEG2RAD,PI2);
    
    om_x = om_x*DEG2RAD;
    om_y = om_y*DEG2RAD;
    om_z = om_z*DEG2RAD;
    
    Mat3x3d ECItoBody, T_ECI2RTN, T_RTN2ECI, T_RTN2Body;
    // Quaternion state + angular velocity vector to be used by the sensor and actuator objects
    VectorNd<7> stateQ = VectorNd<7>::Zero();
    
    #ifdef USE_QUATERNION
    
        VectorNd<7> init_attstateQ = VectorNd<7>::Zero();
        
        if(initstate_in_RTN)
            {
            init_attstateRTN(0) = phi;
            init_attstateRTN(1) = theta;
            init_attstateRTN(2) = psi;
            init_attstateRTN(3) = om_x;
            init_attstateRTN(4) = om_y;
            init_attstateRTN(5) = om_z;
            
            cout << "\nInitial state (RTN): phi = " << phi*RAD2DEG << " theta = " << theta*RAD2DEG << " psi = " << psi*RAD2DEG << " deg," << " om_x = " << om_x*RAD2DEG << " om_y = " << om_y*RAD2DEG << " om_z = " << om_z*RAD2DEG << " deg/s\n" << endl;
        
            T_ECI2RTN = ECI2RTN_Matrix(current_orbstate);
            T_RTN2ECI = T_ECI2RTN.transpose();
            T_RTN2Body = RotationMatrix321(init_attstateRTN(0), init_attstateRTN(1), init_attstateRTN(2));
            
            ECItoBody = T_RTN2Body*T_ECI2RTN;
        
            init_attstateQ.segment(0,4) = RotationMatrix2Quaternion(ECItoBody);
            init_attstateQ.segment(4,3) = init_attstateRTN.segment(3,3);
            }
        else
            {
            ECItoBody = RotationMatrix321(phi, theta, psi);
            init_attstateQ.segment(0,4) = RotationMatrix2Quaternion(ECItoBody);
            init_attstateQ(4) = om_x;
            init_attstateQ(5) = om_y;
            init_attstateQ(6) = om_z;
            }
            
        init_attstate = init_attstateQ;
        stateQ = init_attstateQ;
            
        //cout << "Initial state: " << init_attstate(0) << "   " << init_attstate(1) << "   " << init_attstate(2) << "   " << init_attstate(3) << "   " << init_attstate(4)*RAD2DEG << "   " << init_attstate(5)*RAD2DEG << "   " << init_attstate(6)*RAD2DEG << endl;
        
    #else
    
        Vector6d init_attstate_ECI;
    
        if(initstate_in_RTN)
            {
            init_attstateRTN(0) = phi;
            init_attstateRTN(1) = theta;
            init_attstateRTN(2) = psi;
            init_attstateRTN(3) = om_x;
            init_attstateRTN(4) = om_y;
            init_attstateRTN(5) = om_z;
            
            cout << "\nInitial attitude state (RTN): phi = " << phi*RAD2DEG << " theta = " << theta*RAD2DEG << " psi = " << psi*RAD2DEG << " deg," << " om_x = " << om_x*RAD2DEG << " om_y = " << om_y*RAD2DEG << " om_z = " << om_z*RAD2DEG << " deg/s\n" << endl;
        
            T_ECI2RTN = ECI2RTN_Matrix(current_orbstate);
            T_RTN2ECI = T_ECI2RTN.transpose();
            T_RTN2Body = RotationMatrix321(init_attstateRTN(0), init_attstateRTN(1), init_attstateRTN(2));
            
            ECItoBody = T_RTN2Body*T_ECI2RTN;
            Vec3d EulerAng = EulerAngles321(ECItoBody);
            
            phi = EulerAng(0);
            theta = EulerAng(1);
            psi = EulerAng(2);
            
            init_attstate_ECI(0) = phi;
            init_attstate_ECI(1) = theta;
            init_attstate_ECI(2) = psi;
            init_attstate_ECI(3) = om_x;
            init_attstate_ECI(4) = om_y;
            init_attstate_ECI(5) = om_z;
            }
        else
            {
            ECItoBody = RotationMatrix321(phi, theta, psi);
        
            init_attstate_ECI(0) = phi;
            init_attstate_ECI(1) = theta;
            init_attstate_ECI(2) = psi;
            init_attstate_ECI(3) = om_x;
            init_attstate_ECI(4) = om_y;
            init_attstate_ECI(5) = om_z;
            }
        
        init_attstate = init_attstate_ECI;
        
        stateQ.segment(0,4) = RotationMatrix2Quaternion(ECItoBody);
        stateQ(4) = om_x;
        stateQ(5) = om_y;
        stateQ(6) = om_z;
    
    #endif
    
    //////////////////// Environment models //////////////////////////
    
    bool T_model[4] = {ggrad_on, mag_on, drag_on, srp_on};
    
    //////////////////// Run-start display message //////////////////////////
    
    string sim_duration_string = to_string(SIM_DURATION/86400.0) + " days";
    if(SIM_DURATION/86400.0 < 1.0) sim_duration_string = to_string(SIM_DURATION/3600.0) + " hours";
    if(SIM_DURATION/3600.0 < 1.0) sim_duration_string = to_string(SIM_DURATION/60.0) + " minutes";
    
    string pert_txt = "Simulation duration: " + sim_duration_string + "\n\nPERTURBATIONS\n\n";
    
    if(ggrad_on) pert_txt = pert_txt + "Gravity gradient\n";
    if(mag_on) pert_txt = pert_txt + "Earth's magnetic field: " + magneticfield + " model\n";
    if(drag_on) pert_txt = pert_txt + "Panels atmospheric drag model, Atmosphere: " + atmosphere + " model\n";
    if(srp_on) pert_txt = pert_txt + "Panels solar radiation pressure model\n";
    if( !(ggrad_on || mag_on || srp_on || drag_on) ) pert_txt = "NO PERTURBATIONS\n\n";
    
    cout << pert_txt << "\n" << endl;
    
    ///////////////// Commanded attitude maneuvers /////////////////
    
    // Sort maneuvers by type and put them in vectors
    vector<maneuver> wheels_maneuvers; // Maneuvers executed by reaction/momentum wheels
    vector<maneuver> magneto_maneuvers; // Maneuvers executed by magnetotorquers
    vector<maneuver> thrusters_maneuvers; // Maneuvers executed by thrusters
    
    //vector<maneuver>::iterator man_ind;
    unsigned int man_ind;
    vector<maneuver>::size_type man_size;
    
    string man_name;
    Vec3d ontime = Vec3d::Zero();
    Vec3d unit_ontime = Vec3d::Zero();
    double m_init_time; // Re-used in propagation loop
    maneuver unit_maneuver;
    
    for(man_ind = 0 ; man_ind < all_maneuvers.size(); man_ind++)
        {
        man_name = all_maneuvers[man_ind].name;
        
        if( man_name.compare("WheelsManeuver") == 0 ) wheels_maneuvers.push_back(all_maneuvers[man_ind]);
        else if( man_name.compare("MagnetoManeuver") == 0 ) magneto_maneuvers.push_back(all_maneuvers[man_ind]);
        else if( man_name.compare("ThrustersManeuver") == 0 ) thrusters_maneuvers.push_back(all_maneuvers[man_ind]);
        }
    
    man_size = magneto_maneuvers.size();
    // Process ontimes in magneto_maneuvers and thrusters_maneuvers. If an ontime is larger than  1 s it is splitted in slots of 1 s (e.g. if ontime = 3.25 s it is splitted in ontimes 1 + 1 + 1 + 0.25)
    if( !magneto_maneuvers.empty() )
      {
      for( man_ind = 0; man_ind < man_size; man_ind++ )
          {
          m_init_time = magneto_maneuvers[man_ind].init_time;
          ontime = magneto_maneuvers[man_ind].ManVec;
          ontime = ontime/1000.0; // Conversion of the ontime in seconds
          
          if( (ontime(0) > 1.0) || (ontime(1) > 1.0) || (ontime(2) > 1.0) )
            {
            unit_maneuver = magneto_maneuvers[man_ind];
            
            for(int i = 0 ; i < 3; i++)
                {
                if( ontime(i) > 1.0 )
                  {
                  ontime(i) = ontime(i) - 1.0;
                  m_init_time = m_init_time + 1.0;
                  unit_ontime(i) = 1.0;
                  }  
                }
                
            // Make place for the insertion of a new element
            man_size = man_size + 1;
            magneto_maneuvers.resize(man_size); 
            // Translate up the element over position man_ind
            for( unsigned int i = man_size-1; i > man_ind; i-- ) magneto_maneuvers[i] = magneto_maneuvers[i-1];
            // Change the element that was in position man_ind
            magneto_maneuvers[man_ind+1].init_time = m_init_time;
            magneto_maneuvers[man_ind+1].ManVec = ontime*1000.0; // Re-conversion in milliseconds
            // Insert new element in position man_ind
            unit_maneuver.ManVec = unit_ontime*1000.0;
            magneto_maneuvers[man_ind] = unit_maneuver; 
            }
          }
      }
    
    vector<maneuver>::iterator W_man_ind, M_man_ind, T_man_ind; // To be used in propagation loop
    W_man_ind = wheels_maneuvers.begin();
    M_man_ind = magneto_maneuvers.begin();
    T_man_ind = thrusters_maneuvers.begin();
    
    man_size = magneto_maneuvers.size();
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////// ATTITUDE INITIALIZATION /////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    VectorNd<12> orbstateECI_ECEF;
    orbstateECI_ECEF << current_orbstate, current_orbstateECEF;
        
    SC_Attitude.ggrad_on = T_model[0];
    SC_Attitude.mag_on = T_model[1];
    SC_Attitude.drag_on = T_model[2];
    SC_Attitude.srp_on = T_model[3];
    SC_Attitude.magnetometer_on = Magnetometer1.On || Magnetometer1.On; // True if at least one of the magnetometers is switched-on
    
    SC_Attitude.Setup(SC_prms,envmodels_paths);
    SC_Attitude.Init(ini_GPStime, init_attstate, orbstateECI_ECEF);
    SC_Attitude.ForceModelsSetup();
      
    double eps_abs = 1E-8;
    double eps_rel = 0.0;
    double factor_x = 0.0;
    double factor_dxdt = 0.0;
    
    SC_Attitude.StepperSetup(eps_abs, eps_rel, factor_x, factor_dxdt);
    
    ////////////////////////////////////////////////////////////////////
    ///////////////// ADCS SENSORS INTERFACE-VARIABLES /////////////////
    ////////////////////////////////////////////////////////////////////
        
    // Time
    unsigned int UnixTime = GPS2Unix(ini_GPStime);
    
    // CSS
    VectorNd<6> CssRaw_d;
    VectorNui<6> CssRaw;
    VectorNd<1> Css_out;
    
    Css_out = CSS1.Output(ini_GPStime, stateQ.segment(0,4), posECI);
    CssRaw_d(0) = Css_out(0);
    Css_out = CSS2.Output(ini_GPStime, stateQ.segment(0,4), posECI);
    CssRaw_d(1) = Css_out(0);
    Css_out = CSS3.Output(ini_GPStime, stateQ.segment(0,4), posECI);
    CssRaw_d(2) = Css_out(0);
    Css_out = CSS4.Output(ini_GPStime, stateQ.segment(0,4), posECI);
    CssRaw_d(3) = Css_out(0);
    Css_out = CSS5.Output(ini_GPStime, stateQ.segment(0,4), posECI);
    CssRaw_d(4) = Css_out(0);
    Css_out = CSS6.Output(ini_GPStime, stateQ.segment(0,4), posECI);
    CssRaw_d(5) = Css_out(0);
    for(int i = 0 ; i < 6; i++) CssRaw_d(i) = round(CssRaw_d(i)); // Before casting to unsigned integer we have to round because the casting to integer rounds always down
    CssRaw = CssRaw_d.cast<unsigned int>(); // Cast type of Eigen library variable
    
    // Sun camera
    Vec4d camera_out;
    
    VectorNd<2> SunRaw_d;
    VectorNi<2> SunRaw;
    
    camera_out = SunCamera.Output(ini_GPStime, stateQ.segment(0,4), posECI);
    for(int i = 0 ; i < 2; i++) SunRaw_d(i) = round(camera_out(i));
    SunRaw = SunRaw_d.cast<int>();
    
    unsigned int SunBusy = (unsigned int)camera_out(2);
    unsigned int SunResult = (unsigned int)camera_out(3);
    
    // Earth camera
    VectorNd<2> NadirRaw_d;
    VectorNi<2> NadirRaw;
    
    camera_out = EarthCamera.Output(ini_GPStime, stateQ.segment(0,4), posECI);
    for(int i = 0 ; i < 2; i++) NadirRaw_d(i) = round(camera_out(i));
    NadirRaw = NadirRaw_d.cast<int>();
    
    unsigned int NadirBusy = (unsigned int)camera_out(2);
    unsigned int NadirResult = (unsigned int)camera_out(3);
    
    // Magnetometer
    VectorNd<3> MagRaw_d = VectorNd<3>::Zero();
    VectorNi<3> MagRaw = VectorNi<3>::Zero(); // Magnetic field vector in magnetometer frame nanotesla/2.5
    VectorNd<3> MagRaw_RTN_d = VectorNd<3>::Zero();
    VectorNi<3> MagRaw_RTN = VectorNi<3>::Zero();
    
    if(Magnetometer1.On && Magnetometer2.On) cerr << "WARNING: both magnetometers are on. Only the second magnetometer listed in the XML input file will be considered active" << endl;
    if(Magnetometer1.On) MagRaw_d = Magnetometer1.Output(ini_GPStime, stateQ.segment(0,4), posECEF);
    if(Magnetometer2.On) MagRaw_d = Magnetometer2.Output(ini_GPStime, stateQ.segment(0,4), posECEF);
    for(int i = 0 ; i < 3; i++) MagRaw_d(i) = round(MagRaw_d(i)*1E-3); // Before casting to unsigned integer we have to round because the casting to integer rounds always down. 1E-3: conversion from nanotesla to microtesla
    MagRaw = MagRaw_d.cast<int>(); // Cast type of Eigen library variable
    
    // Rate sensor
    VectorNd<3> RateRaw_d;
    VectorNi<3> RateRaw;
    
    RateRaw_d = RateSensor.Output(ini_GPStime, stateQ, current_orbstate);
    for(int i = 0 ; i < 3; i++) RateRaw_d(i) = round(RateRaw_d(i)*1E3); // 1E3: conversion from deg to millideg
    RateRaw = RateRaw_d.cast<int>();
    
    // Star sensors
    VectorNi<3> Star1Camera = VectorNi<3>::Zero();
    VectorNi<3> Star1Inertial = VectorNi<3>::Zero();
    
    VectorNi<3> Star2Camera = VectorNi<3>::Zero();
    VectorNi<3> Star2Inertial = VectorNi<3>::Zero();
    
    VectorNi<3> Star3Camera = VectorNi<3>::Zero();
    VectorNi<3> Star3Inertial = VectorNi<3>::Zero();
    
    //////////////////////////// Vector for csv file /////////////////////////
    
    VectorNi<27> sensors_output;
    
    sensors_output(0) = ini_GPStime;
    
    sensors_output.segment(1,6) = CssRaw.cast<int>();
    
    sensors_output.segment(7,2) = SunRaw;
    sensors_output(9) = SunBusy;
    sensors_output(10) = SunResult;
    
    sensors_output.segment(11,2) = NadirRaw;
    sensors_output(13) = NadirBusy;
    sensors_output(14) = NadirResult;
    
    sensors_output.segment(15,3) = MagRaw;
    
    //sensors_output.segment(18,3) = VectorNi<3>::Zero()
    
    sensors_output.segment(21,3) = RateRaw;
    
    // WheelRaw: see below
    
    //sensors_output.segment(28,3) = Star1Camera;
    //sensors_output.segment(31,3) = Star1Inertial;
    //
    //sensors_output.segment(34,3) = Star2Camera;
    //sensors_output.segment(37,3) = Star2Inertial;
    //
    //sensors_output.segment(40,3) = Star3Camera;
    //sensors_output.segment(43,3) = Star3Inertial;
    
    ////////////////////////////////////////////////////////////////////
    ///////////////// ADCS ACTUATORS INTERFACE-VARIABLES ///////////////
    ////////////////////////////////////////////////////////////////////
    
    // Magnetotorquer ontimes commanded by ADCS controller
    Vec3d m_ontimes = Vec3d::Zero();
    // Magnetotorquer ontimes commanded manually
    Vec3d m_ontimesCMD = Vec3d::Zero();
    // Vector containing the speeds commanded manually of the 3 reaction/momentum wheels
    Vec3d wheelspeedCMD;
    // Vector containing the speeds commanded by ADCS controller of the 3 reaction/momentum wheels
    Vec3d wheelspeedCTRL;
    // Magnetotorquer torque
    Vec3d Tm = Vec3d::Zero();
    // Reaction/momentum wheels torque
    Vec3d Tw = Vec3d::Zero();
    // Total torque given by actuators
    Vec3d T_act = Vec3d::Zero();
    
    Vec3d T_CTRL = Vec3d::Zero();
    // Vector of reaction wheels angular momentums in spacecraft frame
    Vec3d hw_SC = Vec3d::Zero();
    // Environment perturbation torques
    VectorNd<15> TenvRTN = VectorNd<15>::Zero();
    VectorNd<15> Torque_env = VectorNd<15>::Zero();
    
    Wheel1.wheelspeedCMD = Wheel1.wheelspeed_SYS;
    Wheel2.wheelspeedCMD = Wheel2.wheelspeed_SYS;
    Wheel3.wheelspeedCMD = Wheel3.wheelspeed_SYS;
    
    Vec3d WheelRaw_d; // Velocities of momentum wheels [RPM]
    VectorNi<3> WheelRaw; // Velocities of momentum wheels [RPM]
    Vec3d Wheel_out;
    
    Wheel_out = Wheel1.Output(ini_GPStime, stateQ, current_orbstate);
    WheelRaw_d(0) = Wheel_out(2);
    Wheel_out = Wheel2.Output(ini_GPStime, stateQ, current_orbstate);
    WheelRaw_d(1) = Wheel_out(2);
    Wheel_out = Wheel3.Output(ini_GPStime, stateQ, current_orbstate);
    WheelRaw_d(2) = Wheel_out(2);
    //cout << "Sono qui\n" << endl;
    for(int i = 0 ; i < 3; i++) WheelRaw_d(i) = round(WheelRaw_d(i));
    WheelRaw = WheelRaw_d.cast<int>();
    
    sensors_output.segment(24,3) = WheelRaw;
    
    // Struct for all sensors and actuators (see VarTypes.h)
    sensorTCs SensorReadings;
    sensorTCs* sensorTCs_ptr = &SensorReadings;
    
    //////////////////////////////////////////////////////
    ///////////////// SOLAR PANELS OUTPUT ///////////////
    /////////////////////////////////////////////////////
    
    Vec3d panel1_out, panel2_out, panel3_out;
    
    panel1_out = SolarPanel1.Output(ini_GPStime, stateQ.segment(0,4), posECI);
    panel2_out = SolarPanel2.Output(ini_GPStime, stateQ.segment(0,4), posECI);
    panel3_out = SolarPanel3.Output(ini_GPStime, stateQ.segment(0,4), posECI);
    
    VectorNd<8> panels_output;
    
    panels_output(0) = ini_GPStime;
    
    panels_output.segment(1,2) = panel1_out.segment(0,2);
    panels_output.segment(3,2) = panel2_out.segment(0,2);
    panels_output.segment(5,2) = panel3_out.segment(0,2);
    panels_output(7) = panel1_out(2);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////// ATTITUDE INITIALIZATION /////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    cout << "Start\n" << endl;
    //chrono::time_point<chrono::high_resolution_clock> forstart, forend;
    
    //////////////////////////// CSV files /////////////////////////
    VectorNd<7> attitudeRTN_state_vec;
    Vector6d attstateRTN;
    
    // Attitude state RTN
    ofstream state_file;
    state_file.open(attfile_name);
    
    T_ECI2RTN = ECI2RTN_Matrix(current_orbstate);
    T_RTN2ECI = T_ECI2RTN.transpose();
    T_RTN2Body = ECItoBody*T_RTN2ECI;
    Vec3d euler_ang = EulerAngles321(T_RTN2Body);
    
    attstateRTN << mod(euler_ang(0),PI2)*RAD2DEG, mod(euler_ang(1),PI2)*RAD2DEG, mod(euler_ang(2),PI2)*RAD2DEG, om_x*RAD2DEG, om_y*RAD2DEG, om_z*RAD2DEG;
    
    attitudeRTN_state_vec(0) = ini_GPStime;
    attitudeRTN_state_vec.segment(1,6) = attstateRTN;
    
    state_file << fixed << attitudeRTN_state_vec(0) << "," << attitudeRTN_state_vec(1) << "," << attitudeRTN_state_vec(2) << "," << attitudeRTN_state_vec(3) << "," << attitudeRTN_state_vec(4) << "," << attitudeRTN_state_vec(5) << "," << attitudeRTN_state_vec(6) << endl;
    
    // Attitude state quaternions
    
    string attfile_name_quat = attfile_name.replace(attfile_name.end() - 4,attfile_name.end(),"_Quaternions.csv");
    ofstream state_file_quat;
    state_file_quat.open(attfile_name_quat);
    
    state_file_quat << fixed << ini_GPStime << "," <<stateQ(0) << "," << stateQ(1) << "," << stateQ(2) << "," << stateQ(3) << "," << stateQ(4) << "," << stateQ(5) << "," << stateQ(6) << endl;
    
    // Sensors output
    ofstream sensors_file;
    sensors_file.open(sensors_filename);
    
    sensors_file << "GPS Time [s],CSS1,CSS2,CSS3,CSS4,CSS5,CSS6,SunCam Az,SunCam El,SunCapture,SunDetection,EarthCam Az,EarthCam El,EarthCapture,EarthDetection,MagX_magframe,MagY_magframe,MagZ_magframe,MagR_RTN,MagT_RTN,MagN_RTN,RateX_sensframe,RateY_sensframe,RateZ_sensframe,WheelSpeedX,WheelSpeedY,WheelSpeedZ,SolarPanels1 Pow,SolarPanels1 Curr,SolarPanels2 Pow,SolarPanels2 Curr,SolarPanels3 Pow,SolarPanels3 Curr,Eclipse" << endl;
    
    Mat3x3d Mag_SYS2SC;
    Mag_SYS2SC = Mag_SC2SYS.transpose();
    MagRaw_RTN_d = (T_RTN2Body.transpose())*Mag_SYS2SC*MagRaw_d;
    MagRaw_RTN = MagRaw_RTN_d.cast<int>();
    
    sensors_output.segment(18,3) = MagRaw_RTN;
    
    sensors_file << fixed << sensors_output(0) << "," << sensors_output(1) << "," << sensors_output(2) << "," << sensors_output(3) << "," << sensors_output(4) << "," << sensors_output(5) << "," << sensors_output(6) << "," << sensors_output(7) << "," << sensors_output(8) << "," << sensors_output(9) << "," << sensors_output(10) << "," << sensors_output(11) << "," << sensors_output(12) << "," << sensors_output(13) << "," << sensors_output(14) << "," << sensors_output(15) << "," << sensors_output(16) << "," << sensors_output(17) << "," << sensors_output(18) << "," << sensors_output(19) << "," << sensors_output(20) << "," << sensors_output(21) << "," << sensors_output(22) << "," << sensors_output(23) << "," << sensors_output(24) << "," << sensors_output(25) << "," << sensors_output(26) << "," << panels_output(1) << "," << panels_output(2) << "," << panels_output(3) << "," << panels_output(4) << "," << panels_output(5) << "," << panels_output(6) << "," << panels_output(7) << endl;
    
    // Torques (spacecraft body-fixed frame)
    ofstream torques_file;
    torques_file.open(csv_torques_name);
    
    torques_file << "GPS Time [s],TenvX,TenvY,TenvZ,TenvR,TenvT,TenvN,TmX,TmY,TmZ,TwX,TwY,TwZ,TactX,TactY,TactZ,Mag. ontimes X,Mag. ontimes Y,Mag. ontimes Z,WheelCTRL X,WheelCTRL Y,WheelCTRL Z" << endl;
    
    // Environmental torques in RTN frame
    string csv_envtorques_name = csv_torques_name.replace(csv_torques_name.end() - 4,csv_torques_name.end(),"_EnvironmentRTN.csv");
    
    ofstream env_torques_file;
    
    env_torques_file.open(csv_envtorques_name);
    
    env_torques_file << "GPS Time [s],GravR,GravT,GravN,MagR,MagT,MagN,SRP_R,SRP_T,SRP_N,DragR,DragT,DragN,T_R,T_T,T_N" << endl;
    
    //////////////////////////// Matrices for csv files /////////////////////////
    
    const int output_rows = SIM_DURATION/SIM_STEP;
    int step = 0;
    
    MatrixXd attstate_to_file = MatrixXd::Zero(output_rows+1,11);
    attstate_to_file.row(step).segment(0,7) = attitudeRTN_state_vec;
    attstate_to_file.row(step).segment(7,4) = stateQ.segment(0,4);
    
    MatrixXd sensors_to_file = MatrixXd::Zero(output_rows+1,34);
    sensors_to_file.row(step).segment(0,27) = sensors_output.segment(0,27).cast<double>();
    sensors_to_file.row(step).segment(27,7) = panels_output.segment(1,7);
    
    MatrixXd torques_to_file = MatrixXd::Zero(output_rows+1,22);
    torques_to_file(step,0) = ini_GPStime;
    torques_to_file.row(step).segment(1,3) = Torque_env.segment(12,3);
    torques_to_file.row(step).segment(4,3) = TenvRTN.segment(12,3);
    torques_to_file.row(step).segment(7,3) = Tm;
    torques_to_file.row(step).segment(10,3) = Tw;
    torques_to_file.row(step).segment(13,3) = T_act;
    torques_to_file.row(step).segment(16,3) = m_ontimes;
    torques_to_file.row(step).segment(19,3) = wheelspeedCMD;
    
    MatrixXd envtorquesRTN_to_file = MatrixXd::Zero(output_rows+1,16);
    envtorquesRTN_to_file(step,0) = ini_GPStime;
    envtorquesRTN_to_file.row(step).segment(1,15) = TenvRTN;
    
    step++;
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////// ATTITUDE PROPAGATION LOOP ///////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    VectorXd prop_state;
    
    double part_dur = SIM_DURATION/10.0;
    double sim_done = 10.0;
            
    for( double t = 0.0 ; t < SIM_DURATION ; t += SIM_STEP )
        {
        //////////////// Do propagation step /////////////////
        m_ontimes = Vec3d::Zero();  
          
        //forstart = chrono::high_resolution_clock::now();
        
        ephem_row = loaded_ephem.row(step);
        
        current_orbstate = ephem_row.segment(2,6);
        posECI = current_orbstate.segment(0,3);
        velECI = current_orbstate.segment(3,3);
        //omega_orb = velECI.norm()/posECI.norm();
        
        current_orbstateECEF = ephem_row.segment(8,6);
        posECEF = current_orbstateECEF.segment(0,3);
        GPStime = ephem_row(1);
        
        //////////////// Do propagation step /////////////////
        
        orbstateECI_ECEF << current_orbstate, current_orbstateECEF;
        // Insert orbit state
        SC_Attitude.orbstate = orbstateECI_ECEF;
        // Insert reaction wheels' angular momentum
        SC_Attitude.hw = hw_SC;
        // Insert total torque generated by the actuators
        SC_Attitude.Maneuver(T_act);
        // Integration of the trajectory at time t
        SC_Attitude.Integrate(t,SIM_STEP);
        // Environmental torques
        Torque_env = SC_Attitude.Torque_env;
        // Current attitude state
        prop_state = SC_Attitude.state;
        
        #ifdef VERBOSE
    
            cout << "Step: " << t << "\n" << endl;
            
            #ifdef USE_QUATERNION
                cout << prop_state(0) << "   " << prop_state(1) << "   " << prop_state(2) << "   " << prop_state(3) << "   " << prop_state(4)*RAD2DEG << "   " << prop_state(5)*RAD2DEG << "   " << prop_state(6)*RAD2DEG << endl;
            #else
                cout << prop_state(0)*RAD2DEG << "   " << prop_state(1)*RAD2DEG << "   " << prop_state(2)*RAD2DEG << "   " << prop_state(3)*RAD2DEG << "   " << prop_state(4)*RAD2DEG << "   " << prop_state(5)*RAD2DEG << endl;
            #endif
          
        #endif
        
        //////////////// Display propagation progress /////////////////
        if( (t - part_dur) >= 0)
            {
            cout << (int)sim_done << "% done" << endl;
            part_dur = part_dur + SIM_DURATION/10.0;
            sim_done = sim_done + 10.0;
            }
        
        /////////////// ADCS sensors/actuators inputs/outputs ////////////////
        
        #ifdef USE_QUATERNION
        
            stateQ.segment(0,4) = prop_state.segment(0,4);
            stateQ.segment(4,3) = prop_state.segment(4,3);
        
        #else
        
            ECItoBody = RotationMatrix321(prop_state(0), prop_state(1), prop_state(2));
            stateQ.segment(0,4) = RotationMatrix2Quaternion(ECItoBody);
            stateQ.segment(4,3) = prop_state.segment(3,3);
            
            //cout << "prop_state: " << t << "   " << prop_state(0)*RAD2DEG << "   " << prop_state(1)*RAD2DEG << "   " << prop_state(2)*RAD2DEG << endl;
            //cout << "ECItoBody: " << t << "   " << ECItoBody << endl;
        
        #endif
        
        // Time
        UnixTime = GPS2Unix(GPStime);
        
        // CSS
        //VectorNd<10> CssRaw_d;
        Css_out = CSS1.Output(GPStime, stateQ.segment(0,4), posECI);
        CssRaw_d(0) = Css_out(0);
        Css_out = CSS2.Output(GPStime, stateQ.segment(0,4), posECI);
        CssRaw_d(1) = Css_out(0);
        Css_out = CSS3.Output(GPStime, stateQ.segment(0,4), posECI);
        CssRaw_d(2) = Css_out(0);
        Css_out = CSS4.Output(GPStime, stateQ.segment(0,4), posECI);
        CssRaw_d(3) = Css_out(0);
        Css_out = CSS5.Output(GPStime, stateQ.segment(0,4), posECI);
        CssRaw_d(4) = Css_out(0);
        Css_out = CSS6.Output(GPStime, stateQ.segment(0,4), posECI);
        CssRaw_d(5) = Css_out(0);
        for(int i = 0 ; i < 6; i++) CssRaw_d(i) = round(CssRaw_d(i)); // Before casting to unsigned integer we have to round because the casting to integer rounds always down
        CssRaw = CssRaw_d.cast<unsigned int>(); // Cast type of Eigen library variable
        
        //CssRaw = ( CSS.Output(GPStime, stateQ.segment(0,4), posECI) ).cast<unsigned int>(); // Cast type of Eigen library variable
        
        // Sun camera
        camera_out = SunCamera.Output(GPStime, stateQ.segment(0,4), posECI);
        for(int i = 0 ; i < 2; i++) SunRaw_d(i) = round(camera_out(i));
        SunRaw = SunRaw_d.cast<int>();
        
        //camera_out = SunCamera.Output(GPStime, stateQ.segment(0,4), posECI);
        //SunRaw = ( camera_out.segment(0,2) ).cast<int>();
        
        SunBusy = (unsigned int)camera_out(2);
        SunResult = (unsigned int)camera_out(3);
        
        // Earth camera
        camera_out = EarthCamera.Output(GPStime, stateQ.segment(0,4), posECI);
        for(int i = 0 ; i < 2; i++) NadirRaw_d(i) = round(camera_out(i));
        NadirRaw = NadirRaw_d.cast<int>();
        
        //camera_out = EarthCamera.Output(GPStime, stateQ.segment(0,4), posECI);
        //NadirRaw = ( camera_out.segment(0,2) ).cast<int>();
        
        NadirBusy = (unsigned int)camera_out(2);
        NadirResult = (unsigned int)camera_out(3);
        
        // Magnetometer
        
        if(Magnetometer1.On) MagRaw_d = Magnetometer1.Output(GPStime, stateQ.segment(0,4), posECEF);
        if(Magnetometer2.On) MagRaw_d = Magnetometer2.Output(GPStime, stateQ.segment(0,4), posECEF);
        for(int i = 0 ; i < 3; i++) MagRaw_d(i) = round(MagRaw_d(i)*1E-3); // Before casting to unsigned integer we have to round because the casting to integer rounds always down. 1E-3: conversion from nanotesla to microtesla
        MagRaw = MagRaw_d.cast<int>(); // Cast type of Eigen library variable
        
        // MagRaw = ( Magnetometer1.Output(GPStime, stateQ.segment(0,4), posECEF) ).cast<int>();
        
        // Rate sensor
        RateRaw_d = RateSensor.Output(GPStime, stateQ, current_orbstate);
        for(int i = 0 ; i < 3; i++) RateRaw_d(i) = round(RateRaw_d(i)*1E3); // 1E3: conversion from deg to millideg
        RateRaw = RateRaw_d.cast<int>();
        
        // RateRaw = ( RateSensor.Output(GPStime, stateQ, current_orbstate) ).cast<int>();
    
        // Star sensors
        
        // NOT IMPLEMENTED YET
        
        // Reaction/momentum wheels speed
        
        Wheel_out = Wheel1.Output(GPStime, stateQ, current_orbstate);
        WheelRaw_d(0) = Wheel_out(2);
        Wheel_out = Wheel2.Output(GPStime, stateQ, current_orbstate);
        WheelRaw_d(1) = Wheel_out(2);
        Wheel_out = Wheel3.Output(GPStime, stateQ, current_orbstate);
        WheelRaw_d(2) = Wheel_out(2);
        
        for(int i = 0 ; i < 3; i++) WheelRaw_d(i) = round(WheelRaw_d(i));
        WheelRaw = WheelRaw_d.cast<int>();
        
        sensors_output.segment(24,3) = WheelRaw;
        
        // WheelRaw = ( Wheels.Output(GPStime, stateQ, current_orbstate) ).cast<int>();
        
        ///////////////// IO with hardware /////////////////
        
        if(realtime)
            {
            SensorReadings.UnixTime = UnixTime;
            
            SensorReadings.CssRaw = CssRaw;
            
            SensorReadings.Cam1Raw = SunRaw;
            SensorReadings.Cam1Busy = SunBusy;
            SensorReadings.Cam1Result = SunResult;
            
            SensorReadings.Cam2Raw = NadirRaw;
            SensorReadings.Cam2Busy = NadirBusy;
            SensorReadings.Cam2Result = NadirResult;
            
            SensorReadings.MagRaw = MagRaw;
            SensorReadings.RateRaw = RateRaw;
            SensorReadings.WheelRaw = WheelRaw;
            
            SensorReadings.Star1Camera = Star1Camera;
            SensorReadings.Star1Inertial = Star1Inertial;
            
            SensorReadings.Star2Camera = Star2Camera;
            SensorReadings.Star2Inertial = Star2Inertial;
            
            SensorReadings.Star3Camera = Star3Camera;
            SensorReadings.Star3Inertial = Star3Inertial;
          
            send_receiveTCs(m_ontimes, wheelspeedCTRL, sensorTCs_ptr);
            m_ontimes = m_ontimes*10.0; // Unit of measure of magnetorquer commanded on time is 10ms units
            }
        
        ///////////////// Commanded attitude maneuvers /////////////////
        
        // Magnetotorquers
        if( !magneto_maneuvers.empty() && M_man_ind != magneto_maneuvers.end() ) // If there are 
          {
          m_init_time = M_man_ind->init_time;
          
          if( M_man_ind->maneuver_on && t >= m_init_time )
            {
            m_ontimesCMD = M_man_ind->ManVec;
            m_ontimes = m_ontimes + m_ontimesCMD;
            
            M_man_ind++;
            }
          }
        
        Magnetorquer.ontime = m_ontimes;
        Tm = Magnetorquer.Output(GPStime, stateQ.segment(0,4), posECEF);
        
        // Reaction/momentum wheels
        wheelspeedCMD = wheelspeedCTRL;
        
        if( !wheels_maneuvers.empty() && W_man_ind != wheels_maneuvers.end() ) // If there are 
          {
          m_init_time = W_man_ind->init_time;
          
          if( W_man_ind->maneuver_on && t >= m_init_time )
            {
            wheelspeedCMD = W_man_ind->ManVec;
            wheelspeedCMD = wheelspeedCMD + wheelspeedCTRL;
            W_man_ind++;
            }
          }
        
        // Wheels angular momentum vector
        hw_SC = Wheel1.hw_SC + Wheel2.hw_SC + Wheel3.hw_SC;
        
        // Give to wheels model commanded speeds
        Wheel1.wheelspeedCMD(2) = wheelspeedCMD(0);
        Wheel2.wheelspeedCMD(2) = wheelspeedCMD(1);
        Wheel3.wheelspeedCMD(2) = wheelspeedCMD(2);
        // Wheels torque
        Tw = Wheel1.Tw_SC + Wheel2.Tw_SC + Wheel3.Tw_SC;
        //cout << "wheelspeedCMD: " << wheelspeedCMD << "      WheelRaw: " << WheelRaw << endl;
        
        /////////////////////////////////// PUT YOUR ATTITUDE CONTROL FUNCTION HERE ///////////////////////////////////
        
        //if(attctrl_on)
        //    {
        //    T_CTRL = My_attitude_control_function(prop_state, AttCtrlType, input2, input3);    
        //    }
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        // Actuators' torques in S/C body-fixed frame
        T_act = Tm + Tw + T_CTRL;
        
        //////////////////////////// Vector for mat file /////////////////////////
        //attitude_state(step) = t+1; attitude_state(arr_len+step) = loaded_ephem(step,1); attitude_state(2*arr_len+step) = loaded_ephem(step,2); attitude_state(3*arr_len+step) = mod(prop_state(0),PI2); attitude_state(4*arr_len+step) = mod(prop_state(1),PI2); attitude_state(5*arr_len+step) = mod(prop_state(2),PI2); attitude_state(6*arr_len+step) = prop_state(3); attitude_state(7*arr_len+step) = prop_state(4); attitude_state(8*arr_len+step) = prop_state(5);
        
        //////////////////////////// Vector for csv file /////////////////////////
    
        sensors_output(0) = GPStime;
    
        sensors_output.segment(1,6) = CssRaw.cast<int>();
        
        sensors_output.segment(7,2) = SunRaw;
        sensors_output(9) = SunBusy;
        sensors_output(10) = SunResult;
        
        sensors_output.segment(11,2) = NadirRaw;
        sensors_output(13) = NadirBusy;
        sensors_output(14) = NadirResult;
        
        sensors_output.segment(15,3) = MagRaw;
        
        //sensors_output.segment(18,3) = VectorNi<3>::Zero()
        
        sensors_output.segment(21,3) = RateRaw;
        
        /////////////////////// Attitude state output ///////////////////////
        
        #ifdef USE_QUATERNION
        
            Vec4d q_attstate = prop_state.segment(0,4);
            ECItoBody = Quaternion2RotationMatrix(q_attstate);
        
        #else
        
            //Vec3d EulerAng = prop_state.segment(0,3);
            //
            //phi = EulerAng(0);
            //theta = EulerAng(1);
            //psi = EulerAng(2);
            //
            //ECItoBody = RotationMatrix321(phi, theta, psi);
        
        #endif
        
        T_ECI2RTN = ECI2RTN_Matrix(current_orbstate);
        T_RTN2ECI = T_ECI2RTN.transpose();
        
        T_RTN2Body = ECItoBody*T_RTN2ECI;
        
        Vec3d euler_ang_RTN = EulerAngles321(T_RTN2Body);
        
        //Vector6d attstateRTN;
        attstateRTN << mod(euler_ang_RTN(0),PI2)*RAD2DEG, mod(euler_ang_RTN(1),PI2)*RAD2DEG, mod(euler_ang_RTN(2),PI2)*RAD2DEG, stateQ(4)*RAD2DEG, stateQ(5)*RAD2DEG, stateQ(6)*RAD2DEG;
        
        //cout << "\nPitch - Roll - Yaw: " << attstateRTN(0) << "   " << attstateRTN(1) << "   " << attstateRTN(2) << "   om_x - om_y - om_z (body frame): " << attstateRTN(3) << "   " << attstateRTN(4) << "   " << attstateRTN(5) << endl;
        
        //////////////////// Matrix for csv file ////////////////////
        //attstate_to_file(step,0) = GPStime; attstate_to_file.block<1,6>(step,1) = attstateRTN;
        
        //////////////////////////// Vector for csv file /////////////////////////
        attitudeRTN_state_vec(0) = GPStime;
        attitudeRTN_state_vec.segment(1,6) = attstateRTN;
        
        Mat3x3d T_Body2RTN = T_RTN2Body.transpose();
        
        TenvRTN.segment(0,3) = T_Body2RTN*Torque_env.segment(0,3);
        TenvRTN.segment(3,3) = T_Body2RTN*Torque_env.segment(3,3);
        TenvRTN.segment(6,3) = T_Body2RTN*Torque_env.segment(6,3);
        TenvRTN.segment(9,3) = T_Body2RTN*Torque_env.segment(9,3);
        TenvRTN.segment(12,3) = T_Body2RTN*Torque_env.segment(12,3);
        
        //////////////////////////// Solar panels output /////////////////////////
        
        panel1_out = SolarPanel1.Output(GPStime, stateQ.segment(0,4), posECI);
        panel2_out = SolarPanel2.Output(GPStime, stateQ.segment(0,4), posECI);
        panel3_out = SolarPanel3.Output(GPStime, stateQ.segment(0,4), posECI);
        
        panels_output(0) = GPStime;
        
        panels_output.segment(1,2) = panel1_out.segment(0,2);
        panels_output.segment(3,2) = panel2_out.segment(0,2);
        panels_output.segment(5,2) = panel3_out.segment(0,2);
        panels_output(7) = panel1_out(2);
        
        //////////////////////////// Write csv files //////////////////////////////
        
        Mag_SYS2SC = Mag_SC2SYS.transpose();
        MagRaw_RTN_d = (T_RTN2Body.transpose())*Mag_SYS2SC*MagRaw_d;
        MagRaw_RTN = MagRaw_RTN_d.cast<int>();
        
        sensors_output.segment(18,3) = MagRaw_RTN;
        
        if(realtime) // Write output file inside propagation loop
            {
            // Attitude state
            state_file << fixed << attitudeRTN_state_vec(0) << "," << attitudeRTN_state_vec(1) << "," << attitudeRTN_state_vec(2) << "," << attitudeRTN_state_vec(3) << "," << attitudeRTN_state_vec(4) << "," << attitudeRTN_state_vec(5) << "," << attitudeRTN_state_vec(6) << endl;
            
            state_file_quat << fixed << GPStime << "," << stateQ(0) << "," << stateQ(1) << "," << stateQ(2) << "," << stateQ(3) << "," << stateQ(4) << "," << stateQ(5) << "," << stateQ(6) << endl;
            
            sensors_file << fixed << sensors_output(0) << "," << sensors_output(1) << "," << sensors_output(2) << "," << sensors_output(3) << "," << sensors_output(4) << "," << sensors_output(5) << "," << sensors_output(6) << "," << sensors_output(7) << "," << sensors_output(8) << "," << sensors_output(9) << "," << sensors_output(10) << "," << sensors_output(11) << "," << sensors_output(12) << "," << sensors_output(13) << "," << sensors_output(14) << "," << sensors_output(15) << "," << sensors_output(16) << "," << sensors_output(17) << "," << sensors_output(18) << "," << sensors_output(19) << "," << sensors_output(20) << "," << sensors_output(21) << "," << sensors_output(22) << "," << sensors_output(23) << "," << sensors_output(24) << "," << sensors_output(25) << "," << sensors_output(26) << "," << panels_output(1) << "," << panels_output(2) << "," << panels_output(3) << "," << panels_output(4) << "," << panels_output(5) << "," << panels_output(6) << "," << panels_output(7) << endl;
            
            // Torques (spacecraft body-fixed frame)
            //torques_file << setprecision(15) << Torques(step,0) << "," << Torques(step,1) << "," << Torques(step,2) << "," << Torques(step,3) << "," << Torques(step,4) << "," << Torques(step,5) << "," << Torques(step,6) << "," << Torques(step,7) << "," << Torques(step,8) << "," << Torques(step,9) << "," << Torques(step,10) << "," << Torques(step,11) << "," << Torques(step,12) << "," << Torques(step,13) << "," << Torques(step,14) << "," << Torques(step,15) << "," << m_ontimes(0) << "," << m_ontimes(1) << "," << m_ontimes(2) << "," << wheelspeedCTRL(0) << "," << wheelspeedCTRL(1) << "," << wheelspeedCTRL(2) << endl;
            
            torques_file << setprecision(15) << GPStime << "," << Torque_env(12) << "," << Torque_env(13) << "," << Torque_env(14) << "," << TenvRTN(12) << "," << TenvRTN(13) << "," << TenvRTN(14) << "," << Tm(0) << "," << Tm(1) << "," << Tm(2) << "," << Tw(0) << "," << Tw(1) << "," << Tw(2) << "," << T_act(0) << "," << T_act(1) << "," << T_act(2) << "," << m_ontimes(0) << "," << m_ontimes(1) << "," << m_ontimes(2) << "," << wheelspeedCMD(0) << "," << wheelspeedCMD(1) << "," << wheelspeedCMD(2) << endl;
            
            env_torques_file << setprecision(15) << GPStime << "," << TenvRTN(0) << "," << TenvRTN(1) << "," << TenvRTN(2) << "," << TenvRTN(3) << "," << TenvRTN(4) << "," << TenvRTN(5) << "," << TenvRTN(6) << "," << TenvRTN(7) << "," << TenvRTN(8) << "," << TenvRTN(9) << "," << TenvRTN(10) << "," << TenvRTN(11) << "," << TenvRTN(12) << "," << TenvRTN(13) << "," << TenvRTN(14) << endl;
            
            int dur_mill = realtime_wait*1E3; //SIM_STEP*1E3 - (int)elapsed_millisecs.count();
            this_thread::sleep_for( chrono::milliseconds(dur_mill) );
            }
        else // Put output in matrices
            {
            attstate_to_file.row(step).segment(0,7) = attitudeRTN_state_vec;
            attstate_to_file.row(step).segment(7,4) = stateQ.segment(0,4);
                
            sensors_to_file.row(step).segment(0,27) = sensors_output.segment(0,27).cast<double>();
            sensors_to_file.row(step).segment(27,7) = panels_output.segment(1,7);
            
            torques_to_file(step,0) = GPStime;
            torques_to_file.row(step).segment(1,3) = Torque_env.segment(12,3);
            torques_to_file.row(step).segment(4,3) = TenvRTN.segment(12,3);
            torques_to_file.row(step).segment(7,3) = Tm;
            torques_to_file.row(step).segment(10,3) = Tw;
            torques_to_file.row(step).segment(13,3) = T_act;
            torques_to_file.row(step).segment(16,3) = m_ontimes;
            torques_to_file.row(step).segment(19,3) = wheelspeedCMD;
            
            envtorquesRTN_to_file(step,0) = GPStime;
            envtorquesRTN_to_file.row(step).segment(1,15) = TenvRTN;
            }
        
        step++;
        }
    
    ////////////////////////////////////////////////////////
    ////////////////// WRITE OUTPUT FILES //////////////////
    ////////////////////////////////////////////////////////
    
    //////////////////////////// Write csv files //////////////////////////////
    
    if(!realtime)
        {
    #pragma omp parallel sections
            {
            #pragma omp section
                {
                for(int i = 1 ; i < step; i++)
                    {
                    state_file << fixed << attstate_to_file(i,0) << "," << attstate_to_file(i,1) << "," << attstate_to_file(i,2) << "," << attstate_to_file(i,3) << "," << attstate_to_file(i,4) << "," << attstate_to_file(i,5) << "," << attstate_to_file(i,6) << endl;
                    
                    state_file_quat << fixed << attstate_to_file(i,0) << "," << attstate_to_file(i,7) << "," << attstate_to_file(i,8) << "," << attstate_to_file(i,9) << "," << attstate_to_file(i,10) << "," << attstate_to_file(i,4) << "," << attstate_to_file(i,5) << "," << attstate_to_file(i,6) << endl;
                    }
                }
            #pragma omp section
                {
                for(int i = 1 ; i < step; i++)
                    {
                    sensors_file << fixed << sensors_to_file(i,0) << "," << sensors_to_file(i,1) << "," << sensors_to_file(i,2) << "," << sensors_to_file(i,3) << "," << sensors_to_file(i,4) << "," << sensors_to_file(i,5) << "," << sensors_to_file(i,6) << "," << sensors_to_file(i,7) << "," << sensors_to_file(i,8) << "," << sensors_to_file(i,9) << "," << sensors_to_file(i,10) << "," << sensors_to_file(i,11) << "," << sensors_to_file(i,12) << "," << sensors_to_file(i,13) << "," << sensors_to_file(i,14) << "," << sensors_to_file(i,15) << "," << sensors_to_file(i,16) << "," << sensors_to_file(i,17) << "," << sensors_to_file(i,18) << "," << sensors_to_file(i,19) << "," << sensors_to_file(i,20) << "," << sensors_to_file(i,21) << "," << sensors_to_file(i,22) << "," << sensors_to_file(i,23) << "," << sensors_to_file(i,24) << "," << sensors_to_file(i,25) << "," << sensors_to_file(i,26) << "," << sensors_to_file(i,27) << "," << sensors_to_file(i,28) << "," << sensors_to_file(i,29) << "," << sensors_to_file(i,30) << "," << sensors_to_file(i,31) << "," << sensors_to_file(i,32) << "," << sensors_to_file(i,33) << "," << endl;
                    }
                }
            #pragma omp section
                {
                for(int i = 1 ; i < step; i++)
                    {
                    torques_file << setprecision(15) << torques_to_file(i,0) << "," << torques_to_file(i,1) << "," << torques_to_file(i,2) << "," << torques_to_file(i,3) << "," << torques_to_file(i,4) << "," << torques_to_file(i,5) << "," << torques_to_file(i,6) << "," << torques_to_file(i,7) << "," << torques_to_file(i,8) << "," << torques_to_file(i,9) << "," << torques_to_file(i,10) << "," << torques_to_file(i,11) << "," << torques_to_file(i,12) << "," << torques_to_file(i,13) << "," << torques_to_file(i,14) << "," << torques_to_file(i,15) << "," << torques_to_file(i,16) << "," << torques_to_file(i,17) << "," << torques_to_file(i,18) << "," << torques_to_file(i,19) << "," << torques_to_file(i,20) << "," << torques_to_file(i,21) << endl;
                    }
                }
            #pragma omp section
                {
                for(int i = 1 ; i < step; i++)
                    {
                    env_torques_file << setprecision(15) << envtorquesRTN_to_file(i,0) << "," << envtorquesRTN_to_file(i,1) << "," << envtorquesRTN_to_file(i,2) << "," << envtorquesRTN_to_file(i,3) << "," << envtorquesRTN_to_file(i,4) << "," << envtorquesRTN_to_file(i,5) << "," << envtorquesRTN_to_file(i,6) << "," << envtorquesRTN_to_file(i,7) << "," << envtorquesRTN_to_file(i,8) << "," << envtorquesRTN_to_file(i,9) << "," << envtorquesRTN_to_file(i,10) << "," << envtorquesRTN_to_file(i,11) << "," << envtorquesRTN_to_file(i,12) << "," << envtorquesRTN_to_file(i,13) << "," << envtorquesRTN_to_file(i,14) << "," << envtorquesRTN_to_file(i,15) << endl;
                    }
                }
            }
                    
        state_file.close();
        state_file_quat.close();
        sensors_file.close();
        torques_file.close();
        env_torques_file.close();
        
        clockend = chrono::high_resolution_clock::now();
        
        chrono::duration<double,milli> elapsed_millisecs = clockend - clockstart;
        //cout << "Elapsed_millisecs: " << elapsed_millisecs.count() << endl;
        cout << "Elapsed seconds: " << elapsed_millisecs.count()/1000.0 << endl;
        }
    
  return(0);
  
  } // End of main()
  
