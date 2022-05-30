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

#ifndef VARIABLES_TYPES_H
#define VARIABLES_TYPES_H

#include <string>
#include <map>

#include <Eigen/Core>
#include "boost/multi_array.hpp"

#define GRAV_MAX_SIZE 121
#define FSW_GRAV_MAX_SIZE 31

//------------------------------------------------------------------------------
/**
* Definition of specialized variable types
*/
//------------------------------------------------------------------------------
namespace math
        {
        // Typedef for of N double
        template<unsigned int N>
        using VectorNd = Eigen::Matrix< double, N, 1 >;
        // Typedef for of N unsigned int
        template<unsigned int N>
        using VectorNui = Eigen::Matrix< unsigned int, N, 1 >;
        // Typedef for of N unsigned int
        template<unsigned int N>
        using VectorNi = Eigen::Matrix< int, N, 1 >;
        // Typedef for 3x3 matrix of doubles
        typedef Eigen::Matrix< double, 3, 3 > Mat3x3d;
        // Typedef for 6x6 matrix of doubles
        typedef Eigen::Matrix< double, 6, 6 > Mat6x6d;
        // Typedef for 121x121 matrix of doubles (gravity field model)
        typedef Eigen::Matrix< double, GRAV_MAX_SIZE, GRAV_MAX_SIZE > MatnMAXxnMAXd;
        // Typedef for 121x121 matrix of doubles (gravity field model)
        typedef Eigen::Matrix< double, GRAV_MAX_SIZE, 1 > VectornMAXd;
        // Typedef for 31x31 matrix of doubles (hardcoded gravity field model)
        typedef Eigen::Matrix< double, FSW_GRAV_MAX_SIZE, FSW_GRAV_MAX_SIZE > FSW_MatnMAXxnMAXd;
        // Typedef for 3D vector of doubles
        typedef Eigen::Vector3d Vec3d;
        // Typedef for 4D vector of doubles
        typedef Eigen::Vector4d Vec4d;
        // Typedef for Vector6d
        typedef Eigen::Matrix< double, 6, 1 > Vector6d;
        // Typedef for Vector6i
        typedef Eigen::Matrix< int, 6, 1 > Vector6i;
        // Typedef for Vector6f
        typedef Eigen::Matrix< float, 6, 1 > Vector6f;
        // Typedef for Matrix6d
        typedef Eigen::Matrix< double, 6, 6 > Matrix6d;
        // Typedef for Matrix6i
        typedef Eigen::Matrix< int, 6, 6 > Matrix6i;
        // Typedef for Matrix6f
        typedef Eigen::Matrix< float, 6, 6 > Matrix6f;
        // Typedef for Matrix3D
        typedef boost::multi_array<double, 3> Matrix3D;
        // Typedef for Matrix3D_index
        typedef Matrix3D::index Matrix3D_index;
        };
        
namespace SC
        {
        using namespace std;
        using namespace math;
        // Struct defining the parameters of one spacecraft's face
        struct Face
                  {
                  Vec3d n;         // Surface's unit normal vector in SC-body-fixed coordinates
                  Vec3d cP;        // Position of the surface's center of solar radiation pressure in SC-body-fixed coordinates
                  Vec3d cA;        // Position of the surface's center of aerodynamic pressure in SC-body-fixed coordinates
                  string Material; // Material composing the surface
                  double Area;     // Surface area
                  };
        // Struct containing spacecraft physical parameters
        struct SC_params
                       {
                        map<string, Face> Segment; // Map defining the surface of the spacecraft associated with a certain axis //Vec3d Faces;   // Vector containing the areas of the spacecraft's faces (simmetric spacecraft)  [m^2]
                        Mat3x3d MoI;               // Moments of inertia matrix
                        Vec3d Mdip;                // Spacecraft magnetic dipole moment vector
                        double SC_mass;            // Spacecraft mass
                        double CD;                 // Drag coefficient
                        double C_SRP;              // Radiation pressure coefficient
						double Area_D;             // Drag area to be used with atmospheric drag simple model
                        double Area_R;             // Radiation area to be used with solar radiation pressure simple model
                       };
        // Struct containing space environment models files paths
        struct EnvModels
                        {
                        string datapath;     // Planets ephemerides
                        string sunmoon;     // Planets ephemerides
                        string magneticfield;   // Magnetic field
                        string atmosphere;      // Atmospheric model
                        string gravityfield;         // Gravity field model
                        };             
        // Struct containing subsystem physical parameters
        struct SYS_params
                        {
                        EnvModels SpaceEnv;
                        Eigen::VectorXd ConstPrm;    // Constant parameters
						Eigen::VectorXd AuxPrm;     // Accuracy (e.g. estimation, actuation, etc.)
						Eigen::VectorXd Accuracy;     // Accuracy (e.g. estimation, actuation, etc.)
						Eigen::VectorXd OPS_limits;  // Operability limits (e.g. field of view limitations, etc.)
						Mat3x3d SC2SYS;       // Transformation matrix from SC body-fixed frame subsystem frame
                        Vec3d Position;       // Position of subsystem geometric center with respect to the spacecraft body-fixed frame center
                        string Name;
                        double Range;         // Operability range
                        int MaxUpdateRate;    // Rate of the output [Hz]
                        bool on_off;
                        };
                       
        struct maneuver
                      {
                      Vec3d ManVec;
                      string name;
                      double init_time;
                      double duration;
                      bool maneuver_on;
                      };
        };

namespace ground
    {
	// Target        
	struct TG
		{
		std::string name;
		double lon; // Longitude
		double lat;	// Latitude
		double alt;	// Altitude
		};
	// Target        
	struct GS
		{
		std::string name;
		double lon;	// Longitude
		double lat;	// Latitude
		double alt;	// Altitude
		double minelev;	// Minimum spacecraft elevation for acquisition of signal
		};
	};
	
namespace events
	{
	struct Pass
		{
		std::string Location_name, Epoch_in, Epoch_out; // GS name, start epoch UTC, end epoch UTC
		int GPSsecs_in; // Start epoch in GPS seconds
		double duration; // Pass duration
        double elev_in, elev_out; // Star elevation, end elevation 
        std::string maxel_time; // max elevation epoch UTC
        double maxel; // max elevation
		double Az_in, Az_out, Az_maxel;	// Azimuth at AOS, azimuth at LOS, azimuth at max elevation
		double lon, lat;	// Longitude and latitude of pass location
		double PP, SS; // Number of orbital plane, number of spacecraft
        std::string pass_type;
		
		bool operator<(Pass const& other_Pass) { return GPSsecs_in < other_Pass.GPSsecs_in; }
		};
        
//    struct TG_contact
//		{
//		std::string TG, Epoch_in, LOS; // GS name, AOS epoch UTC, LOS epoch UTC
//		int GPSsecs; // AOS GPS seconds
//		double dur; // Contact duration
//        std::string maxel_time; // max elevation epoch UTC
//        double maxel; // max elevation
//		double AOS_Az, LOS_Az, maxel_Az;	// Azimuth at AOS, azimuth at LOS, azimuth at max elevation
//		double lon, lat;	// Longitude and latitude of ground station
//		double PP, SS; // Number of orbital plane, number of spacecraft
//		
//		bool operator<(GS_contact const& other_GS_contact) { return GPSsecs < other_GS_contact.GPSsecs; }
//		};
	};
        
using namespace math;
// Sensor readings TC for hardware-in-the-loop simulations
struct sensorTCs
        {
        unsigned int UnixTime;
        
        VectorNui<6> CssRaw;
        
        VectorNi<2> Cam1Raw;
        unsigned int Cam1Busy;
        unsigned int Cam1Result;
        
        VectorNi<2> Cam2Raw;
        unsigned int Cam2Busy;
        unsigned int Cam2Result;
        
        VectorNi<3> MagRaw;
        VectorNi<3> RateRaw;
        VectorNi<3> WheelRaw;
        
        VectorNi<3> Star1Camera;
        VectorNi<3> Star1Inertial;
        
        VectorNi<3> Star2Camera;
        VectorNi<3> Star2Inertial;
        
        VectorNi<3> Star3Camera;
        VectorNi<3> Star3Inertial;
        };

#endif
