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

#ifndef TRANSFORMATIONS_H_
#define TRANSFORMATIONS_H_

#include <string>
#include <Eigen/Core>
#include <Eigen/Geometry> 
#include <VarTypes.h>

extern "C"
      {
      #include <SpiceUsr.h>
      }
      
using namespace std;
using namespace math;

//------------------------------------------------------------------------------
// Group of functions for coordinate systems transformations
//------------------------------------------------------------------------------

// Compute rotation matrix from Euler angles (rotation a-b-c)
Mat3x3d RotationMatrix(double alpha, double beta, double gamma, int a, int b, int c);
// Compute rotation matrix from Euler angles (rotation 3-2-1)
Mat3x3d RotationMatrix321(double phi, double theta, double psi);
// Compute Euler angles (rotation 3-2-1) from rotation matrix
Vec3d EulerAngles321(Mat3x3d& Mat);
// Compute Euler angles (rotation 3-2-1) from quaternion
Vec3d EulerAngles321(Vec4d& q);
// Compute quaternion from rotation matrix
Vec4d RotationMatrix2Quaternion(Mat3x3d& Mat);
// Compute rotation matrix from quaternion
Mat3x3d Quaternion2RotationMatrix(Vec4d& q);
// Rotation by means of quaternion
Vec3d TransbyQ(Vec3d& v3D, Vec4d& q);
// Inversion of quaternion
Vec4d q_inv(Vec4d& q);
// Convert from rectangular coordinates to geodetic coordinates
Vec3d ECEF2lonlath(Vec3d& posECEF);
// Convert from rectangular coordinates to geodetic coordinates
Vec3d lonlath2ECEF(Vec3d& lonlath);
// Compute the elevation with respect to a reference point on Earth of a point given its geodetic coordinates
double lonlath2El(Vec3d& lonlath, double ref_lon, double ref_lat);
// Convert from topocentric horizon coordinates (SEZ) to ECEF coordinates
Vec3d SEZ2ECEF(Vec3d& posSEZ, double lon, double lat);
// Convert from ECEF coordinates to topocentric horizon coordinates (SEZ)
Vec3d ECEF2SEZ(Vec3d& posECEF, double lon, double lat);
// Compute the azimuth, elevation and altitude with respect to a reference point on Earth of a point vector in ECEF coordinates
Vec3d ECEF2AzElAlt(Vector6d& stateECEF, double ref_lon, double ref_lat);
// State transformation matrix from ECI to RTN frame
Mat3x3d ECI2RTN_Matrix(Vector6d& ECIstate);
// State transformation from ECI to RTN frame
Vector6d ECI2RTN(Vector6d& ECIstate);
// Vector transformation from RTN to ECI frame
Vec3d RTN2ECI(Vec3d& RTNvec, Vector6d& ECIstate);
// State transformation from ECEF (ITRF93) to ECI (J2000) frame
Vector6d ECEF2ECI(double GPStime, Vector6d& ECEFstate);
// State transformation from ECI (J2000) to ECEF (ITRF93) frame
Vector6d ECI2ECEF(double GPStime, Vector6d& ECIstate);
// Convert position vector from rectangular coordinates ECI to right ascension, declination and range ECI
Vec3d ECI2RAD(Vec3d& posECI);
// Transformation of a generic 3-dimensional vector between two frames
Vec3d v3D_transform(double GPStime, Vec3d& v3D, const string frame_from, const string frame_to);
// Gives unit vector of axis x1, y1 or z1-direction, expressed in x2, y2 and z2 coordinates given transformation matrix T1toT2 and the desired axis name x, y or z
Vec3d u_axis(Mat3x3d& T1toT2, const string axis_name);
// Conversion of Cartesian state vector to osculating Keplerian elements for non-circular non-equatorial orbits
Vector6d rv2Kepler(Vector6d& ECIstate, bool& valid);
// Conversion of Cartesian state vector to osculating orbital elements for non-circular non-equatorial orbits
Vector6d rv2oe(Vector6d& ECIstate, bool& valid);
// Conversion of osculating orbital elements to mean orbital elements
Vector6d osc2mean(Vector6d& osc_elem, bool& valid);
// Computation of eccentric anomaly for elliptic orbits
double EccAnomaly(double M, double e, const int maxit, const double eps);
// Compute coordinates of a point on Earth in spacecraft's ground-track coordinates given longitude and latitudes of the spacecraft and the point
VectorNd<2> lonlat2satcart(double lonS, double latS, double lonT, double latT, double inc, double ki);
// Conversion of GPS seconds to ephemeris time (TDT)
double GPS2ET(double GPSsecs);
// Conversion of GPS seconds to local solar time given the geodetic longitude
Vec3d GPS2LST(double GPSsecs, double lon);
// Conversion of GPS seconds to UTC string yyyy-mm-ddThh:mm:ss.sss
string GPS2UTCstr(double GPSsecs, const string format);
// Conversion of GPS seconds to UTC yyyy mm dd hh mm ss.sss
Vector6d GPS2UTCdate(double GPSsecs, const string format);
// Conversion of UTC yyyy mm dd hh mm ss.sss to GPS seconds
double UTCdate2GPSsecs(Vector6d& UTCdate);
// Conversion of GPS seconds to day of year
double GPS2DOY(double GPSsecs);
// Gives leap seconds at epoch GPS seconds
double leapsec(double GPSsecs);
// Conversion of GPS seconds to Unix time
double GPS2Unix(double GPSsecs);
// Conversion of GPS seconds to Julian date
double GPS2JD(double GPSsecs);
// Conversion of GPS seconds to modified Julian date
double GPS2MJD(double GPSsecs);
// Modulo function
double mod(double x, double y);
// Correct a vector into another vector which is on the surface of a cone whose axis is the initial vector and with aperture alpha 
Vec3d vectcorr_cone(Vec3d& v_in, double alpha, double theta);
// Correct a vector into another vector which has same direction and pointing but different norm
Vec3d vectcorr_norm(Vec3d& v_in, double newnorm);
// Compute quaternion from Euler angles
//Vec4d toQ(double phi, double theta, double psi);

#endif // PROPAGATION_H_