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
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

//#include <iomanip>

#include <Transformations.h>
#include <Constants.h>
#include <Eigen/Core>
#include <Eigen/Geometry> 
#include <VarTypes.h>
#include <NutCoeff.h>

//#include <boost/math/quaternion.hpp>

#ifdef USE_SPICE
extern "C"
      {
      #include <SpiceUsr.h> //#include "extlib/cspice/include/SpiceUsr.h"
      }
#endif

using namespace std;
using namespace math;
using namespace Eigen;
using namespace constants;
using namespace mathconst;
using namespace astro;



//------------------------------------------------------------------------------
// Mat3x3d Rot_x(double alpha)
//------------------------------------------------------------------------------
/**
 * Compute elementary x-rotation matrix
 *
 * @param alpha   Rotation around x-axis
 *
 * @return Rotation matrix
 */
//------------------------------------------------------------------------------  
Mat3x3d Rot_x(double alpha)
            {
            Mat3x3d Tx;
            double calpha = cos(alpha);
            double salpha = sin(alpha);
            
            Tx(0,0) = 1.0;     Tx(0,1) = 0.0;         Tx(0,2) = 0.0;
                        
            Tx(1,0) = 0.0;     Tx(1,1) = +calpha;     Tx(1,2) = +salpha;
            
            Tx(2,0) = 0.0;     Tx(2,1) = -salpha;     Tx(2,2) = +calpha;
      
            return(Tx);
            };
//------------------------------------------------------------------------------
// Mat3x3d Rot_y(double alpha)
//------------------------------------------------------------------------------
/**
 * Compute elementary y-rotation matrix
 *
 * @param alpha   Rotation around y-axis
 *
 * @return Rotation matrix
 */
//------------------------------------------------------------------------------  
Mat3x3d Rot_y(double alpha)
            {
            Mat3x3d Ty;
            double calpha = cos(alpha);
            double salpha = sin(alpha);
            
            Ty(0,0) = +calpha;    Ty(0,1) = 0.0;    Ty(0,2) = -salpha;
                        
            Ty(1,0) = 0.0;        Ty(1,1) = 1.0;    Ty(1,2) = 0.0;
            
            Ty(2,0) = +salpha;    Ty(2,1) = 0.0;    Ty(2,2) = +calpha;
      
            return(Ty);
            };
//------------------------------------------------------------------------------
// Mat3x3d Rot_z(double alpha)
//------------------------------------------------------------------------------
/**
 * Compute elementary z-rotation matrix
 *
 * @param alpha   Rotation around z-axis
 *
 * @return Rotation matrix
 */
//------------------------------------------------------------------------------  
Mat3x3d Rot_z(double alpha)
            {
            Mat3x3d Tz;
            double calpha = cos(alpha);
            double salpha = sin(alpha);
            
            Tz(0,0) = +calpha;    Tz(0,1) = +salpha;    Tz(0,2) = 0.0;
                        
            Tz(1,0) = -salpha;    Tz(1,1) = +calpha;    Tz(1,2) = 0.0;
            
            Tz(2,0) = 0.0;        Tz(2,1) = 0.0;        Tz(2,2) = 1.0;
      
            return(Tz);
            };
//-------------------------------------------------------------------------------------
// Mat3x3d RotationMatrix(double alpha, double beta, double gamma, int a, int b, int c)
//-------------------------------------------------------------------------------------
/**
 * Compute rotation matrix from Euler angles (rotation a-b-c) 
 *
 * @param phi, theta, psi   Input Euler angles
 * @param a, b, c           Integers number representing Euler rotation convention
 *
 * @return 3x3 rotation matrix
 */
//-------------------------------------------------------------------------------------             
Mat3x3d RotationMatrix(double alpha,
                       double beta,
                       double gamma,
                       int a,
                       int b,
                       int c)
                    {
                    Mat3x3d T, T1, T2, T3;
                    
                    switch(a)
                         {
                          case 1: T1 = AngleAxisd(alpha, Vector3d::UnitX()).inverse();
                          case 2: T1 = AngleAxisd(alpha, Vector3d::UnitY()).inverse();
                          case 3: T1 = AngleAxisd(alpha, Vector3d::UnitZ()).inverse();
                         }
                         
                    switch(b)
                         {
                          case 1: T2 = AngleAxisd(beta, Vector3d::UnitX()).inverse();
                          case 2: T2 = AngleAxisd(beta, Vector3d::UnitY()).inverse();
                          case 3: T2 = AngleAxisd(beta, Vector3d::UnitZ()).inverse();
                         }
                         
                    switch(c)
                         {
                          case 1: T3 = AngleAxisd(gamma, Vector3d::UnitX()).inverse();
                          case 2: T3 = AngleAxisd(gamma, Vector3d::UnitY()).inverse();
                          case 3: T3 = AngleAxisd(gamma, Vector3d::UnitZ()).inverse();
                         }
                    
                    T = T3*T2*T1;
                        
                    return(T);    
                    };
//------------------------------------------------------------------------------
// RotationMatrix321(double phi, double theta, double psi)
//------------------------------------------------------------------------------
/**
 * Compute rotation matrix from Euler angles (rotation 3-2-1) 
 *
 * @param phi, theta, psi   Input Euler angles
 *
 * @return 3x3 rotation matrix
 */
//------------------------------------------------------------------------------                      
Mat3x3d RotationMatrix321(double phi,
                          double theta,
                          double psi)
                        {
                        Mat3x3d T;
                        
                        double c_phi, s_phi, c_theta, s_theta, c_psi, s_psi;
                        
                        c_phi = cos(phi);      s_phi = sin(phi);
                        c_theta = cos(theta);    s_theta = sin(theta);
                        c_psi = cos(psi);      s_psi = sin(psi);
                        
                        T(0,0) = c_theta*c_phi;                         T(0,1) = c_theta*s_phi;                         T(0,2) = -s_theta;
                        
                        T(1,0) = s_psi*s_theta*c_phi - c_psi*s_phi;     T(1,1) = s_psi*s_theta*s_phi + c_psi*c_phi;     T(1,2) = s_psi*c_theta;
                        
                        T(2,0) = c_psi*s_theta*c_phi + s_psi*s_phi;     T(2,1) = c_psi*s_theta*s_phi - s_psi*c_phi;     T(2,2) = c_psi*c_theta;
                        
                        return(T);
                        };
//------------------------------------------------------------------------------
// Vec3d EulerAngles321(Mat3x3d& Mat) (overloaded function)
//------------------------------------------------------------------------------
/**
 * Compute Euler angles (rotation 3-2-1) from rotation matrix
 *
 * @param Mat   Input rotation matrix
 *
 * @return Euler angles for Euler rotation 3-2-1
 */
//------------------------------------------------------------------------------  
Vec3d EulerAngles321(Mat3x3d& Mat)
                    {
                    double phi, theta, psi;
                    double c2;
                    
                    psi = atan2( Mat(1,2) , Mat(2,2) );
                    
                    c2 = sqrt( Mat(0,0)*Mat(0,0) + Mat(0,1)*Mat(0,1) );
                    
                    theta = atan2( -Mat(0,2) , c2 );
                    
                    phi = atan2( Mat(2,0)*sin(psi) - Mat(1,0)*cos(psi) , Mat(1,1)*cos(psi) - Mat(2,1)*sin(psi) );
                       
                    Vec3d angles(phi, theta, psi);
                    
                    return(angles);
                    };
//------------------------------------------------------------------------------
// Vec3d EulerAngles321(Vec4d& q) (overloaded function)
//------------------------------------------------------------------------------
/**
 * Compute Euler angles (rotation 3-2-1) from quaternion
 *
 * @param q   Input quaternion (q = q1*i + q2*j + q3*k + q4)
 *
 * @return Euler angles for Euler rotation 3-2-1
 */
//------------------------------------------------------------------------------  
Vec3d EulerAngles321(Vec4d& q)
                    {
                    double q1, q2, q3, q4, phi, theta, psi;
                    
                    for(int i = 0; i < 4; i++) if( fabs(q(i)) < 1.0E-6 ) q(i) = 0.0;
                
                    q1 = q(0);  q2 = q(1);  q3 = q(2);  q4 = q(3);
                    
                    phi = atan2( 2.0*(q1*q2 + q3*q4) , 1.0 - 2.0*(q2*q2 + q3*q3) );
                    
                    theta = asin(-2.0*(q1*q3 - q2*q4));
                    
                    psi = atan2( 2.0*(q2*q3 + q1*q4) , 1.0 - 2.0*(q1*q1 + q2*q2) );
                    
                    //cout << "tanphi = " << 2.0*(q2*q3 + q1*q4)/(1.0 - 2.0*(q1*q1 + q2*q2)) << endl;
                    //cout << "sintheta = " << -2.0*(q1*q3 - q2*q4) << endl;
                    //cout << "tanpsi = " << 2.0*(q1*q2 + q3*q4)/(1.0 - 2.0*(q2*q2 + q3*q3)) << endl;
                    
                    //if(exp2 >= 0.499)
                    //  { // singularity at north pole
                    //   phi = 0.0;
                    //   theta = mathconst::PI/2;
                    //   psi = 2.0*atan2(q1,q4);
                    //  }
                    //  
                    //if(exp2 <= -0.499)
                    //  { // singularity at south pole
                    //   phi = 0.0;
                    //   theta = -mathconst::PI/2;
                    //   psi = -2.0*atan2(q1,q4);                       
                    //  }
                    
                    Vec3d angles(phi, theta, psi);
                    
                    return(angles);
                    };
//------------------------------------------------------------------------------
// Vec3d TransbyQ(Vec3d& v3D, Vec4d& q)
//------------------------------------------------------------------------------
/**
 * Rotation by means of quaternion
 *
 * @param v3D         3-dimensional input vector
 * @param q           Quaternion representing the desired rotation
 *
 * @return Rotated vector
 */
//------------------------------------------------------------------------------  
Vec3d TransbyQ(Vec3d& v3D, Vec4d& q)
                    {
                    Vec3d v3D_rotated;
                    Quaterniond q_eigen, v3D_q, v3D_rotated_q;
                    
                    // Quaternion in Eigen library format. The eigen constructor for Quaterniond is Quaterniond(const Scalar &w, const Scalar &x, const Scalar &y, const Scalar &z) (w+xi+yj+zk)
                    q_eigen = Quaterniond(q(3), q(0), q(1), q(2)); 
                    
                    v3D_q.w() = 0.0;
                    v3D_q.vec() = v3D;
                    
                    //v3D_rotated_q = q_eigen*v3D_q*q_eigen.inverse(); 
                    v3D_rotated_q = q_eigen.inverse()*v3D_q*q_eigen; 
                    v3D_rotated = v3D_rotated_q.vec();
                    
                    return(v3D_rotated);
                    };
//Vec3d TransbyQ(Vec3d& v3D, Vec4d& q)
//                    {
//                    Vec3d v3D_rotated;
//                    boost::math::quaternion<double> v3D_rotated_q;
//                    boost::math::quaternion<double> q_attstate_conj;
//                    
//                    // Quaternion in Eigen library format. The eigen constructor for Quaterniond is Quaterniond(const Scalar &w, const Scalar &x, const Scalar &y, const Scalar &z) (w+xi+yj+zk)
//                    boost::math::quaternion<double> q_attstate(q(3), q(0), q(1), q(2));
//                    
//                    boost::math::quaternion<double> v3D_q(0.0, v3D(0), v3D(1), v3D(2));
//                    
//                    q_attstate_conj = boost::math::conj(q_attstate);
//                    
//                    v3D_rotated_q = q_attstate*v3D_q*q_attstate_conj;
//                    
//                    v3D_rotated(0) = v3D_rotated_q.R_component_2();
//                    v3D_rotated(1) = v3D_rotated_q.R_component_3();
//                    v3D_rotated(2) = v3D_rotated_q.R_component_4();
//                    
//                    return(v3D_rotated);
//                    };
//------------------------------------------------------------------------------
// Vec4d q_inv(Vec4d& q)
//------------------------------------------------------------------------------
/**
 * Rotation by means of quaternion
 *
 * @param q           Quaternion representing the desired rotation
 *
 * @return Inverted quaternion
 */
//------------------------------------------------------------------------------  
Vec4d q_inv(Vec4d& q)
                    {
                    Quaterniond q_eigen, q_eigen_inv;
                    Vec4d q_inverse;
                    
                    // Quaternion in Eigen library format. The eigen constructor for Quaterniond is Quaterniond(const Scalar &w, const Scalar &x, const Scalar &y, const Scalar &z) (w+xi+yj+zk)
                    q_eigen = Quaterniond(q(3), q(0), q(1), q(2));
                    
                    q_eigen_inv = q_eigen.inverse();
                    
                    q_inverse.segment(0,3) = q_eigen_inv.vec();
                    q_inverse(3) = q_eigen_inv.w();
                    
                    return(q_inverse);
                    };                       
//------------------------------------------------------------------------------
// Vec4d RotationMatrix2Quaternion(Mat3x3d& Mat)
//------------------------------------------------------------------------------
/**
 * Compute quaternion from rotation matrix
 *
 * @param Mat   Input rotation matrix
 *
 * @return Quaternion (q = q1*i + q2*j + q3*k + q4)
 */
//------------------------------------------------------------------------------                      
Vec4d RotationMatrix2Quaternion(Mat3x3d& Mat)
                    {
                    //double q1 = 0.0, q2 = 0.0, q3 = 0.0, q4 = 0.0, alpha;
                    //
                    //double testvec[4] = { Mat(0,0) + Mat(1,1) + Mat(2,2), Mat(0,0), Mat(1,1), Mat(2,2) };
                    //
                    //int maxelement = max_element(testvec, testvec + 4) - testvec;
                    //
                    //switch(maxelement)
                    //    {
                    //    case 0:
                    //          alpha = sqrt( 1.0 + Mat(0,0) + Mat(1,1) + Mat(2,2) );
                    //          
                    //          q4 = 0.5*alpha;
                    //
                    //          q1 = 0.5*( Mat(2,1) - Mat(1,2) )/alpha;
                    //                    
                    //          q2 = 0.5*( Mat(0,2) - Mat(2,0) )/alpha;
                    //                    
                    //          q3 = 0.5*( Mat(1,0) - Mat(0,1) )/alpha;
                    //          
                    //          break;
                    //          
                    //    case 1:
                    //          alpha = sqrt( 1.0 + Mat(0,0) - Mat(1,1) - Mat(2,2) );
                    //          
                    //          q4 = 0.5*( Mat(2,1) - Mat(1,2) )/alpha;
                    //
                    //          q1 = 0.5*alpha;
                    //                    
                    //          q2 = 0.5*( Mat(0,1) + Mat(1,0) )/alpha;
                    //                    
                    //          q3 = 0.5*( Mat(2,0) + Mat(0,2) )/alpha;
                    //          
                    //          break;
                    //          
                    //    case 2:
                    //          alpha = sqrt( 1.0 - Mat(0,0) + Mat(1,1) - Mat(2,2) );
                    //          
                    //          q4 = 0.5*( Mat(0,2) - Mat(2,0) )/alpha;
                    //
                    //          q1 = 0.5*( Mat(0,1) + Mat(1,0) )/alpha;
                    //                    
                    //          q2 = 0.5*alpha;
                    //                    
                    //          q3 = 0.5*( Mat(1,2) + Mat(2,1) )/alpha;
                    //          
                    //          break;
                    //          
                    //    case 3:
                    //          alpha = sqrt( 1.0 - Mat(0,0) - Mat(1,1) + Mat(2,2) );
                    //          
                    //          q4 = 0.5*( Mat(1,0) - Mat(0,1) )/alpha;
                    //
                    //          q1 = 0.5*( Mat(2,0) + Mat(0,2) )/alpha;
                    //                    
                    //          q2 = 0.5*( Mat(2,1) + Mat(1,2) )/alpha;
                    //                    
                    //          q3 = 0.5*alpha;
                    //          
                    //          break;     
                    //    }
                    //
                    //Vec4d q(q1, q2, q3, q4);
                    //
                    //return(q);
                    
                    Quaterniond q(Mat);
                    
                    Vec4d q_out;
                    
                    q_out(0) = q.x();
                    q_out(1) = q.y();
                    q_out(2) = q.z();
                    q_out(3) = q.w();
                    
                    //return(q.normalized());
                    return(q_out);
                    };
//------------------------------------------------------------------------------
// Mat3x3d Quaternion2RotationMatrix(Vec4d& q)
//------------------------------------------------------------------------------
/**
 * Compute rotation matrix from quaternion
 *
 * @param q   Input quaternion (q = q1*i + q2*j + q3*k + q4)
 *
 * @return Rotation matrix
 */
//------------------------------------------------------------------------------                         
Mat3x3d Quaternion2RotationMatrix(Vec4d& q)
                  {
                  Mat3x3d RotMat = Mat3x3d::Zero();
                  
                  double q1, q2, q3, q4;
                  //Mat3x3d SkewMat = Mat3x3d::Zero();
                  //Mat3x3d IdMat = Mat3x3d::Identity();
                  //Vec3d rotvec = Vec3d::Zero();
                  
                  q1 = q(0);  q2 = q(1);  q3 = q(2);  q4 = q(3);
                  
                  Quaterniond quat(0.0,0.0,0.0,0.0);
                  
                  quat.x() = q1;
                  quat.y() = q2;
                  quat.z() = q3;
                  quat.w() = q4;
                  
                  RotMat = quat.normalized().toRotationMatrix();
                  
                  //// Rotation vector
                  //rotvec << q1, q2, q3;
                  //
                  //// Skew-symmetric cross-product matrix SkewMat(rotvec)
                  //SkewMat(0,0) = 0.0;           SkewMat(0,1) = -rotvec(2);    SkewMat(0,2) = +rotvec(1);
                  //SkewMat(1,0) = +rotvec(2);    SkewMat(1,1) = 0.0;           SkewMat(1,2) = -rotvec(0);
                  //SkewMat(2,0) = -rotvec(1);    SkewMat(2,1) = +rotvec(0);    SkewMat(2,2) =   0.0;
                  //
                  //RotMat = IdMat + 2.0*q4*SkewMat + 2.0*SkewMat*SkewMat;
                  
                  return(RotMat);
                  };
#ifdef USE_SPICE
//------------------------------------------------------------------------------
// Vec3d ECEF2lonlath(Vec3d& posECEF)
//------------------------------------------------------------------------------
/**
 * Convert from rectangular coordinates to geodetic coordinates
 *
 * @note This function is based on the NASA SPICE library's function recgeo_c
 *
 * @param posECEF   3-dimensional position vector in ECEF frame
 *
 * @return 3-dimensional vector whose components are geodetic longitude [rad] [-pi, pi], geodetic
 *         latitude [rad] [-pi/2, pi/2] and altitude [m] of the point above the reference spheroid.
 *         The reference spheroid is defined by constants R_EARTH and F_EARTH
 *         @see Constants.h
 */
//------------------------------------------------------------------------------                   
Vec3d ECEF2lonlath(Vec3d& posECEF)
                    {
                    double pos[3], lon, lat, h;
                    for( int i = 0; i < 3 ; i++ ) pos[i] = posECEF(i);
                    
                    recgeo_c(pos, astro::R_EARTH, astro::F_EARTH, &lon, &lat, &h);
                    
                    //double x, y, z, r_xy, lon, lat, h;
                    //x = posECEF(0);  y = posECEF(1);  z = posECEF(2);
                    //lon = atan2(y,x);
                    //r_xy = sqrt(x*x + y*y);
                    //lat = atan2(z,r_xy);                    
                    //h = posECEF.norm() - astro::R_EARTH;
                    
                    Vec3d lonlath(lon,lat,h);
                    
                    return(lonlath);
                    };
//------------------------------------------------------------------------------
// Vec3d lonlath2ECEF(Vec3d& lonlath)
//------------------------------------------------------------------------------
/**
 * Convert from rectangular coordinates to geodetic coordinates
 *
 * @note This function is based on the NASA SPICE library's function georec_c
 *
 * @param lonlath   3-dimensional vector whose components are geodetic longitude [rad], geodetic
 *         latitude [rad] and altitude [m] of the point above the reference spheroid.
 *         The reference spheroid is defined by constants R_EARTH and F_EARTH
 *         @see Constants.h
 *
 * @return 3-dimensional position vector in ECEF frame
 */
//------------------------------------------------------------------------------                   
Vec3d lonlath2ECEF(Vec3d& lonlath)
                    {
                    double pos[3], lon, lat, h;
                    
                    lon = lonlath(0);
                    lat = lonlath(1);
                    h = lonlath(2);
                    
                    georec_c(lon, lat, h, astro::R_EARTH, astro::F_EARTH, pos);

                    Vec3d posECEF;
                    for( int i = 0; i < 3 ; i++ ) posECEF(i) = pos[i];
                    
                    return(posECEF);
                    };
#endif
//------------------------------------------------------------------------------
// Vec3d lonlath2El(Vec3d& lonlath, double ref_lon, double ref_lat);
//------------------------------------------------------------------------------
/**
 * Convert from rectangular coordinates to geodetic coordinates
 *
 * @param lonlath   3-dimensional vector whose components are geodetic longitude [rad], geodetic
 *                  latitude [rad] and altitude [m] of the point above the reference spheroid.
 * @param ref_lon   geodetic longitude of reference point [rad]
 * @param ref_lat   geodetic latitude of reference point [rad]
 *
 * @return 3-dimensional position vector in ECEF frame
 */
//------------------------------------------------------------------------------                   
double lonlath2El(Vec3d& lonlath,
                  double ref_lon,
                  double ref_lat)
                    {
                    double El;
                    double lon, lat, h, sinrho, deltaL, cos_lambda, sin_lambda;
                    
                    lon = lonlath(0);
                    lat = lonlath(1);
                    h = lonlath(2);
                    
                    sinrho = astro::R_EARTH/(astro::R_EARTH + h);
                    
                    deltaL = fabs( mod(lon,PI2) - mod(ref_lon,PI2) );
                                            
                    cos_lambda = sin(lat)*sin(ref_lat) + cos(lat)*cos(ref_lat)*cos(deltaL);
                    sin_lambda = sqrt( 1.0 - cos_lambda*cos_lambda);
                    
                    El = atan2( sinrho*sin_lambda , (1.0 - sinrho*cos_lambda) );
                    
                    return(El);
                    };                 
//------------------------------------------------------------------------------
// Vec3d SEZ2ECEF(Vec3d& posSEZ)
//------------------------------------------------------------------------------
/**
 * Convert from topocentric horizon coordinates (SEZ) to ECEF coordinates
 *
 * @note This function is based on the formulation given in book Fundamental of Astrodynamics and Applications by David A. Valldo
 *
 * @param posSEZ   3-dimensional position vector in SEZ frame: S points to South from the site considered, E points East and
 *                 is undefined for the North and South poles, Z is the zenith and points radially outward fro the site along
 *                 the local vertical
 * @param          geodetic longitude [rad]
 * @param          geodetic latitude [rad]
 *
 * @return         3-dimensional ECEF vector
 */
//------------------------------------------------------------------------------                   
Vec3d SEZ2ECEF(Vec3d& posSEZ,
               double lon,
               double lat)
              {
              Mat3x3d RotMat = Mat3x3d::Zero();
              Vec3d posECEF = Vec3d::Zero();
              
              double sinlon, coslon, sinlat, coslat;
              
              sinlon = sin(lon);  coslon = cos(lon);
              sinlat = sin(lat);  coslat = cos(lat);
              
              RotMat(0,0) = sinlat*coslon;    RotMat(0,1) = -sinlon;    RotMat(0,2) = coslat*coslon;
              
              RotMat(1,0) = sinlat*sinlon;    RotMat(1,1) = coslon;    RotMat(1,2) = coslat*sinlon;
              
              RotMat(2,0) = -coslat;    RotMat(2,1) = 0;    RotMat(2,2) = sinlat;
              
              posECEF = RotMat*posSEZ;
              
              return(posECEF);
              };
//------------------------------------------------------------------------------
// Vec3d ECEF2SEZ(Vec3d& posECEF, double lon, double lat)
//------------------------------------------------------------------------------
/**
 * Convert from ECEF coordinates to topocentric horizon coordinates (SEZ) 
 *
 * @note This function is based on the formulation given in book Fundamental of Astrodynamics and Applications by David A. Valldo
 *
 * @param posSEZ   3-dimensional ECEF vector
 * @param          geodetic longitude [rad]
 * @param          geodetic latitude [rad]
 *
 * @return         3-dimensional position vector in SEZ frame: S points to South from the site considered, E points East and
 *                 is undefined for the North and South poles, Z is the zenith and points radially outward from the site along
 *                 the local vertical
 */
//------------------------------------------------------------------------------                   
Vec3d ECEF2SEZ(Vec3d& posECEF,
               double lon,
               double lat)
              {
              Mat3x3d RotMat = Mat3x3d::Zero();
              Vec3d posSEZ = Vec3d::Zero();
              
              double sinlon, coslon, sinlat, coslat;
              
              sinlon = sin(lon);  coslon = cos(lon);
              sinlat = sin(lat);  coslat = cos(lat);
              
              RotMat(0,0) = sinlat*coslon;    RotMat(0,1) = -sinlon;    RotMat(0,2) = coslat*coslon;
              
              RotMat(1,0) = sinlat*sinlon;    RotMat(1,1) = coslon;    RotMat(1,2) = coslat*sinlon;
              
              RotMat(2,0) = -coslat;    RotMat(2,1) = 0;    RotMat(2,2) = sinlat;
              
              posSEZ = (RotMat.transpose())*posECEF;
              
              return(posSEZ);
              };
#ifdef USE_SPICE              
//------------------------------------------------------------------------------
// Vec3d ECEF2AzElAlt(Vector6d& stateECEF, double ref_lon, double ref_lat)
//------------------------------------------------------------------------------
/**
 * Compute the azimuth, elevation and altitude with respect to a reference point on Earth of a point vector in ECEF coordinates.
 * Elevation is defined as
 * Azimuth is defined as
 * Altitude is defined as the z component of the point vector in topocentric horizon coordinates (SEZ)
 *
 * @note This function is based on the formulation given in book Fundamental of Astrodynamics and Applications by David A. Valldo
 *
 * @param posECEF   3-dimensional ECEF vector
 * @param ref_lon   geodetic longitude of reference point [rad]
 * @param ref_lat   geodetic latitude of reference point [rad]
 *
 * @return         3-dimensional vector whose components are azimuth, elevation and altitude [rad, rad, m]
 */
//------------------------------------------------------------------------------                   
Vec3d ECEF2AzElAlt(Vector6d& stateECEF,
                  double ref_lon,
                  double ref_lat)
                 {
                  Vec3d posECEF, velECEF, rho_ECEF, rho_SEZ, d_rho_ECEF, d_rho_SEZ, rho_ref_ECEF;
                  double Az, El, Alt, sinAz, cosAz, rhoSE, rho, d_rhoSE;
                  Az = 0.0;
                  El = 0.0;
                  Alt= 0.0;
                  
                  Vec3d lonlath(ref_lon,ref_lat,0.0);
                  posECEF = stateECEF.segment(0,3);
                  velECEF = stateECEF.segment(3,3);
                  
                  //Compute ECEF vector representing the reference point on the Earth's surface
                  rho_ref_ECEF = lonlath2ECEF(lonlath);
                  //Vector from reference point on Earth to point vector (ECEF)
                  rho_ECEF = posECEF - rho_ref_ECEF;
                  d_rho_ECEF = velECEF;
                  //Transform rho_ECEF in SEZ coordinates system
                  rho_SEZ = ECEF2SEZ(rho_ECEF, ref_lon, ref_lat);
                  d_rho_SEZ = ECEF2SEZ(d_rho_ECEF, ref_lon, ref_lat);
                  
                  //Compute elevation
                  rho = rho_SEZ.norm();
                  El = asin(rho_SEZ(2)/rho);
                  
                  //Compute azimuth
                  if(El == PI/2.0)
                    {
                    d_rhoSE = sqrt( d_rho_SEZ(0)*d_rho_SEZ(0) + d_rho_SEZ(1)*d_rho_SEZ(1) );
                  
                    sinAz = d_rho_SEZ(1)/d_rhoSE;
                    cosAz = -d_rho_SEZ(0)/d_rhoSE;
                    
                    if(cosAz > 0.0 && sinAz > 0.0) Az = asin(sinAz); // North-East
                    if(cosAz > 0.0 && sinAz < 0.0) Az = PI2 + asin(sinAz); // South-East
                    if(cosAz < 0.0 && sinAz < 0.0) Az = PI - asin(sinAz); // South-West
                    if(cosAz < 0.0 && sinAz > 0.0) Az = acos(cosAz);; // North-West 
                    }
                  else
                    {
                    rhoSE = sqrt( rho_SEZ(0)*rho_SEZ(0) + rho_SEZ(1)*rho_SEZ(1) );
                    
                    sinAz = rho_SEZ(1)/rhoSE;
                    cosAz = -rho_SEZ(0)/rhoSE;
                    
                    if(cosAz > 0.0 && sinAz > 0.0) Az = asin(sinAz); // North-East
                    if(cosAz > 0.0 && sinAz < 0.0) Az = PI2 + asin(sinAz); // South-East
                    if(cosAz < 0.0 && sinAz < 0.0) Az = PI - asin(sinAz); // South-West
                    if(cosAz < 0.0 && sinAz > 0.0) Az = acos(cosAz);; // North-West
                    }
                  
                  //Compute altitude
                  Alt = rho_SEZ(2);
                  
                  Vec3d AzElAlt(Az,El,Alt);
                  
                  return(AzElAlt);
                  };
#endif
//------------------------------------------------------------------------------
// Mat3x3d ECI2RTN_Matrix(Vector6d& ECIstate)
//------------------------------------------------------------------------------
/**
 * State transformation matrix from ECI to RTN frame
 *
 * @param ECIstate     6-dimensional ECI state vector
 *
 * @return 3x3 transformation matrix from ECI to RTN
 */
//------------------------------------------------------------------------------
Mat3x3d ECI2RTN_Matrix(const Vector6d& ECIstate)
                {
                Vec3d e1, e2, e3, r_ECI, v_ECI;
                Mat3x3d T;
                Vector6d RTNstate;
                
                r_ECI = ECIstate.segment(0,3);
                v_ECI = ECIstate.segment(3,3);
                
                e1 = r_ECI/r_ECI.norm(); // Radial
                e3 = r_ECI.cross(v_ECI); e3 = e3/e3.norm(); // Cross-track
                e2 = e3.cross(e1); // Along-track
                
                T.row(0) = e1;
                T.row(1) = e2;
                T.row(2) = e3;

                return(T);
                };
//------------------------------------------------------------------------------
// Mat3x3d RTN_Matrix(const Vector6d& orbstate)
//------------------------------------------------------------------------------
/**
 * State transformation matrix from RTN frame to given reference frame
 *
 * @param ECIstate     6-dimensional state vector in given frame
 *
 * @return 3x3 transformation matrix from RTN to given frame
 */
//------------------------------------------------------------------------------
Mat3x3d RTN_Matrix(const Vector6d& orbstate) 
      {
      Vec3d e1, e2, e3, r, v;
      Mat3x3d T;
    
      // RTN frame definition
      r = orbstate.segment(0,3);
      v = orbstate.segment(3,3);
      
      // Radial
      e1 = r.normalized();
      
      // Cross-track
      e3 = r.cross(v);
      e3.normalize();
      
      // Along-track
      e2 = e3.cross(e1);
      
      // Rotation matrix
      T.col(0) = e1;
      T.col(1) = e2;
      T.col(2) = e3;
    
      return(T);
      }           
//------------------------------------------------------------------------------
// Vector6d ECI2RTN(Vector6d& ECIstate)
//------------------------------------------------------------------------------
/**
 * State transformation from ECI to RTN frame
 *
 * @param ECIstate     6-dimensional ECI state vector
 *
 * @return 6-dimensional RTN state vector
 */
//------------------------------------------------------------------------------
Vector6d ECI2RTN(Vector6d& ECIstate)
                {
                Vec3d r_ECI, v_ECI, r_RTN, v_RTN;
                Mat3x3d T;
                Vector6d RTNstate;
                
                r_ECI = ECIstate.segment(0,3);
                v_ECI = ECIstate.segment(3,3);    
                    
                T = ECI2RTN_Matrix(ECIstate);

                r_RTN = T*r_ECI;
                v_RTN = T*v_ECI;
                
                RTNstate.segment(0,3) = r_RTN;
                RTNstate.segment(3,3) = v_RTN;
                
                return(RTNstate);
                };
//------------------------------------------------------------------------------
// Vec3d RTN2ECI(Vec3d& RTNvec, Vector6d& ECIstate)
//------------------------------------------------------------------------------
/**
 * Vector transformation from RTN to ECI frame
 *
 * @param RTNvec     3-dimensional RTN vector
 * @param ECIstate   6-dimensional ECI state vector
 *
 * @return 3-dimensional ECI vector
 */
//------------------------------------------------------------------------------
Vec3d RTN2ECI(Vec3d& RTNvec, Vector6d& ECIstate)
                {
                Vec3d ECI_vec;
                Mat3x3d T, T_inv;
                
                T = ECI2RTN_Matrix(ECIstate);
                T_inv = T.transpose();
                
                ECI_vec = T_inv*RTNvec;
                
                return(ECI_vec);
                };
#ifdef USE_SPICE                
//------------------------------------------------------------------------------
// Vector6d ECEF2ECI(double GPStime, Vector6d& ECEFstate)
//------------------------------------------------------------------------------
/**
 * State transformation from ECEF (ITRF93) to ECI (J2000) frame
 *
 * @note This function is based on the NASA SPICE library's function sxform_c
 *
 * @param GPStime       GPS epoch (seconds) of the input state
 * @param ECEFstate     6-dimensional ECEF state vector
 *
 * @return 6-dimensional ECI state vector
 */
//------------------------------------------------------------------------------
Vector6d ECEF2ECI(double GPStime, Vector6d& ECEFstate)
                    {
                    Vector6d ECIstate;
                    double T[6][6];
                    Mat6x6d T_ECEF2ECI;
                    const string frame_from = "ITRF93";
                    const string frame_to = "J2000";
                    
                    double time = GPS2ET(GPStime);
                    
                    sxform_c(frame_from.c_str( ),frame_to.c_str( ),time,T);
                    
                    for( int i = 0; i < 6 ; i++ )
                        for( int j = 0; j < 6 ; j++ ) T_ECEF2ECI(i,j) = T[i][j];
                        
                    ECIstate = T_ECEF2ECI*ECEFstate;
                    
                    return(ECIstate);
                    };
//------------------------------------------------------------------------------
// Vector6d ECI2ECEF(double GPStime, Vector6d& ECIstate)
//------------------------------------------------------------------------------
/**
 * State transformation from ECI (J2000) to ECEF (ITRF93) frame
 *
 * @note This function is based on the NASA SPICE library's function sxform_c
 *
 * @param GPStime      GPS epoch (seconds) of the input state
 * @param ECIstate     6-dimensional ECI state vector
 *
 * @return 6-dimensional ECEF state vector
 * 
 */
//------------------------------------------------------------------------------
Vector6d ECI2ECEF(double GPStime, Vector6d& ECIstate)
                    {
                    Vector6d ECEFstate;
                    double T[6][6];
                    Mat6x6d T_ECI2ECEF;
                    const string frame_from = "J2000";
                    const string frame_to = "ITRF93";
                    
                    double time = GPS2ET(GPStime);
                    
                    //int lenout = 35;
                    //char utcstr[25];
                    //et2utc_c(time, "C", 0, lenout, utcstr);
                    //cout << utcstr << endl;
                    
                    sxform_c(frame_from.c_str( ),frame_to.c_str( ),time,T);
                    
                    for( int i = 0; i < 6 ; i++ )
                        for( int j = 0; j < 6 ; j++ ) T_ECI2ECEF(i,j) = T[i][j];
                        
                    ECEFstate = T_ECI2ECEF*ECIstate;
                    
                    return(ECEFstate);
                    };
//------------------------------------------------------------------------------
// Mat3x3d ECI2ECEF_Mat(double GPStime)
//------------------------------------------------------------------------------
/**
 * Position transformation matrix from ECI (J2000) to ECEF (ITRF93) frame
 *
 * @note This function is based on the NASA SPICE library's function sxform_c
 *
 * @param GPStime      GPS epoch (seconds) of the input state
 *
 * @return 3x3 state transformation matrix
 * 
 */
//------------------------------------------------------------------------------
Mat3x3d ECI2ECEF_Mat(double GPStime)
                    {
                    double T[3][3];
                    Mat3x3d T_ECI2ECEF;
                    const string frame_from = "J2000";
                    const string frame_to = "ITRF93";
                    
                    double time = GPS2ET(GPStime);
                    
                    pxform_c(frame_from.c_str( ),frame_to.c_str( ),time,T);
                    
                    for( int i = 0; i < 3 ; i++ )
                        for( int j = 0; j < 3 ; j++ ) T_ECI2ECEF(i,j) = T[i][j];
                        
                    return(T_ECI2ECEF);
                    };
#endif
//------------------------------------------------------------------------------
// Vector6d ICRF2ITRF(double GPStime, Vec3d eop, double TAI_UTC, Vector6d& ICRFstate, double dt)
//------------------------------------------------------------------------------
/**
 * State transformation from ICRF to ITRF
 *
 * @param GPStime       GPS epoch (seconds) of the input state
 * @param eop           3-dimensional vector containing Earth orientation parameters
 *                      eop(0) = CEP x-coordinate xp [1e-6 arcsec]
 *                      eop(1) = CEP y-coordinate yp [1e-6 arcsec]
 *                      eop(1) = UT1 - UTC [1e-7 sec]
 *                      The measurement units are the same of Accumulated Ultra Rapid IGS erp files (IGU)
 * @param TAI_UTC       Current TAI - UTC [s]
 * @param ICRFstate     6-dimensional ICRF state state vector
 * @param dt            Integration step for computatation of matrix d(ICRF2ITRF)/dt
 *
 * @return 6-dimensional ITRF state vector
 * 
 */
//------------------------------------------------------------------------------
Vector6d ICRF2ITRF(double GPStime, Vec3d eop, double TAI_UTC, Vector6d& ICRFstate, double dt)
                    {
                    Vector6d ITRFstate;
                    
                    Vec3d r, v, posITRF, velITRF;
                    
                    r = ICRFstate.segment(0,3); // Position
                    v = ICRFstate.segment(3,3); // Velocity
                    
                    Mat3x3d T_ICRF2ITRF, dT_ICRF2ITRF, dT_ICRF2ITRF_2, dT_ICRF2ITRF_1;
                    // Compute transformation matrices
                    T_ICRF2ITRF = ICRF2ITRF_Mat(GPStime, eop, TAI_UTC);
                    
                    dT_ICRF2ITRF_2 = ICRF2ITRF_Mat(GPStime + dt, eop, TAI_UTC);
                    dT_ICRF2ITRF_1 = ICRF2ITRF_Mat(GPStime - dt, eop, TAI_UTC);
                    
                    dT_ICRF2ITRF = (dT_ICRF2ITRF_2 - dT_ICRF2ITRF_1)/(2.0*dt);
                    
                    // Transform state vector
                    posITRF = T_ICRF2ITRF*r;
                    velITRF = dT_ICRF2ITRF*r + T_ICRF2ITRF*v;
                    
                    ITRFstate.segment(0,3) = posITRF;
                    ITRFstate.segment(3,3) = velITRF;
                    
                    return(ITRFstate);
                    };
//------------------------------------------------------------------------------
// Vector6d ITRF2ICRF(double GPStime, Vec3d eop, double TAI_UTC, Vector6d& ITRFstate, double dt)
//------------------------------------------------------------------------------
/**
 * State transformation from ITRF to ICRF
 *
 * @param GPStime       GPS epoch (seconds) of the input state
 * @param eop           3-dimensional vector containing Earth orientation parameters
 *                      eop(0) = CEP x-coordinate xp [1e-6 arcsec]
 *                      eop(1) = CEP y-coordinate yp [1e-6 arcsec]
 *                      eop(1) = UT1 - UTC [1e-7 sec]
 *                      The measurement units are the same of Accumulated Ultra Rapid IGS erp files (IGU)
 * @param TAI_UTC       Current TAI - UTC [s]
 * @param ITRFstate     6-dimensional ITRF state state vector
 * @param dt            Integration step for computatation of matrix d(ICRF2ITRF)/dt
 *
 * @return 6-dimensional ICRF state vector
 * 
 */
//------------------------------------------------------------------------------
Vector6d ITRF2ICRF(double GPStime, Vec3d eop, double TAI_UTC, Vector6d& ITRFstate, double dt)
                    {
                    Vector6d ICRFstate;
                    
                    Vec3d r, v, posICRF, velICRF;
                    
                    r = ITRFstate.segment(0,3); // Position
                    v = ITRFstate.segment(3,3); // Velocity
                    
                    Mat3x3d T_ICRF2ITRF, dT_ICRF2ITRF, dT_ICRF2ITRF_2, dT_ICRF2ITRF_1;
                    // Compute transformation matrices
                    T_ICRF2ITRF = ICRF2ITRF_Mat(GPStime, eop, TAI_UTC);
                    
                    dT_ICRF2ITRF_2 = ICRF2ITRF_Mat(GPStime + dt, eop, TAI_UTC);
                    dT_ICRF2ITRF_1 = ICRF2ITRF_Mat(GPStime - dt, eop, TAI_UTC);
                    
                    dT_ICRF2ITRF = (dT_ICRF2ITRF_2 - dT_ICRF2ITRF_1)/(2.0*dt);
                    
                    // Transform state vector
                    posICRF = T_ICRF2ITRF.transpose()*r;
                    velICRF = dT_ICRF2ITRF.transpose()*r + T_ICRF2ITRF.transpose()*v;
                    
                    ICRFstate.segment(0,3) = posICRF;
                    ICRFstate.segment(3,3) = velICRF;
                    
                    return(ICRFstate);
                    };
//------------------------------------------------------------------------------
// Mat6x6d ICRF2ITRF_Mat(double GPStime, Vec3d eop, double TAI_UTC)
//------------------------------------------------------------------------------
/**
 * State (position) transformation matrix from ICRF to ITRF
 *
 * @see Montenbruck, O., and Gill, E.,“Satellite Orbits - Model, Methods and Applications”,
 *      Springer Verlag, Heidelberg, Germany, 2000, ISBN:3-540-67280-X and
 *      https://gssc.esa.int/navipedia/index.php/Transformation_between_Celestial_and_Terrestrial_Frames
 *
 * @param GPStime      GPS epoch (seconds) of the input state
 * @param eop          3-dimensional vector containing Earth orientation parameters
 *                     eop(0) = CEP x-coordinate xp [1e6 arcsec]
 *                     eop(1) = CEP y-coordinate yp [1e6 arcsec]
 *                     eop(2) = UT1 - UTC [1e7 sec]
 *                     The measurement units are the same of Accumulated Ultra Rapid IGS erp files (IGU)
 * @param TAI_UTC      Current TAI - UTC [s]
 * 
 * @return 3-dimensional ICRF2ITRF transformation matrix
 * 
 */
//------------------------------------------------------------------------------
Mat3x3d ICRF2ITRF_Mat(double GPStime, Vec3d eop, double TAI_UTC)
                      {
                      Mat3x3d T_ICRF2ITRF;
                      
                      double xp = 1.0e-6*ARCS2RAD*eop(0); // [rad]
                      double yp = 1.0e-6*ARCS2RAD*eop(1); // [rad]
                      double UT1_UTC = 1.0e-7*eop(2); // [s]
                      
                      VectorNd<2> GPSws;
                      double TTsecs, UTCsecs, UT1secs, MJD_TT, MJD_UT1;
                      Mat3x3d T_Prec, T_Nut, T_GHA, T_Pole;
                      
                      // Compute TT seconds from GPS seconds
                      GPSws = VectorNd<2>::Zero();
                      TTsecs = GPS2TT(GPStime);
                      GPSws = GPS2GPSws(TTsecs);
                      // Compute MJD at TT seconds
                      MJD_TT = GPSws2MJD(GPSws(0), GPSws(1));
                      
                      // Compute UTC seconds from GPS seconds
                      GPSws = VectorNd<2>::Zero();
                      UTCsecs = GPS2UTC(GPStime, TAI_UTC);
                      GPSws = GPS2GPSws(UTCsecs);
                      // Compute MJD at UTC seconds
                      //MJD_UTC = GPSws2MJD(GPSws(0), GPSws(1));
                      
                      // Compute UT1 seconds from UTC seconds
                      GPSws = VectorNd<2>::Zero();
                      UT1secs = UTCsecs + UT1_UTC;
                      GPSws = GPS2GPSws(UT1secs);
                      // Compute MJD at UT1 seconds
                      MJD_UT1 = GPSws2MJD(GPSws(0), GPSws(1));
                      
                      // Precession transformation of equatorial coordinates
                      T_Prec = PrecMat(timescales::MJD_J2000, MJD_TT);
                      // Transformation matrix from mean to true equator and equinox
                      T_Nut = NutMat(MJD_TT); // dpsi is computed in function NutMat and used as input of function GHAMat
                      // Transformation matrix from true equator and equinox to Earth equator and Greenwich meridian system
                      T_GHA = GHAMat(MJD_UT1);
                      // Transformation matrix from pseudo Earth-fixed to Earth-fixed coordinates for a given date
                      T_Pole = PoleMat(xp, yp);
                      
                      // Transformation matrix ICRF2ITRF
                      T_ICRF2ITRF = T_Pole*T_GHA*T_Nut*T_Prec;
                      
                      return(T_ICRF2ITRF);
                      };
//------------------------------------------------------------------------------
// Mat3x3d ICRF2ITRF_Matrix(const double MJD_TT, const Vec3d eop, const double UTC_TAI)
//------------------------------------------------------------------------------
/**
 * State (position) transformation matrix from ICRF to ITRF
 *
 * @see Montenbruck, O., and Gill, E.,“Satellite Orbits - Model, Methods and Applications”,
 *      Springer Verlag, Heidelberg, Germany, 2000, ISBN:3-540-67280-X and
 *      https://gssc.esa.int/navipedia/index.php/Transformation_between_Celestial_and_Terrestrial_Frames
 *
 * @param Mjd_TT        Terrestrial time
 * @param eop           3-dimensional vector containing Earth orientation parameters
 *                      eop(0) = CEP x-coordinate xp [1e6 arcsec]
 *                      eop(1) = CEP y-coordinate yp [1e6 arcsec]
 *                      eop(2) = UT1 - UTC [1e7 sec]
 *                      The measurement units are the same of Accumulated Ultra Rapid IGS erp files (IGU)
 * @param UTC_TAI       Current UTC - TAI [s]
 * 
 * @return 3-dimensional ICRF2ITRF transformation matrix
 * 
 */
//------------------------------------------------------------------------------
Mat3x3d ICRF2ITRF_Matrix(const double MJD_TT, const Vec3d eop, const double UTC_TAI)
                        {
                        Mat3x3d T_ICRF2ITRF;
                        double TT_UTC, MJD_UT1, xp, yp, UT1_UTC;
                        Mat3x3d T_Prec, T_Nut, T_GHA, T_Pole;
                        
                        xp = 1.0e-6*ARCS2RAD*eop(0); // [rad]
                        yp = 1.0e-6*ARCS2RAD*eop(1); // [rad]
                        UT1_UTC = 1.0e-7*eop(2); // [s]
                        
                        TT_UTC = timescales::TT_TAI - UTC_TAI;
                        MJD_UT1 = MJD_TT + (UT1_UTC - TT_UTC)/86400.0;
                        
                        // Precession transformation of equatorial coordinates
                        T_Prec = PrecMat(timescales::MJD_J2000, MJD_TT);
                        // Transformation matrix from mean to true equator and equinox
                        T_Nut = NutMat(MJD_TT); // dpsi is computed in function NutMat and used as input of function GHAMat
                        // Transformation matrix from true equator and equinox to Earth equator and Greenwich meridian system
                        T_GHA = GHAMat(MJD_UT1);
                        // Transformation matrix from pseudo Earth-fixed to Earth-fixed coordinates for a given date
                        T_Pole = PoleMat(xp, yp);
                        
                        // Transformation matrix ICRF2ITRF
                        T_ICRF2ITRF = T_Pole*T_GHA*T_Nut*T_Prec;
                        
                        return(T_ICRF2ITRF);
                        };
//------------------------------------------------------------------------------
// Mat3x3d PrecMat(double MJD1_TT, double MJD2_TT)
//------------------------------------------------------------------------------
/**
 * Precession matrix for equatorial coordinates
 *
 * @see Montenbruck, O., and Gill, E.,“Satellite Orbits - Model, Methods and Applications”,
 *      Springer Verlag, Heidelberg, Germany, 2000, ISBN:3-540-67280-X and
 *      https://gssc.esa.int/navipedia/index.php/Transformation_between_Celestial_and_Terrestrial_Frames
 *
 * @param MJD1_TT    Epoch given (Modified Julian Date TT)
 * @param MJD2_TT    Epoch to precess to (Modified Julian Date TT)
 *
 * @return 3-dimensional precession transformation matrix
 * 
 */
//------------------------------------------------------------------------------
Mat3x3d PrecMat(double MJD1_TT, double MJD2_TT)
                {
                Mat3x3d T_Prec;
                Mat3x3d T_z, T_y, T_zeta;
                double zeta,z,theta;
                
                // Conversion from MJD to Julian centuries Terrestrial Time since J2000 TT
                const double T = (MJD1_TT - timescales::MJD_J2000)/timescales::JULIAN_DAYS_CENTURY;
                const double t = (MJD2_TT - MJD1_TT)/timescales::JULIAN_DAYS_CENTURY;
                
                // Precession angles                
                zeta = ( ( 2306.2181 + (1.39656 - 0.000139*T)*T ) + ( (0.30188 - 0.000344*T) + 0.017998*t )*t )*t*ARCS2RAD;
                
                z = zeta + ( (0.79280 + 0.000411*T) + 0.000205*t )*t*t*ARCS2RAD;
                
                theta = ( ( 2004.3109 - (0.85330 + 0.000217*T)*T ) - ( (0.42665 + 0.000217*T) + 0.041833*t )*t )*t*ARCS2RAD;
              
                // Precession matrix
                T_z = Rot_z(-z);
                T_y = Rot_y(theta);
                T_zeta = Rot_z(-zeta);
                
                T_Prec = T_z*T_y*T_zeta;
                  
                return(T_Prec);
                };
//------------------------------------------------------------------------------
// Mat3x3d NutMat(double MJD_TT)
//------------------------------------------------------------------------------
/**
 * Nutation matrix (transformation from mean-of-date to true-of-date coordinates)
 *
 * @see Montenbruck, O., and Gill, E.,“Satellite Orbits - Model, Methods and Applications”,
 *      Springer Verlag, Heidelberg, Germany, 2000, ISBN:3-540-67280-X and
 *      https://gssc.esa.int/navipedia/index.php/Transformation_between_Celestial_and_Terrestrial_Frames
 *
 * @param MJD_TT    Modified Julian Date TT
 *
 * @return 3-dimensional nutation transformation matrix
 * 
 */
//------------------------------------------------------------------------------
Mat3x3d NutMat(double MJD_TT)
                {
                Mat3x3d T_nut;
                Mat3x3d T_epsdeps, T_dpsi, T_eps;
                Matrix<long, 106, 9> NutCoeff = Matrix<long, 106, 9>::Zero();
                double T, eps, deps, dpsi;                
                deps = 0.0;
                dpsi = 0.0;
                
                // Conversion from MJD to Julian centuries Terrestrial Time since J2000 TT
                T = (MJD_TT - timescales::MJD_J2000)/timescales::JULIAN_DAYS_CENTURY;
                
                ////////////// Mean obliquity of the ecliptic ////////////
                eps = ( 23.43929111 - ( 46.8150 + (0.00059 - 0.001813*T)*T )*T/3600.0 )*DEG2RAD;
                
                ////////////// Nutation in longitude and obliquity ////////////
                // Get Coefficients
                get_NutCoeff(NutCoeff);
                // Compute nutation angles
                NutationAngles(MJD_TT, NutCoeff, dpsi, deps);
                
                //////////// Nutation matrix ////////////                
                T_epsdeps = Rot_x(-eps -deps);                
                T_dpsi = Rot_z(-dpsi);                
                T_eps = Rot_x(+eps);
                
                T_nut = T_epsdeps*T_dpsi*T_eps;
                  
                return(T_nut);
                };
//------------------------------------------------------------------------------
// void NutationAngles(double MJD, Matrix<long, 106, 9> NutCoeff, double& dpsi, double& deps)
//------------------------------------------------------------------------------
/**
 * Nutation in longitude and obliquity
 *
 * @see Montenbruck, O., and Gill, E.,“Satellite Orbits - Model, Methods and Applications”,
 *      Springer Verlag, Heidelberg, Germany, 2000, ISBN:3-540-67280-X and
 *      https://gssc.esa.int/navipedia/index.php/Transformation_between_Celestial_and_Terrestrial_Frames
 *
 * @param MJD       Modified Julian Date
 * @param NutCoeff  Matrix of IAU 1980 nutation theory coefficients
 *
 * @return dpsi and deps   Nutation in longitude and obliquity [rad]
 * 
 */
//------------------------------------------------------------------------------
void NutationAngles(double MJD, Matrix<long, 106, 9> NutCoeff, double& dpsi, double& deps)
                {
                deps = 0.0;
                dpsi = 0.0;
                
                //Matrix<long, 106, 9> NutCoeff = Matrix<long, 106, 9>::Zero();
                const int N = 106;
                const double rev = 360.0*3600.0; // [arcsec/revolution]
                
                double T, T2, T3;
                double  l, lp, F, D, Om;
                double  arg;
                
                // Conversion from MJD to Julian centuries Terrestrial Time since J2000 TT
                T = (MJD - timescales::MJD_J2000)/timescales::JULIAN_DAYS_CENTURY;
                T2 = T*T;
                T3 = T2*T;
                
                // Mean anomaly of the Moon
                l = mod(485866.733 + (1325.0*rev + 715922.633)*T + 31.310*T2 + 0.064*T3, rev);
                // Mean anomaly of the Sun
                lp = mod(1287099.804 + ( 99.0*rev + 1292581.224)*T - 0.577*T2 - 0.012*T3, rev);
                // Mean argument of latitude
                F = mod(335778.877 + (1342.0*rev +  295263.137)*T - 13.257*T2 + 0.011*T3, rev);
                // Mean longitude elongation of the Moon from the Sun
                D = mod(1072261.307 + (1236.0*rev + 1105601.328)*T - 6.891*T2 + 0.019*T3, rev );
                // Mean longitude of the ascending node
                Om = mod(450160.280 - ( 5.0*rev +  482890.539)*T + 7.455*T2 + 0.008*T3, rev);
                
                //get_NutCoeff(NutCoeff);

                // Nutation in longitude and obliquity [rad]
                for(int i = 0; i < N; i++)
                    {
                    arg = ( NutCoeff(i,0)*l + NutCoeff(i,1)*lp + NutCoeff(i,2)*F + NutCoeff(i,3)*D + NutCoeff(i,4)*Om )*ARCS2RAD;
                    
                    dpsi += ( NutCoeff(i,5) + NutCoeff(i,6)*T )*sin(arg);
                  
                    deps += ( NutCoeff(i,7) + NutCoeff(i,8)*T )*cos(arg);
                    };
                    
                dpsi = 1.0e-5*dpsi*ARCS2RAD;
                deps = 1.0e-5*deps*ARCS2RAD;
                };
//------------------------------------------------------------------------------
// Mat3x3d GHAMat(double MJD_UT1)
//------------------------------------------------------------------------------
/**
 * Transformation from true equator and equinox to Earth equator and Greenwich meridian system
 *
 * @see Montenbruck, O., and Gill, E.,“Satellite Orbits - Model, Methods and Applications”,
 *      Springer Verlag, Heidelberg, Germany, 2000, ISBN:3-540-67280-X and
 *      https://gssc.esa.int/navipedia/index.php/Transformation_between_Celestial_and_Terrestrial_Frames
 *
 * @param MJD    Modified Julian Date TT
 *
 * @return 3-dimensional Greenwich Hour Angle transformation matrix
 * 
 */
//------------------------------------------------------------------------------
Mat3x3d GHAMat(double MJD_UT1)
                {
                Mat3x3d T_GHA;
                
                double MJD0, UT1, T0 , T, gmst, GMST_MJD_UT1_;//intpart
                double GMST_MJD_UT1, MeanObliquity, EqnEquinox, GAST_MJD_UT1;
                double deps = 0.0, dpsi = 0.0;
                Matrix<long, 106, 9> NutCoeff = Matrix<long, 106, 9>::Zero();

                // Compute Greenwich Mean Sidereal Time
                MJD0 = floor(MJD_UT1);
                UT1 = timescales::JULIAN_DAY*(MJD_UT1 - MJD0); // [s]
                T0 = (MJD0 - timescales::MJD_J2000)/timescales::JULIAN_DAYS_CENTURY; 
                T = (MJD_UT1 - timescales::MJD_J2000)/timescales::JULIAN_DAYS_CENTURY; 

                gmst  = 24110.54841 + 8640184.812866*T0 + 1.002737909350795*UT1 + (0.093104 - 6.2e-6*T)*T*T;// [s]
                
                GMST_MJD_UT1_ = gmst/timescales::JULIAN_DAY;
                GMST_MJD_UT1 = PI2*(GMST_MJD_UT1_ - floor(GMST_MJD_UT1_));//PI2*modf(GMST_MJD_UT1_, &intpart);// [rad], [0,2PI]
                
                // Computation of the equation of the equinoxes: the equation of the equinoxes dpsi*cos(eps) is the right ascension
                // of the mean equinox referred to the true equator and equinox and is equal to the difference between apparent and
                // mean sidereal time.                
                MeanObliquity = ( 23.43929111 - ( 46.8150 + (0.00059 - 0.001813*T)*T )*T/3600.0 )*DEG2RAD;
                
                // Get Coefficients
                get_NutCoeff(NutCoeff);
                // Compute nutation angles
                NutationAngles(MJD_UT1, NutCoeff, dpsi, deps);
                
                // Equation of the equinoxes: the equation of the equinoxes dpsi*cos(eps) is the right ascension of the mean equinox
                // referred to the true equator and equinox and is equal to the difference between apparent and mean sidereal time
                EqnEquinox = dpsi*cos(MeanObliquity);
                
                // Computation of Greenwich Apparent Sidereal Time
                GAST_MJD_UT1 = mod( GMST_MJD_UT1 + EqnEquinox, PI2 );
                
                // Greenwich Hour Angle transformation matrix
                T_GHA = Rot_z(GAST_MJD_UT1);
                
                return(T_GHA);
                };
//------------------------------------------------------------------------------
// Mat3x3d PoleMat(double xp, double yp)
//------------------------------------------------------------------------------
/**
 * Transformation from true equator and equinox to Earth equator and Greenwich meridian system
 *
 * @see Montenbruck, O., and Gill, E.,“Satellite Orbits - Model, Methods and Applications”,
 *      Springer Verlag, Heidelberg, Germany, 2000, ISBN:3-540-67280-X and
 *      https://gssc.esa.int/navipedia/index.php/Transformation_between_Celestial_and_Terrestrial_Frames
 *
 * @param   xp, yp define the position of the Conventional Ephemeris Pole (CEP) with respect to
 *          the Conventional Terrestrial Pole (CTP)
 *
 * @return 3-dimensional polar motion matrix
 * 
 */
//------------------------------------------------------------------------------
Mat3x3d PoleMat(double xp, double yp)
              {                
              Mat3x3d T_Pole;
              //Mat3x3d T_xp, T_yp;
              //
              //T_xp = Rot_y(-xp);
              //T_yp = Rot_x(-yp);
              //
              //T_Pole = T_xp*T_yp;
              
              T_Pole = Rot_y(-xp)*Rot_x(-yp);
              
              return(T_Pole);
              };
#ifdef USE_SPICE              
//------------------------------------------------------------------------------
// Vec3d ECI2RAD(Vec3d& posECI)
//------------------------------------------------------------------------------
/**
 * Convert position vector from rectangular coordinates ECI to right ascension, declination and range ECI
 *
 * @param posECI   3-dimensional position vector in ECI
 *
 * @return 3-dimensional position vector as right ascension, declination and range ECI
 */
//------------------------------------------------------------------------------                   
Vec3d ECI2RAD(Vec3d& posECI)
            {
            Vec3d posRAD = Vec3d::Zero();
            double posECI_arr[3];
            double ra, dec, range;
                          
            for( int i = 0; i < 3 ; i++ ) posECI_arr[i] = posECI(i);
            
            // Convert in right ascension, declination and range, 
            recrad_c(posECI_arr, &range, &ra, &dec);
            
            posRAD(0) = ra;
            posRAD(1) = dec;
            posRAD(2) = range;
                            
            return(posRAD);
            };                    
//------------------------------------------------------------------------------------------------
// Vec3d v3D_transform(double GPStime, Vec3d& v3D, const string frame_from, const string frame_to)
//------------------------------------------------------------------------------------------------
/**
 * Transformation of a generic 3-dimensional vector between two frames
 *
 * @note This function is based on the NASA SPICE library's function pxform_c
 *
 * @param GPStime       GPS epoch (seconds) of the input vector
 * @param v3D           3-dimensional input vector
 * @param frame_from    String containing the name of the starting frame. See SPICE documentation
 *                      for the names of the built-in inertial and body-fixed frames
 * @param frame_to      String containing the name of the end frame
 *
 * @return 6-dimensional ECEF state vector
 */
//------------------------------------------------------------------------------------------------
Vec3d v3D_transform(double GPStime, Vec3d& v3D, const string frame_from, const string frame_to)
                    {
                    Vec3d v3D_transformed;
                    double T[3][3];
                    Mat3x3d Transmat;
                    
                    double time = GPS2ET(GPStime);
                    
                    pxform_c(frame_from.c_str( ),frame_to.c_str( ),time,T);
                    
                    for( int i = 0; i < 3 ; i++ )
                        for( int j = 0; j < 3 ; j++ ) Transmat(i,j) = T[i][j];
                        
                    v3D_transformed = Transmat*v3D;
                    
                    return(v3D_transformed);
                    };
#endif
//------------------------------------------------------------------------------------------------
// Vec3d u_axis(Mat3x3d& T1toT2, const string axis)
//------------------------------------------------------------------------------------------------
/**
 * Gives unit vector of axis x1, y1 or z1-direction, expressed in x2, y2 and z2 coordinates given
 * matrix T1toT2 and the desired axis name x, y or z
 *
 * @param T1toT2       Transformation matrix from frame 1 to frame 2
 * @param axis_name   Name of desired axis "x", "y" or "z"
 *
 * @return Unit vecor of frame 1's principal direction x, y or z of expressed in frame 2's coordinates
 */
//------------------------------------------------------------------------------------------------
Vec3d u_axis(Mat3x3d& T1toT2, const string axis_name)
                    {
                    Vec3d axis;
                    
                    if(axis_name.compare("x") == 0) axis = T1toT2.col(0);
                    else if(axis_name.compare("y") == 0) axis = T1toT2.col(1);
                    else if(axis_name.compare("z") == 0) axis = T1toT2.col(2);
                    else if(axis_name.compare("x") || axis_name.compare("y") || axis_name.compare("z"))
                            {
                            cout << "Function parameter 2 has to be string \"x\", \"y\" or \"z\"" << endl;
                            return(Vec3d::Zero());
                            }
                            
                    return(axis);
                    };
//------------------------------------------------------------------------------
// Vector6d rv2Kepler(Vector6d& ECIstate)
//------------------------------------------------------------------------------
/**
 * Conversion of Cartesian state vector into osculating orbital elements for non-circular non-equatorial orbits
 * @note The function cannot be used with state vectors describing a circular or non-inclined orbit.
 *
 * @param ECIstate   6-dimensional ECI state vector
 *
 * @return Validity flag
 * @return 6-dimensional Keplerian elements (a,e,i,Omega,omega,M) with
 *         a      Semimajor axis [m]
 *         e      Eccentricity 
 *         i      Inclination [rad]
 *         Omega  Longitude of the ascending node [rad]
 *         omega  Argument of pericenter  [rad]
 *         M      Mean anomaly  [rad]
 *         
 */
//------------------------------------------------------------------------------
Vector6d rv2Kepler(Vector6d& ECIstate, bool& valid)
                  {
                  Vector6d KepElem = Vector6d::Zero();
                  valid = true;
                  
                  Vec3d r, v, h;
                  double H, u, R;
                  double eCosE, eSinE, e2, E, nu;
                  double a, e, i, Omega, omega, M;
                
                  // Position
                  r = ECIstate.segment(0,3);
                  
                  // Velocity
                  v = ECIstate.segment(3,3);
                  
                  // Areal velocity
                  h = r.cross(v);
                  H = h.norm();
                  
                  // Longitude of ascending node
                  Omega = atan2( h(0), -h(1) );
                  Omega = mod(Omega,PI2);
                  
                  // Inclination
                  i = atan2( sqrt(h(0)*h(0) + h(1)*h(1)), h(2) );
                  
                  // Argument of latitude
                  u = atan2( r(2)*H, -r(0)*h(1) + r(1)*h(0) );   
                
                  // Distance
                  R = r.norm();
                  
                  // Error handling
                  if( ( 2.0/R - v.dot(v)/GM_EARTH ) <= 0.0 )
                    {
                    valid = false;
                    return(KepElem); // Null vector
                    }
  
                  // Semi-major axis 
                  a = 1.0/( 2.0/R - v.dot(v)/GM_EARTH );
                
                  // e*cos(E) 
                  eCosE = 1.0 - R/a;
                  
                  // e*sin(E)
                  eSinE = r.dot(v)/sqrt(GM_EARTH*a);
                  
                  // Eccentricity
                  e2 = eCosE*eCosE + eSinE*eSinE;
                  e = sqrt(e2);
                  
                  // Eccentric anomaly
                  E = atan2(eSinE, eCosE);
                  
                  // Mean anomaly
                  M = mod(E - eSinE,PI2);
                
                  // True anomaly
                  nu = atan2( sqrt(1.0 - e2)*eSinE, eCosE - e2 );
                  
                  // Argument of perihelion
                  omega = mod(u - nu,PI2);
                 
                  // Keplerian elements vector
                  KepElem << a, e, i, Omega, omega, M;
                  
                  return(KepElem);
                  };          
//------------------------------------------------------------------------------
// Vector6d Kepler2rv(Vector6d& KepElem, bool& valid)
//------------------------------------------------------------------------------
/**
 * Conversion of osculating Keplerian elements for non-circular non-equatorial orbits to Cartesian state vector
 * @note The function cannot be used with state vectors describing a circular or non-inclined orbit.
 *
 * @param  6-dimensional Keplerian elements (a,e,i,Omega,omega,M) with
 *         a        Semimajor axis [m]
 *         e        Eccentricity 
 *         i        Inclination [rad]
 *         Omega    Longitude of the ascending node [rad]
 *         omega    Argument of pericenter  [rad]
 *         M        Mean anomaly  [rad]
 * @param  dt       Time since epoch [s]
 *         
 * @return Validity flag
 * @return 6-dimensional ECI state vector
 *         
 */
//------------------------------------------------------------------------------
Vector6d Kepler2rv(Vector6d& KepElem, double dt, bool& valid)
                    {
                    Vector6d ECIstate = Vector6d::Zero();
                    valid = true;
                    
                    double  a, e, i, Omega, omega, M, M0, n;
                    double  E, cosE, sinE, fac, R, V;
                    Vec3d r, v;
                    Mat3x3d PQW;
  
                    // Keplerian elements at epoch
                    a = KepElem(0);
                    e = KepElem(1);
                    i = KepElem(2);
                    Omega = KepElem(3);
                    omega = KepElem(4);
                    M0 = KepElem(5);
  
                    // Mean anomaly
                    if(dt == 0.0) M = M0;
                    else
                        {
                        n = sqrt( GM_EARTH/(a*a*a) );
                        M = M0 + n*dt;
                        };

                    // Eccentric anomaly
                    E = EccAnomaly(M, e, 15, 1.0E-14);
                  
                    cosE = cos(E); 
                    sinE = sin(E);

                    // Perifocal coordinates
                    fac = sqrt( (1.0 - e)*(1.0 + e) );
                    
                    // Distance
                    R = a*(1.0 - e*cosE);
                    // Velocity
                    V = sqrt(GM_EARTH*a)/R;
                  
                    r << a*(cosE - e), a*fac*sinE , 0.0;
                    v << -V*sinE, +V*fac*cosE, 0.0;

                    // Transformation to reference system (Gaussian vectors)
                    PQW = Rot_z(-Omega)*Rot_x(-i)*Rot_z(-omega);
                  
                    r = PQW*r;
                    v = PQW*v;
                    
                    ECIstate.segment(0,3) = r;
                    ECIstate.segment(3,3) = v;
                  
                    return(ECIstate);
                    };          
//------------------------------------------------------------------------------
// Vector6d posvecs2Kepler(Vec4d& timepos1, Vec4d& timepos2, bool& valid)
//------------------------------------------------------------------------------
/**
 * Computation of Keplerian orbital elements from two given position vectors and associated times
 * @note The function cannot be used with state vectors describing a circular or non-inclined orbit.
 *
 * @param timepos1   4-dimensional vector: timepos1(0) = epoch of pos1 (Modified Julian Date), timepos1(1) to timepos1(3) = position vector components
 * @param timepos2   4-dimensional vector: timepos2(0) = epoch of pos2 (Modified Julian Date), timepos2(1) to timepos2(3) = position vector components
 * 
 * @return Validity flag
 * @return 6-dimensional Keplerian elements (a,e,i,Omega,omega,M) at epoch 1 with
 *         a      Semimajor axis [m]
 *         e      Eccentricity 
 *         i      Inclination [rad]
 *         Omega  Longitude of the ascending node [rad]
 *         omega  Argument of pericenter  [rad]
 *         M      Mean anomaly  [rad]
 *         
 */
//------------------------------------------------------------------------------
Vector6d posvecs2Kepler(Vec4d& timepos1, Vec4d& timepos2, bool& valid)
                        {
                        Vector6d KepElem = Vector6d::Zero();
                        valid = true;
                        
                        double t1_MJD, t2_MJD;
                        Vec3d r_1, r_2;
                        
                        t1_MJD = timepos1(0);
                        r_1 = timepos1.segment(1,3);
                        
                        t2_MJD = timepos2(0);
                        r_2 = timepos2.segment(1,3);
                        
                        double tau, eta, p;
                        double nu, E, u;
                        double s_1, s_2, s_0, fac, sinhH;
                        double cos_dnu, sin_dnu, ecos_nu, esin_nu;
                        double a, e, i, Omega, omega, M;
                        Vec3d e_1, r_0, e_0, W;
  
                        // Compute vector r_0 (fraction of r_2 perpendicular to r_1) the norms of r_1, r_2 and r_0
                        s_1 = r_1.norm();
                        e_1 = r_1/s_1;
                        
                        s_2 = r_2.norm();
                        fac = r_2.dot(e_1);
                        
                        r_0 = r_2 - fac*e_1;
                        s_0 = r_0.norm();
                        e_0 = r_0/s_0;
  
                        // Inclination and ascending node
                        W = e_1.cross(e_0);
                        
                        // Longitude of ascending node
                        Omega = atan2( W(0), -W(1) );
                        Omega = mod(Omega, PI2);
                        
                        // Inclination
                        i = atan2( sqrt( W(0)*W(0) + W(1)*W(1) ), W(2) );
                        if (i == 0.0) u = atan2( r_1(1), r_1(0) );
                        else u = atan2 ( +e_1(2) , -e_1(0)*W(1) + e_1(1)*W(0) );
  
                        // Semilatus rectum
                        tau = sqrt(GM_EARTH)*86400.0*fabs(t2_MJD - t1_MJD);   
                        eta = posvecs2Eta(r_1, r_2, tau, valid);
                        p   = (s_1*s_0*eta/tau)*(s_1*s_0*eta/tau);//pow ( s_1*s_0*eta/tau, 2 );   

                        // Eccentricity, true anomaly and argument of perihelion
                        cos_dnu = fac/s_2;    
                        sin_dnu = s_0/s_2;
                      
                        ecos_nu = p/s_1 - 1.0;  
                        esin_nu = ( ecos_nu*cos_dnu - (p/s_2 - 1.0) )/sin_dnu;
                      
                        e = sqrt( ecos_nu*ecos_nu + esin_nu*esin_nu );
                        nu = atan2(esin_nu,ecos_nu);
                      
                        omega = mod(u - nu,PI2);

                        // Perihelion distance, semimajor axis and mean motion
                        a = p/(1.0 - e*e);
                        //n = sqrt( GM_EARTH/fabs(a*a*a) );

                        // Mean anomaly and time of perihelion passage
                        if(e < 1.0)
                          {
                          E = atan2( sqrt( (1.0 - e)*(1.0 + e) )*esin_nu, ecos_nu + e*e );
                          M = mod( E - e*sin(E), PI2 );
                          }
                        else 
                          {
                          sinhH = sqrt( (e - 1.0)*(e + 1.0) )*esin_nu/( e + e*ecos_nu );
                          M = e*sinhH - log( sinhH + sqrt(1.0 + sinhH*sinhH) );
                          }

                        // Keplerian elements vector
                        KepElem << a, e, i, Omega, omega, M;
                        
                        return(KepElem);
                        };
// Local function to be used inside function posvecs2Eta
double secant_eta(double eta, double m, double l)
                    {
                    const double delta = 100.0*EPS_MACHINE_D;
                    const int maxit = 30;
                    double F = 0.0;
                    double w, W, a, n, g;
                    int i = 0;
  
                    w = m/(eta*eta) - l; 
  
                    if(fabs(w) < 0.1) // Series expansion
                      { 
                      a = 4.0/3.0;
                      W = a;
                      n = 0.0;
                      do
                        {
                        n += 1.0;
                        a *= w*(n + 2.0)/(n + 1.5);
                        W += a;
                        
                        ++i;
                      
                        if( i == maxit )
                              {
                              //valid = false; // Convergence problems in FindEta
                              break;
                              }
                        }
                        while(fabs(a) >= delta);
                      }
                    else
                      {
                      if(w > 0.0)
                        {
                        g = 2.0*asin(sqrt(w));  
                        W = (2.0*g - sin(2.0*g))/( sin(g)*sin(g)*sin(g) );//pow(sin(g), 3);
                        }
                      else
                        {
                        g = 2.0*log(sqrt(-w) + sqrt(1.0-w));  // =2.0*arsinh(sqrt(-w))
                        W = (sinh(2.0*g) - 2.0*g)/( sin(g)*sin(g)*sin(g) );//pow(sinh(g), 3);
                        }
                      }
                      
                    F = 1.0 - eta + (w + l)*W;
                    
                    return(F);
                    };                        
//------------------------------------------------------------------------------
// double posvecs2Eta(Vec3d& r_1, Vec3d& r_2, double tau, bool& valid)
//------------------------------------------------------------------------------
/**
 * Computation of the sector-triangle ratio from two position vectors and the intermediate time
 *
 * @param r_1     Position vector at epoch t_1
 * @param r_2     Position vector at epoch t_2
 * @param tau     Normalized time sqrt(GM_EARTH)*(t_2 - t_1) where t_2 - t_1 is in seconds
 * 
 * @return Validity flag
 * @return Sector-triangle ratio
 * 
 */
//------------------------------------------------------------------------------
double posvecs2Eta(Vec3d& r_1, Vec3d& r_2, double tau, bool& valid)
                    {
                    double eta2 = 0.0;
                    valid = true;
                    
                    const int maxit = 30;
                    const double delta = 100.0*EPS_MACHINE_D;
  
                    int i;
                    double kappa, m, l, s_1, s_2, eta_min, eta1, F1, F2, d_eta;
  
                    // Auxiliary quantities
                    s_1 = r_1.norm();
                    s_2 = r_2.norm();
                  
                    kappa = sqrt( 2.0*( s_1*s_2 + r_1.dot(r_2) ) );
                  
                    m = tau*tau/(kappa*kappa*kappa);  
                    l = (s_1 + s_2)/(2.0*kappa) - 0.5;
                  
                    eta_min = sqrt( m/(l + 1.0) );
  
                    // Starting values computed with Hansen's approximation
                    eta2 = ( 12.0 + 10.0*sqrt( 1.0 + (44.0/9.0)*m/(l + 5.0/6.0) ) )/22.0;
                    eta1 = eta2 + 0.1;
                    
                    // Secant method
                    F1 = secant_eta(eta1, m, l);   
                    F2 = secant_eta(eta2, m, l);
                    
                    i = 0;
                  
                    while( fabs(F2 - F1) > delta )
                        {
                        d_eta = -F2*(eta2 - eta1)/(F2 - F1);  
                        eta1 = eta2;
                        F1 = F2;
                        
                        while( (eta2 + d_eta) <= eta_min ) d_eta *= 0.5;
                        eta2 += d_eta;  
                        F2 = secant_eta(eta2,m,l);
                        
                        ++i;
                      
                        if( i == maxit )
                          {
                          valid = false; // Convergence problems in FindEta
                          break;
                          }
                        }
                    
                    return(eta2);
                    };
//------------------------------------------------------------------------------
// Vector6d rv2oe(Vector6d& ECIstate, bool& valid)
//------------------------------------------------------------------------------
/**
 * Conversion of Cartesian state vector into osculating orbital elements for non-circular non-equatorial orbits
 *
 * @param ECIstate   6-dimensional ECI state vector
 * @param valid      Validity flag
 *
 * @return 6-dimensional orbital elements (a,ex,ey,i,Omega,u) with
 *         a      Semimajor axis [m]
 *         ex     x-component eccentricity vector
 *         ey     y-component eccentricity vector
 *         i      Inclination [rad]
 *         Omega  Longitude of the ascending node [rad]
 *         u      Mean argument of latitude [rad]
 */
//------------------------------------------------------------------------------
Vector6d rv2oe(Vector6d& ECIstate, bool& valid)
                {
                Vector6d orbel = Vector6d::Zero();
                Vector6d orbelout = Vector6d::Zero();
                double ex,ey;       // Components of the eccentricity vector
                double u;           // Mean argument of latitude
                bool   valid_elem;  // Validity flag
            
                // Call "rv2Kepler" function (not valid for circular orbits, e=0)
                orbel = rv2Kepler(ECIstate,valid_elem);
            
                // Error handling
                if(!valid_elem)
                  {
                  valid = false;
                  return(orbelout); // Null vector
                  };
            
                // Compute Eccentricity Vector
                ex = orbel(1)*cos(orbel(4));
                ey = orbel(1)*sin(orbel(4));
            
                // Compute mean argument of latitude
                u = orbel(5)+orbel(4);
                u = mod(u,PI2);
            
                // Compute orbital elements vector
                orbelout << orbel(0), ex, ey, orbel(2), orbel(3), u;
            
                valid = true;
                
                return(orbelout);    
                };
//------------------------------------------------------------------------------
// Vector6d osc2mean(Vector6d& osc_elem, bool& valid)
//------------------------------------------------------------------------------
/**
 * Conversion of osculating orbital elements into mean orbital elements
 * @note A first-order mapping algorithm is applied developed by Brouwer and Lyddane
 * @see Brouwer D., Solution of the Problem of Artificial Satellite Theory without Drag,
 *      Astronautical Journal,Vol. 64,No. 1274, 1959,pp. 378-397.
 * @see Lyddane R.H., Small Eccentricities or Inclinations in the Brouwer Theory of the
 *      Artificial Satellite, Astronautical Journal, Vol. 68, No. 8, 1963, pp. 555-558.
 *
 * @param ECIstate   6-dimensional osculating orbital elements (a,ex,ey,i,Omega,u) vector
 * @param valid      Validity flag
 *
 * @return 6-dimensional mean orbital elements (a,ex,ey,i,Omega,u) vector with
 *         a      Semimajor axis [m]
 *         ex     x-component eccentricity vector
 *         ey     y-component eccentricity vector
 *         i      Inclination [rad]
 *         Omega  Longitude of the ascending node [rad]
 *         u      Mean argument of latitude [rad]
 */
//------------------------------------------------------------------------------
Vector6d osc2mean(Vector6d& osc_elem, bool& valid)
                    {
                    Vector6d mean_elem = Vector6d::Zero();
                    
                    double a,ex,ey,i,O,u,e_2,e,w,M,gamma2,eta,eta_2,eta_3,eta_4,eta_6;
                    double gammap2,E,f,c_f,c_f_2,c_f_3,s_f,c_i,c_i_2,c_i_4,c_i_6,t_i;
                    double c_2w2f,c_2w,c_2wf,c_2w3f,s_2w2f,s_2wf,s_2w3f;
                    double c_M,s_M,c_O,s_O,s_ix2,c_ix2,axr,axr_2,axr_3;
                    double threeci2m1,onemci2,ac1,ac2,ac,am,de1c,de1,sef,sef1,sef2;
                    double dec1,dec2,de,sewf,di,onem5ci2,onem5ci22,fmMpesf,threem5ci2,cewf;
                    double L1,L21,L22,L23,L2,L31,L32,L3,L41,L4,L5,L;
                    double saxr1,saxr2,saxr3,edM1,edM21,edM22,edM2,edM,dO,d1,d2,Mm,em;
                    double scix2,d3,d4,Om,im,wm,exm,eym,um;

                    //  Get osculating orbital elements from input vector
                    a  = osc_elem(0); 
                    ex = osc_elem(1); 
                    ey = osc_elem(2); 
                    i  = osc_elem(3);
                    O  = osc_elem(4);
                    u  = osc_elem(5);
                
                    //  Compute standard orbital elements
                    e_2 = ex*ex + ey*ey;
                    e   = sqrt(e_2);
                    w   = atan2(ey,ex);
                    w   = mod(w,PI2);
                    M   = u - w;
                    M   = mod(M,PI2);

                    //  Error handling
                    if(a==0.0 || e>=1.0 || e<=0.0)
                      {
                      valid = false;
                      return(mean_elem); // Null vector
                      };

                    //  Algorithm parameters
                    gamma2  = -(J2/2.0)*(R_EARTH*R_EARTH)/(a*a);
                    eta     = sqrt(1.0 - (e*e));
                    eta_2   = 1.0 - (e*e); 
                    eta_4   = eta_2*eta_2;
                    eta_3   = eta*eta*eta;
                    eta_6   = eta_3*eta_3;
                    gammap2 = gamma2/eta_4;

                    //  Eccentric anomaly
                    E = EccAnomaly(M, e, 15, 1.0E-14);

                    //  True anomaly
                    f = 2.0*atan(sqrt((1.0 + e)/(1.0 - e))*tan(E/2.0));
                    f = mod(f,PI2);

                    //  Used cosines (c_), sins (s_), tangents (t_) and powers (*_*)
                    c_f    = cos(f);
                    c_f_2  = c_f*c_f; 
                    c_f_3  = c_f_2*c_f; 
                    s_f    = sin(f);
                    c_i    = cos(i); 
                    c_i_2  = c_i*c_i; 
                    c_i_4  = c_i_2*c_i_2; 
                    c_i_6  = c_i_4*c_i_2;
                    t_i    = tan(i); 
                    c_2w2f = cos(2.0*w + 2.0*f); 
                    c_2w   = cos(2.0*w); 
                    c_2wf  = cos(2.0*w + f);
                    c_2w3f = cos(2.0*w + 3.0*f); 
                    s_2w2f = sin(2.0*w + 2.0*f);
                    s_2wf  = sin(2.0*w + f); 
                    s_2w3f = sin(2.0*w + 3.0*f);
                    c_M    = cos(M); 
                    s_M    = sin(M); 
                    c_O    = cos(O); 
                    s_O    = sin(O);
                    s_ix2  = sin(i/2.0); 
                    c_ix2  = cos(i/2.0);

                    //  Error handling
                    if( c_i_2 == 1.0/5.0 || eta == -1.0)
                      {
                      valid = false;
                      return(mean_elem); // Null vector
                      };

                    //  a/r
                    axr   = (1.0 + (e*c_f))/eta_2;
                    axr_2 = axr*axr; 
                    axr_3 = axr_2*axr;
                
                    //  Mean semimajor axis
                    threeci2m1 = (3.0*c_i_2) - 1.0;
                    onemci2    = 1.0 - c_i_2;
                    ac1        = threeci2m1*(axr_3 - (1.0/eta_3));
                    ac2        = 3.0*onemci2*axr_3*c_2w2f;
                    ac         = ac1 + ac2;
                    am         = a + a*gamma2*ac;

                    //  Intermediate parameters
                    de1c = 1.0 - (11.0*c_i_2) - (40.0*(c_i_4/(1.0 - (5.0*c_i_2))));
                    de1  = (1.0/8.0)*gammap2*e*eta_2*de1c*c_2w;
                
                    sef  = (3.0*c_f) + (3.0*e*c_f_2) + (e_2*c_f_3);
                    sef1 = (e*eta) + (e/(1.0+eta)) + sef;
                    sef2 = e + sef;
                    dec1 = ((threeci2m1/eta_6)*sef1) + ((3.0*onemci2/eta_6)*sef2*c_2w2f);
                    dec2 = onemci2*(3.0*c_2wf + c_2w3f);
                    de   = de1 + (eta_2/2.0)*(gamma2*dec1 - gammap2*dec2);

                    sewf = (3.0*c_2w2f) + (3.0*e*c_2wf) + (e*c_2w3f);
                    di   = -(e*de1/(eta_2*t_i)) + (0.5*gammap2*c_i*sqrt(onemci2)*sewf);
                
                    onem5ci2   = 1.0 - (5.0*c_i_2);
                    onem5ci22  = onem5ci2*onem5ci2;
                    fmMpesf    = f - M + (e*s_f);
                    threem5ci2 = 3.0 - (5.0*c_i_2);
                    cewf       = (3*s_2w2f)+(3*e*s_2wf) + (e*s_2w3f);
                    L1         = (1.0/8.0)*gammap2*eta_3*de1c;
                    L21        = 11.0*(2.0 + 3.0*e_2)*c_i_2;
                    L22        = 40.0*(2.0 + 5.0*e_2)*(c_i_4/(onem5ci2));
                    L23        = 400.0*e_2*(c_i_6/onem5ci22);
                    L2         = -(gammap2/16.0)*(2.0 + e_2 - L21 - L22 - L23);
                    L31        = -6.0*onem5ci2*fmMpesf;
                    L32        = threem5ci2*cewf;
                    L3         = (1.0/4.0)*gammap2*(L31 + L32);
                    L41        = 11.0 + (80.0*c_i_2/onem5ci2) + (200.0*c_i_4/onem5ci22);
                    L4         = -(1.0/8.0)*gammap2*e_2*c_i*L41;
                    L5         = -0.5*gammap2*c_i*((6.0*fmMpesf) - cewf);
                    L          = M + w + O + L1 + L2 + L3 + L4 + L5;

                    saxr1 = (axr_2*eta_2) + axr + 1.0;
                    saxr2 = -(axr_2*eta_2) - axr + 1.0;
                    saxr3 = (axr_2*eta_2) + axr + (1.0/3.0);
                    edM1  = (1.0/8.0)*gammap2*e*eta_3*de1c;
                    edM21 = 2.0*threeci2m1*saxr1*s_f;
                    edM22 = 3.0*onemci2*(saxr2*s_2wf + saxr3*s_2w3f);
                    edM2  = -(1.0/4.0)*gammap2*eta_3*(edM21 + edM22);
                    edM   = edM1 + edM2;
                
                    dO  = L4 + L5;

                    //  Mean mean anomaly and eccentricity
                    d1 = (e + de)*s_M + edM*c_M;
                    d2 = (e + de)*c_M - edM*s_M;
                    Mm = atan2(d1,d2);
                    Mm = mod(Mm,PI2);
                    em = sqrt(d1*d1+d2*d2);

                    //  Mean right ascension of the ascending node and inclination
                    scix2 = s_ix2 + (c_ix2*di/2.0);
                    d3    = scix2*s_O + s_ix2*dO*c_O;
                    d4    = scix2*c_O - s_ix2*dO*s_O;
                    Om    = atan2(d3,d4);
                    Om    = mod(Om,PI2);
                    im    = 2.0*asin(sqrt(d3*d3+d4*d4));
                
                    //  Mean argument of perigee
                    wm = L - Mm - Om;
                    wm = mod(wm,PI2);
                
                    //  Mean eccentricity vector
                    exm = em*cos(wm);
                    eym = em*sin(wm);
                
                    //  Mean mean argument of latitude
                    um = wm + Mm;
                    um = mod(um,PI2);

                    //  Mean orbital elements vector
                    mean_elem << am,exm,eym,im,Om,um;

                    valid = true;

                    return(mean_elem);    
                    };
//------------------------------------------------------------------------------
// double EccAnomaly(double M, double e)
//------------------------------------------------------------------------------
/**
 * Computes the eccentric anomaly for elliptic orbits
 *
 * @param M   Mean anomaly [rad]
 * @param e   Eccentricity of the orbit
 *
 * @return Eccentric anomaly [rad]
 */
//------------------------------------------------------------------------------                
double EccAnomaly(double M,
                  double e,
                  const int maxit,
                  const double eps)
                  {
                  int i = 0;
                  double E, f;
            
                  // Starting value
            
                  //M = mod(M, PI2);
                  M = PI2*( M/PI2 - floor(M/PI2) );
                  if(e < 0.8) E = M;
                  else E = PI;
            
                  // Iteration
                  do
                    {
                    f = E - e*sin(E) - M;
                    E = E - f/( 1.0 - e*cos(E) );
                    i++;//++i;
                    if(i > maxit) break;
                    }
                  while( fabs(f) > eps );
                  
                  return(E);
                  };
//------------------------------------------------------------------------------
// VectorNd<2> lonlat2satcart(double lonS, double latS, double lonT, double latT, double i, double ki)
//------------------------------------------------------------------------------
/**
 * Compute coordinates of a point on Earth (target T) in spacecraft's ground-track coordinates given
 * longitude and latitudes of the spacecraft's sub-satellite point (S) and the point
 * 
 * @note The spacecraft's ground-track-frame is defined having the centre in the sub-satellite point,
 *       x-axis perpendicular to the ground-track and y-axis in the direction of the ground-track and
 *       pointing in the direction of motion of the su-satellite point, z in the direction of the line
 *       joining the centre of the Earth with the spacecraft
 * @note The Vincenty formula is used for the computation of the great circle distance between
 *       the sub-satellite point and the target
 *
 * @param lonS   Spacecraft longitude [deg]
 * @param latS   Spacecraft latitude [deg]
 * @param lonT   Target longitude [deg]
 * @param latT   Target latitude [deg]
 * @param inc    Orbit's inclination [deg]
 * @param ki     Daily recurrence frequency of the orbit
 *               (@see Handbook of Satellite Orbits by Michel Capderou)
 *
 * @return   xT, yT coordinates (spacecraft's ground-track-frame) [deg]
*/
//------------------------------------------------------------------------------                      
VectorNd<2> lonlat2satcart(double lonS,
                           double latS,
                           double lonT,
                           double latT,
                           double inc,
                           double ki)
                          {
                          VectorNd<2> xyT;
                          
                          lonS = lonS*DEG2RAD;
                          latS = latS*DEG2RAD;
                          lonT = lonT*DEG2RAD;
                          latT = latT*DEG2RAD;
                          inc = inc*DEG2RAD;
                          
                          double s1, s2, c1, c2, ci, sd, cd, y, x, sj, cj, xT, yT;
                          double delta_lon, delta_sigma, SMang, delta_theta, alpha;
                          
                          delta_lon = lonT - lonS;
                          
                          s1 = sin(latS);
                          c1 = cos(latS);
                          s2 = sin(latT);
                          c2 = cos(latT);
                          ci = cos(inc);
                          
                          sd = sin(delta_lon);
                          cd = cos(delta_lon);
                          
                          // Great circle distance between the sub-satellite point and the target
                          y = sqrt( (c2*sd)*(c2*sd) + (c1*s2 - s1*c2*cd)*(c1*s2 - s1*c2*cd) );
                          x = s1*s2 + c1*c2*cd;
                          
                          delta_sigma = atan2(y,x);
                          //delta_sigma = delta_sigma*RAD2DEG;
                          
                          //Angle between ground-track and local meridian
                          //y = fabs(ci) - (1.0/ki)*c1*c1;
                          //x = sqrt( c1*c1 - ci*ci );
                          //
                          //SMang = atan2(y,x);
                          //SMang = mod(SMang,PI2);
                          //double SMang1 = SMang;
                          
                          sj = fabs(ci)/c1;
                          double cj2 = 1.0 - sj*sj;
                          if(cj2 < 0.0) cj2 = 0.0;
                          cj = sqrt(cj2);
                          
                          y = sj - (1.0/ki)*c1;
                          x = cj;
                          
                          SMang = atan2(y,x);
                          SMang = mod(SMang,PI2);
                          //cout << "SMang1: " << SMang1*RAD2DEG << "    " << "SMang: " << SMang*RAD2DEG << endl;
                         
                          // Angle between ground-track and local meridian measured clocwise from the North direction
                          SMang = PI - SMang;
                          SMang = mod(SMang,PI2);
                          
                          //cout << "SMang: " << SMang*RAD2DEG << "   y: " << y << "   x: " << x << "    sj: " << sj << endl;
                          // Bearing going from sub-satellite point to the target on a great circle path (orthodrome).
                          // delta_theta is the angle between the great circle ST and the meridian passing through S and evaluated clockwise from the Noth direction 
                          y = sd*c2;
                          x = c1*s2 - s1*c2*cd;
                          
                          delta_theta = atan2(y,x); 
                          delta_theta = mod(delta_theta,PI2);
                          // Angle between great circle ST and ground-track measured counter-clockwise from -yT
                          alpha = SMang - delta_theta + PI2;
                          alpha = mod(alpha,PI2);
                          // Angle between great circle ST and ground-track measured counter-clockwise from xT
                          alpha = alpha - PI/2.0;
                          alpha = mod(alpha,PI2);
                          
                          // Cartesian coordinates of T in spacecraft's ground-track-frame
                          xT = (delta_sigma*cos(alpha))*RAD2DEG;
                          yT = (delta_sigma*sin(alpha))*RAD2DEG;
                          
                          xyT(0) = xT;
                          xyT(1) = yT;
                          
                          return(xyT);  
                          };
//------------------------------------------------------------------------------
// double GPS2ET(double GPSsecs)
//------------------------------------------------------------------------------
/**
 * Conversion of GPS seconds to ephemeris time (TDT). Ephemeris time seconds are seconds past from
 * the J2000 system reference epoch which is Greenwich noon on January 1, 2000 Barycentric Dynamical Time.
 * The currently-used standard epoch "J2000" is defined by international agreement to be equivalent to:
 * The Gregorian date January 1, 2000 at 12:00 TT (Terrestrial Time)
 * The Julian date 2451545.0 TT (Terrestrial Time)
 * January 1, 2000, 11:59:27.816 TAI (International Atomic Time)
 * January 1, 2000, 11:58:55.816 UTC (Coordinated Universal Time)
 *
 * @param GPSsecs   GPS seconds
 *
 * @return ET [s]
 */
//------------------------------------------------------------------------------                
double GPS2ET(double GPSsecs)
                  {
                  double ETsecs;
                  
                  ETsecs = GPSsecs - timescales::J2000_GPSSECS;
                  //ETsecs = GPSsecs - floor(timescales::J2000_GPSSECS);
                  //cout << fixed << "ETsecs: " << ETsecs << endl;
                  
                  return(ETsecs);
                  };
//------------------------------------------------------------------------------
// double GPS2TT(double GPSsecs)
//------------------------------------------------------------------------------
/**
 * Conversion of GPS seconds to terrestrial time (TT) seconds
 *
 * @param GPSsecs   GPS seconds
 *
 * @return TT [s]
 */
//------------------------------------------------------------------------------                
double GPS2TT(double GPSsecs)
                  {
                  double TTsecs;
                  
                  TTsecs = GPSsecs + timescales::TT_GPS;
                  
                  return(TTsecs);
                  };
//------------------------------------------------------------------------------
// double GPS2UTC(double GPSsecs, double leapsec)
//------------------------------------------------------------------------------
/**
 * Conversion of GPS seconds to UTC seconds
 *
 * @param GPSsecs   GPS seconds
 *
 * @return UTC [s]
 */
//------------------------------------------------------------------------------                
double GPS2UTC(double GPSsecs, double leapsec)
                  {
                  double UTCsecs;
                  
                  UTCsecs = GPSsecs + timescales::TAI_GPS - leapsec;
                  
                  return(UTCsecs);
                  };
//------------------------------------------------------------------------------
// VectorNd<2> GPS2GPSws(double GPSsecs)
//------------------------------------------------------------------------------
/**
 * Conversion of GPS seconds to GPS week and seconds of week
 *
 * @param GPSsecs   GPS seconds
 *
 * @return  2-dimensional vector containing GPS week count (week 0 starts at 1980/01/06.0 GPS time)
 *          and GPS seconds of week ([0s,604800s[)
 */
//------------------------------------------------------------------------------                
VectorNd<2> GPS2GPSws(double GPSsecs)
                  {
                  VectorNd<2> GPSws;
                  double GPSweek, GPSsecs_w;
                  
                  GPSweek = floor(GPSsecs/timescales::WEEK2SEC);
                  GPSsecs_w = GPSsecs - GPSweek*timescales::WEEK2SEC;
                  
                  GPSws(0) = GPSweek;
                  GPSws(1) = GPSsecs_w;
                  
                  return(GPSws);
                  };
//------------------------------------------------------------------------------
// double GPSws2GPSsecs(int week, double secs)
//------------------------------------------------------------------------------
/**
 * Conversion of GPS week and seconds of week to GPSsecs
 *
 * @param GPSws   2-dimensional vector containing GPS week count (week 0 starts at 1980/01/06.0 GPS time)
 *          and GPS seconds of week ([0s,604800s[)
 *
 * @return  GPS seconds
 */
//------------------------------------------------------------------------------                
double GPSws2GPSsecs(int week, double secs)
                  {
                  double GPSseconds = 0.0;
                  const double  CUC_TIME_RES  = 1E-6;
                  
                  double dw = floor(secs/604800.0);
                  int GPSweek = week + (int) dw;
                  double GPSsecs = secs - 604800.0*dw; // [0s,604800s[
                  
                  // CUC: Int secs since 6th Jan. 1980
                  unsigned int sec = (unsigned int)( GPSweek*604800.0 + floor(GPSsecs) );
                  if ( (unsigned int)(floor( (( GPSsecs - floor(GPSsecs) )/CUC_TIME_RES ) + 0.5 ) ) == 1e6 ) sec++;
                  
                  // CUC: Micro seconds.
                  unsigned int msec = (unsigned int)( floor( ( (GPSsecs - floor(GPSsecs))/CUC_TIME_RES) + 0.5 ) );
                  if(msec == 1e6) msec = 0;
                  
                  GPSseconds = (double)(sec + msec*1e-6);
                  
                  //cout << fixed << "GPSseconds = " << GPSseconds << endl;
                  //
                  //GPSsecs = 0.0;
                  //
                  //// CUC seconds
                  //int CUCsec = week*timescales::WEEK2SEC + floor(secs);
                  //if( floor( ( secs - floor(secs) )/1e-6 + 0.5 ) == 1e6) CUCsec++;
                  //
                  //// CUC microseconds
                  //int CUCmusec = floor( ( secs - floor(secs) )/1e-6 + 0.5 );
                  //if (CUCmusec == 1e6) CUCmusec = 0;
                  //
                  //GPSsecs = (double)(CUCsec + CUCmusec*1e-6);
                  //
                  //cout << fixed << "diff = " << GPSseconds - GPSsecs << endl;
                  //cout << "\n" << endl;
                  
                  return(GPSseconds);
                  };
                  
#ifdef USE_SPICE
//------------------------------------------------------------------------------
// Vec3d GPS2LST(double GPSsecs, double lon)
//------------------------------------------------------------------------------
/**
 * Conversion of GPS seconds to local solar time given the geodetic longitude
 *
 * @param GPSsecs   GPS seconds
 *
 * @return Local solar time vector LST[0] = hh, LST[1 = mm, LST[2] = ss
 */
//------------------------------------------------------------------------------                
Vec3d GPS2LST(double GPSsecs, double lon)
                  {
                  Vec3d LST = Vec3d::Zero();
                  double ETsecs;
                  int hh, mm, ss;
                  int body = 399; // SPICE ID-code for "EARTH";
                  const string lontype = "PLANETOCENTRIC";
                  int timlen = 51;
                  int ampmlen = 51;
                  SpiceChar ampm[51];
                  SpiceChar time[51];
                  
                  ETsecs = GPSsecs - timescales::J2000_GPSSECS;
                  
                  et2lst_c(ETsecs, body, lon, lontype.c_str(), timlen, ampmlen, &hh, &mm, &ss, time, ampm);
                  
                  LST << (double)hh, (double)mm, (double)ss;
                  
                  return(LST);
                  };
//------------------------------------------------------------------------------
// string GPS2UTCstr(double GPSsecs, const string format)
//------------------------------------------------------------------------------
/**
 * Conversion of GPS seconds to UTC string yyyy-mm-ddThh:mm:ss.sss
 *
 * @param GPSsecs       GPS seconds
 * @param format        "C"      "1986 APR 12 16:31:09.814"
                        "D"      "1986-102 // 16:31:12.814"
                        "J"      "JD 2446533.18834276"
                        "ISOC"   "1987-04-12T16:31:12.814"
                        "ISOD"   "1987-102T16:31:12.814"
 *
 * @return UTC date string in ISOC format (yyyy-mm-ddThh:mm:ss.sss)
 */
//------------------------------------------------------------------------------                
string GPS2UTCstr(double GPSsecs,
                  const string format)
                  {
                  //Vector6d UTCdate = Vector6d::Zero();
                  double sec_J2000;
                  SpiceChar UTCdate_char[35];
                  //ConstSpiceChar* format;
                  
                  sec_J2000 = GPS2ET(GPSsecs);
                  
                  //format = "DD-MM-YYYY HR:MN:SC";
                  //timout_c(sec_J2000, format, 35, UTCdate_char);
                  // ISOC format is yyyy-mm-ddThh:mm:ss.sss
                  et2utc_c(sec_J2000, format.c_str(), 3, 35, UTCdate_char);
                  
                  string UTCdate_string(UTCdate_char);
                  
                  return(UTCdate_string);
                  };
//------------------------------------------------------------------------------
// Vector6d GPS2UTCdate(double GPSsecs, const string format)
//------------------------------------------------------------------------------
/**
 * Conversion of GPS seconds to UTC yyyy mm dd hh mm ss.sss
 *
 * @param GPSsecs       GPS seconds
 * @param format        "ISOC"
                        "ISOD"
 *
 * @return Vector [yyyy mm dd hh mm ss.sss] if format is "ISOC" and [yyyy doy hh mm ss.sss 0] if format is "ISOD"
 */
//------------------------------------------------------------------------------                
Vector6d GPS2UTCdate(double GPSsecs,
                     const string format)
                  {     
                  Vector6d UTCdate = Vector6d::Zero();
                  string UTCdate_string;
                  string year, doy, month, day, h, m, s;
                  
                  UTCdate_string = GPS2UTCstr(GPSsecs, format);
                  
                  if( format.compare("ISOC") == 0 )
                        {
                        // "ISOC"   "1987-04-12T16:31:12.814"
                        year = UTCdate_string.substr(0,4);
                        month = UTCdate_string.substr(5,2);
                        day = UTCdate_string.substr(8,2);
                        h = UTCdate_string.substr(11,2);
                        m = UTCdate_string.substr(14,2);
                        s = UTCdate_string.substr(17,6);
                        
                        UTCdate(0) = atof(year.c_str());
                        UTCdate(1) = atof(month.c_str());
                        UTCdate(2) = atof(day.c_str());
                        UTCdate(3) = atof(h.c_str());
                        UTCdate(4) = atof(m.c_str());
                        UTCdate(5) = atof(s.c_str());
                        }
                  else if( format.compare("ISOD") == 0 )
                        {
                        // "ISOD"   "1987-102T16:31:12.814"
                        year = UTCdate_string.substr(0,4);
                        doy = UTCdate_string.substr(5,3);
                        h = UTCdate_string.substr(9,2);
                        m = UTCdate_string.substr(12,2);
                        s = UTCdate_string.substr(15,6);
                        
                        UTCdate(0) = atof(year.c_str());
                        UTCdate(1) = atof(doy.c_str());
                        UTCdate(2) = atof(h.c_str());
                        UTCdate(3) = atof(m.c_str());
                        UTCdate(4) = atof(s.c_str());
                        }
                  
                  return(UTCdate);
                  };
//------------------------------------------------------------------------------
// double UTCdate2GPSsecs(Vector6d& UTCdate, const string format)
//------------------------------------------------------------------------------
/**
 * Conversion of UTC yyyy mm dd hh mm ss.sss to GPS seconds
 *
 * @param UTCdate   Vector of integers [yyyy mm dd hh mm ss.sss]. The time is given on a 24-hour clock
 *
 * @return GPS seconds 
 */
//------------------------------------------------------------------------------                
double UTCdate2GPSsecs(Vector6d& UTCdate)
                      {
                      double ETsecs, GPSsecs;
                      string UTCdate_string;
                      string year, month, day, h, m, s;
                      
                      year = to_string( (int)UTCdate(0) );
                      month = to_string( (int)UTCdate(1) );
                      day = to_string( (int)UTCdate(2) );
                      h = to_string( (int)UTCdate(3) );
                      m = to_string( (int)UTCdate(4) );
                      s = to_string( UTCdate(5) );
                      
                      UTCdate_string = year + "/" + month + "/" + day + " " + h + ":" + m + ":" + s;
                      
                      str2et_c(UTCdate_string.c_str(), &ETsecs);  // In SPICE function str2et_c default interpretation of the input string is the time of day to be a time on a 24-hour clock in the UTC time system
                                                                  // An identical result is obtained using utc2et_c(UTCdate_string.c_str(), &ETsecs);
                      
                      GPSsecs = round(ETsecs + timescales::J2000_GPSSECS);
                      
                      return(GPSsecs);
                      };            
//------------------------------------------------------------------------------
// double GPS2DOY(double GPSsecs)
//------------------------------------------------------------------------------
/**
 * Conversion of GPS seconds to day of year
 *
 * @param  GPS seconds 
 *
 * @return DOY
 */
//------------------------------------------------------------------------------                
double GPS2DOY(double GPSsecs)
            {
            double sec_J2000, doy;
            SpiceChar doy_chr[35];
            string doy_str;
            ConstSpiceChar* format;
            format = "DOY";
            
            sec_J2000 = GPS2ET(GPSsecs);
            timout_c( sec_J2000, format, 35, doy_chr );
            
            doy_str = string(doy_chr);
            
            doy = stod(doy_str);
            
            return(doy);
            };
//------------------------------------------------------------------------------
// double leapsec(double GPSsecs)
//------------------------------------------------------------------------------
/**
 * Gives leap seconds at epoch GPS seconds
 *
 * @param GPSsecs   GPS seconds
 *
 * @return Leap seconds [s]
 */
//------------------------------------------------------------------------------                
double leapsec(double GPSsecs)
                  {
                  double leapsecs;
                  double time = GPS2ET(GPSsecs);
                  double delta;
                  
                  deltet_c(time,"ET",&delta);
                  leapsecs = round(delta - 32.184);
                  
                  return(leapsecs);
                  };
//------------------------------------------------------------------------------
// double GPS2Unix(double GPSsecs)
//------------------------------------------------------------------------------
/**
 * Conversion of GPS seconds to Unix time
 *
 * @param GPSsecs   GPS seconds
 *
 * @return Unix time [s]
 */
//------------------------------------------------------------------------------                
double GPS2Unix(double GPSsecs)
                  {
                  double Unixsecs;
                  
                  double leapsecs = leapsec(GPSsecs);
                  
                  Unixsecs = GPSsecs + 315964819 - leapsecs;
                  
                  return(Unixsecs);
                  };
//------------------------------------------------------------------------------
// double GPS2JD(double GPSsecs)
//------------------------------------------------------------------------------
/**
 * Conversion of GPS seconds to Julian date
 *
 * @param GPSsecs   GPS seconds
 *
 * @return Julian date [d]
 */
//------------------------------------------------------------------------------                
double GPS2JD(double GPSsecs)
                  {
                  double JD = 0.0;
                  double time = GPS2ET(GPSsecs);
                  
                  int lenout = 35;
                  char JDstring[35];
                  
                  et2utc_c(time, "J", 6, lenout, JDstring);
                  
                  // Convert char array into string
                  string JDstr(JDstring);
                  // Delete "JD" from julian date string
                  size_t found = JDstr.find("JD");
                  JDstr.erase(found,2);
                  // Convert Julian date string into double
                  JD = stod(JDstr);
                  //cout << "JD: " << JDstr << endl;
                  //cout << "JD: " << JD << endl;
                  
                  return(JD);
                  };
//------------------------------------------------------------------------------
// GPS2MJD(double GPSsecs)
//------------------------------------------------------------------------------
/**
 * Conversion of GPS seconds to modified Julian date
 *
 * @param GPSsecs   GPS seconds
 *
 * @return Modified Julian date [d]
 */
//------------------------------------------------------------------------------                
double GPS2MJD(double GPSsecs)
                  {
                  double MJD = 0.0;
                  double JD = 0.0;
                  
                  JD = GPS2JD(GPSsecs);
                  MJD = JD - 2400000.5;
                  
                  return(MJD);
                  };
#endif
//------------------------------------------------------------------------------
// double GPSws2MJD(double GPSweek, double GPSsecs_w)
//------------------------------------------------------------------------------
/**
 * GPS week and seconds of week to modified Julian date
 *
 * @param   GPS week count (week 0 starts at 1980-01-06 GPS time)
 * @param   GPS seconds of week ([0s,604800s[)
 *
 * @return Modified Julian date [d]
 */
//------------------------------------------------------------------------------                
double GPSws2MJD(double GPSweek, double GPSsecs_w)
                  {
                  double MJD;
                  
                  MJD = timescales::MJD_GPST0 + 7.0*GPSweek + GPSsecs_w/timescales::SEC2DAY;
                  
                  return(MJD);
                  };
//------------------------------------------------------------------------------
// double GPSws2JD(double GPSweek, double GPSsecs_w)
//------------------------------------------------------------------------------
/**
 * GPS week and seconds of week to modified Julian date
 *
 * @param   GPS week count (week 0 starts at 1980-01-06 GPS time)
 * @param   GPS seconds of week ([0s,604800s[)
 *
 * @return Julian date [d]
 */
//------------------------------------------------------------------------------                
double GPSws2JD(double GPSweek, double GPSsecs_w)
                  {
                  double JD = 0.0;
                  double MJD;
                  
                  MJD = timescales::MJD_GPST0 + 7.0*GPSweek + GPSsecs_w/timescales::SEC2DAY;
                  
                  JD = MJD + 2400000.5;
                  
                  return(JD);
                  };    
//------------------------------------------------------------------------------
// Vector6d MJD2UTCdate(double MJD)
//------------------------------------------------------------------------------
/**
 * Conversion of modified Julian date to UTC yyyy mm dd hh mm ss.sss
 *
 * @see Montenbruck, O., and Gill, E.,“Satellite Orbits - Model, Methods and Applications”,
 *      Springer Verlag, Heidelberg, Germany, 2000, ISBN:3-540-67280-X
 *
 * @param MJD       Modified Julian date
 *
 * @return Vector [yyyy mm dd hh mm ss.sss] if format is "ISOC" and [yyyy doy hh mm ss.sss 0] if format is "ISOD"
 */
//------------------------------------------------------------------------------                
Vector6d MJD2UTCdate(double MJD)
                              {     
                              Vector6d UTCdate = Vector6d::Zero();
                              
                              double year = 0.0, month = 0.0, day = 0.0, hour = 0.0, min = 0.0, sec = 0.0;
                              double hours, h_diff;      
                              long a, b, c, d, e, f, q;
              
                              // Convert Julian days number to calendar date
                              a = long(MJD) + 2400001.0;
                              q = MJD - long(MJD);
                            
                              if(a < 2299161) // Julian calendar
                                {  
                                b = 0;
                                c = a + 1524;
                                }
                              else // Gregorian calendar
                                {                
                                b = long( (a - 1867216.25)/36524.25 );
                                c = a + b - long(b/4) + 1525;
                                }
                            
                              d = long( (c - 121.1)/365.25 );
                              e = long(365.25*d);//e = 365*d + d/4;
                              f = long( (c - e)/30.6001 );
                            
                              day = c - e - long(30.6001*f) + q;//day = c - e - int(30.6001*f);
                              
                              month = f - 1 - 12*long(f/14);
                              
                              year = d - 4715 - long( (7 + month)/10 );
                            
                              hours = 24.0*( MJD - floor(MJD) );
                              hour = int(hours);
                              
                              h_diff = (hours - hour)*60.0;
                              
                              min = int(h_diff);
                              
                              sec = (h_diff - min)*60.0;
                              
                              UTCdate << year, month, day, hour, min, sec;
                              
                              return(UTCdate);
                              };
//------------------------------------------------------------------------------
// double UTCdate2MJD(Vector6d& UTCdate)
//------------------------------------------------------------------------------
/**
 * Conversion of UTC yyyy mm dd hh mm ss.sss to modified Julian date
 *
 * @see Montenbruck, O., and Gill, E.,“Satellite Orbits - Model, Methods and Applications”,
 *      Springer Verlag, Heidelberg, Germany, 2000, ISBN:3-540-67280-X
 *
 * @param UTCdate   Vector of integers [yyyy mm dd hh mm ss.sss]. The time is given on a 24-hour clock
 *
 * @return Modified Julian date
 */
//------------------------------------------------------------------------------                
double UTCdate2MJD(Vector6d& UTCdate)
                  {     
                  double MJD = 0.0;
                  
                  int year, month, day, hour, min, B;
                  double sec;
                  long MJD_Midnight;
                  double Fraction_of_day;
                  
                  year = (int)UTCdate(0);
                  month = (int)UTCdate(1);
                  day = (int)UTCdate(2);
                  hour = (int)UTCdate(3);
                  min = (int)UTCdate(4);
                  sec = UTCdate(5);

                  if(month <= 2)
                    {
                    month += 12;
                    --year;
                    }
                  
                  if( (10000L*year + 100L*month + day) <= 15821004L ) B = -2 + int( (year + 4716)/4 ) - 1179; // Julian calendar 
                  else B = int(year/400) - int(year/100) + int(year/4); // Gregorian calendar 
                    
                  MJD_Midnight = 365L*year - 679004L + B + int( 30.6001*(month + 1) ) + day;
                  Fraction_of_day = ( hour + min/60.0 + sec/3600.0 )/24.0; 
                
                  MJD = MJD_Midnight + Fraction_of_day;
                              
                  return(MJD);
                  };              
//------------------------------------------------------------------------------
// double mod(double x, double y)
//------------------------------------------------------------------------------
/**
 * Modulo function
 *
 * @param x   Value to be normalized
 * @param y   Reference value for the normalization
 *
 * @return Normalized value of input
 */
//------------------------------------------------------------------------------                
double mod(double x, double y)
        {
        double z;
        z = x - y*floor(x/y);
        
        return(z);
        };
//------------------------------------------------------------------------------------------------
// Vec3d vectcorr_cone(Vec3d& v_in, double alpha, double theta)
//------------------------------------------------------------------------------------------------
/**
 * Correct a vector into another vector which is on the surface of a cone whose axis is the
 * initial vector and with aperture alpha. The final point vector is positioned on the circular
 * of the cone. The specific position is determined parameter theta of the circle parametric
 * representation in a 3D space.
 * A typical application of this function is the introduction of the attitude error into the
 * dv of an orbital maneuver.
 *
 * @param v_in       3-dimensional input vector
 * @param alpha      Aperture of the conical error [rad]
 * @param theta      Parameter [rad]. 0-0 <= theta <= 2.0*PI
 *
 * @return 3-dimensional corrected vector
 */
//------------------------------------------------------------------------------------------------
Vec3d vectcorr_cone(Vec3d& v_in,
               double alpha,
               double theta)
                  {
                  Vec3d v_out = Vec3d::Zero();
                  Vec3d u = Vec3d::Zero();
                  Vec3d c = Vec3d::Zero();
                  Vec3d P = Vec3d::Zero();
                  Vec3d PC = Vec3d::Zero();
                  Vec3d i1 = Vec3d::Zero();
                  Vec3d j1 = Vec3d::Zero();
                  
                  double vx, vy, vz, vxy, Pxy, v, rho;
                  double sinalpha, cosalpha, singamma, cosgamma, sinbeta, cosbeta, sinphi, cosphi;
                  
                  sinalpha = sin(alpha);
                  cosalpha = cos(alpha);
                  
                  vx = v_in(0);
                  vy = v_in(1);
                  vz = v_in(2);
                  
                  v = v_in.norm();
                  
                  u = v_in/v;
                  c = v*cosalpha*u; // Centre of the spatial circle
                  rho = v*sinalpha; // Radius of spatial circle
                  
                  if(v == vx) // v_in aligned with x-axis
                        {
                        i1 << 0.0, 1.0, 0.0;      
                        j1 << 0.0, 0.0, 1.0;
                        }
                  else if(v == vy)
                        {
                        i1 << 0.0, 0.0, 1.0;      
                        j1 << 1.0, 0.0, 0.0;
                        }
                  else if(v == vz)
                        {
                        i1 << 1.0, 0.0, 0.0;      
                        j1 << 0.0, 1.0, 0.0;
                        }
                  else
                        {
                        vxy = sqrt(vx*vx + vy*vy);
      
                        singamma = vz/v;
                        cosgamma = vxy/v;
                        //cosgamma = sqrt(1.0 - singamma*singamma);
                        
                        cosphi = vx/vxy;
                        sinphi = vy/vxy;
                        //sinphi = sqrt(1.0 - cosphi*cosphi);
                        
                        sinbeta = singamma*cosalpha - cosgamma*sinalpha;
                        cosbeta = cosgamma*cosalpha + singamma*sinalpha;
                        
                        Pxy = v*cosbeta;
                        
                        P(0) = Pxy*cosphi;
                        P(1) = Pxy*sinphi;
                        P(2) = v*sinbeta;
                        
                        PC = P - c;
                        
                        i1 = PC/PC.norm();
                        j1 = u.cross(i1);
                        }
                  
                  v_out = c + rho*cos(theta)*i1 + rho*sin(theta)*j1;
                  
                  return(v_out);
                  };
//------------------------------------------------------------------------------------------------
// Vec3d vectcorr_norm(Vec3d& v_in, double newnorm)
//------------------------------------------------------------------------------------------------
/**
 * Correct a vector into another vector which has same direction and pointing but different norm
 * A typical application of this function is the introduction of the performance error of a propulsion
 * system inexecuting a commanded dv of an orbital maneuver.
 *
 * @param v_in       3-dimensional input vector
 * @param newnorm    New vector norm
 *
 * @return 3-dimensional corrected vector
 */
//------------------------------------------------------------------------------------------------
Vec3d vectcorr_norm(Vec3d& v_in,
                    double newnorm)
                    {
                    Vec3d v_out = Vec3d::Zero();
                    
                    double vx, vy, vz, v;
                    
                    vx = v_in(0);
                    vy = v_in(1);
                    vz = v_in(2);
                    
                    v = v_in.norm();
                    
                    v_out(0) = newnorm*vx/v;
                    v_out(1) = newnorm*vy/v;
                    v_out(2) = newnorm*vz/v;
                    
                    return(v_out);
                    };
//------------------------------------------------------------------------------
// Vec4d& toQ(double phi, double theta, double psi)
//------------------------------------------------------------------------------
/**
 * Compute quaternion from Euler angles
 *
 * @param phi, theta, psi   Input Euler angles
 *
 * @return Quaternion for Euler rotation 3-2-1
 */
//------------------------------------------------------------------------------  
//Vec4d& toQ(double phi, double theta, double psi)
//                    {
//                    Vec4d q;
//                    
//                    return(q);
//                    };
        