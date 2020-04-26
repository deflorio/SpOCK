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

#ifndef ATMO_H
#define ATMO_H

#include <SpaceEnvironment.h>
#include <Solarsys.h>
#include <Transformations.h>
#include <VarTypes.h>
// External libraries
#include <Eigen/Core>

extern "C"
      {
       #include "extlib/cspice/include/SpiceUsr.h"
      }

using namespace spaceenvironment;
using namespace math;
using namespace std;
using namespace solarsystem;

//typedef Eigen::Matrix<double,1,32> SpW_row;

namespace atmosphere
   {
   //------------------------------------------------------------------------------
   //! Class ATMO
   //------------------------------------------------------------------------------
   /*!
      Class derived from SPACEENV which implements the solar radiation pressure model
    */
   //------------------------------------------------------------------------------ 
   class ATMO : public SPACEENV
      {
      public:
      //! Constructor.
      /*!
          Using class SPACEENV constructor
        */
      using SPACEENV::SPACEENV;
      //! Destructor.
      ~ATMO();
      //------------------------------------------------------------------------------
      //
      // Class methods specification
      //
      //------------------------------------------------------------------------------
      // Compute atmospheric density with a chosen model
      void AtmosphericDensity(double time, const Ref<const VectorXd>& orbstate);
      // Set up surface parameters and orientation required to compute the solar radiation pressure
      void SetSurfaceParameters(SC::Face face_prms, Vec4d& q_attitude);
      // Compute the force applied on a surface by the solar radiation pressure using a multiple-surfaces model
      Vec3d field_vec(double time, const Ref<const VectorXd>& orbstate);
      // Get surface optical coefficients
      void getmodel_coeff();
      
      public:
      
      /** Matrix containing space weather indices. */
      MatrixXf SpaceWeather_idx;
      /** Vector containing nrlmsise00 model's atmospheric gases number densities and total mass density (nrlmsise00_rho_vec[5]). */
      float nrlmsise00_rho_vec[9];
      /** Vector containing JB2008 atmospheric model exospheric temperature and temperature at altitude. */
      float nrlmsise00_temperatures[2];
      /** Vector containing JB2008 atmospheric model exospheric temperature and temperature at altitude. */
      double jb2008_temperatures[2];
      /** Initial epoch to read space weather indices file. */
      double init_epoch;
      /** Atmospheric density.*/ 
      double rho_atm;
      /** Atmospheric drag model used.*/
      string Drag_Model;
      /** Propagation duration.*/
      int simduration;
      
      private:
      /** Object of class SOLSYS type.*/
      SOLSYS Solar;
      //SOLSYS* Solar_ptr = &Solar;
      /** Row of space weather indices matrix currently in use. */
      Eigen::Matrix<float,1,27> Ap_4days;
      /** Unit vector normal to a surface (ECI).*/
      Vec3d n;
      /** Indices for JB2008 atmospheric model.*/
      static double f10, s10, y10, f10b, s10b, y10b, m10, m10b;
      static int dstdtc;
      /** Space weather indices files for JB2008 atmospheric model.*/
      static char c_dtc_file[200];
      static char c_dst_file[200];
      static char c_solfsmy_file[200];
      /** Day of year at propagation start epoch.*/
      static double init_doy;
      /** Flag to lock space weather file lines to propagation time.*/
      static bool idx_locked;
      /** Line number of space weather indices matrix.*/
      static int swind;
      };

   }; // End of namespace magnetic

#endif
