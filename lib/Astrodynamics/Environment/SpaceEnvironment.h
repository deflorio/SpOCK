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

#ifndef SPACEENV_H
#define SPACEENV_H

#include <optional>
#include <Transformations.h>

// External libraries: Eigen
#include <Eigen/Core>

#include <VarTypes.h>

using namespace math;
using namespace std;
using namespace Eigen;

//
// Namespace spaceenvironment
//

namespace spaceenvironment
   {
    //------------------------------------------------------------------------------
    //! Class SPACEENV
    //------------------------------------------------------------------------------
    /*!
       Base class for space environment modeling
     */
    //------------------------------------------------------------------------------ 
    class SPACEENV
         {
         public:
         //! Constructor.
         /*!
             Using class SPACEENV constructor
           */
         SPACEENV();
         //! Destructor.
         virtual ~SPACEENV();
         // Set model's file path
         void setfilespath(const string model_filename);
         // Set model file's folder path and name
         void setfilespath(const string model_filename, const string model_filepath);
         //------------------------------------------------------------------------------
         // Abstract method Vec3d field_vec(double time, const Ref<const VectorXd>& orbstate, const Ref<const VectorXd>& auxparam = MatrixXd::Zero(4,1) )
         //------------------------------------------------------------------------------
         /**
           * Compute 3D model vector (field, acceleration, force, etc.)
           *
           * @param time       Epoch
           * @param orbstate   Orbital state vector
           * @param auxparam   Optional auxiliary parameters vector
           *
           * @return 3D model vector (field, acceleration, force, etc.)
           */
         //------------------------------------------------------------------------------  
         virtual Vec3d field_vec(double time, const Ref<const VectorXd>& orbstate) = 0;
         // Set the reference frame in which the field's vector is computed
         void SetReferenceFrame(const string ref_sys);
         // Get model coefficients
         virtual void getmodel_coeff() = 0;
                       
         public:
         /** Maximum limit of Generic auxiliary coefficient ravity field max degree.*/
         static int n_max;
         /** Area required to compute forces depending on surface.*/
         double Area;
         /** Generic force coefficient (e.g. CD = drag coefficient, etc.).*/
         double CF;
         /** 3-dimensional vector containing Earth orientation parameters. */
         static Vec3d eop;
         /** Current leap second. */
         static double tai_utc;
         /** Transformation matrix from ICRF to ITRF. */
         Mat3x3d T_ICRF2ITRF;
         
         //private:
         //static void DynModel(const state_type x , state_type &dxdt , const double t); // Definition of dynamic model
     
         protected:
         /** Model name.*/
         string modelname;
         /** Model's file path or model file's folder path.*/
         string modelfilepath;
         /** Reference frame of output of field_vec.*/
         string refsys;
         /** Generic auxiliary coefficient.*/
         double rho1;
         /** Generic auxiliary coefficient.*/
         double rho2;
         /** Generic auxiliary coefficient.*/
         double rho3;
         };

   }; // End of namespace spaceenvironment

#endif
