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

#ifndef GRAV_H
#define GRAV_H

#include <SpaceEnvironment.h>
#include <Solarsys.h>
#include <Transformations.h>
#include <VarTypes.h>
// External libraries
#include <Eigen/Core>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/factorials.hpp>
//#include <GeographicLib/include/MagneticModel.hpp>
extern "C"
      {
       #include "extlib/cspice/include/SpiceUsr.h"
      }

using namespace spaceenvironment;
using namespace math;
using namespace std;
//using namespace GeographicLib;
using namespace solarsystem;

namespace gravity
   {
    //------------------------------------------------------------------------------
    //! Class GRAV
    //------------------------------------------------------------------------------
    /*!
       Class derived from SPACEENV which implements the solar radiation pressure model
     */
    //------------------------------------------------------------------------------ 
    class GRAV : public SPACEENV
        {
        public:
        //! Constructor.
        /*!
            Using class SPACEENV constructor
          */
        using SPACEENV::SPACEENV;
        //! Destructor.
        ~GRAV();
        //------------------------------------------------------------------------------
        //
        // Class methods specification
        //
        //------------------------------------------------------------------------------
        // Compute the central body gravitational acceleration vector
        Vec3d field_vec(double time, const Ref<const VectorXd>& orbstate);
        // Get gravity field model coefficients
        void getmodel_coeff();
        
        private:
		  // Normalized associated Legendre functions (c++17 function)
        double Legendre_norm(unsigned int n, unsigned int m, double x);
        // First derivative of normalized associated Legendre functions (c++17 function)
        double d_Legendre_norm(unsigned int n, unsigned int m, double x);
        // Normalized associated Legendre functions and first derivatives (forward columns recursive method)
        static void Legendre_norm_FC(double u, double w);
        
        //MatrixXd Legendre_norm_FR(double u, double w);
        //MatrixXd d_Legendre_norm_FR(const Ref<const MatrixXd>& L_FR, double u, double w);
        
        private:
        // Matrix of normalized C gravity field model coefficients
        static MatnMAXxnMAXd C;//Eigen::MatrixXd C;
        // Matrix of normalized S gravity field model coefficients
        static MatnMAXxnMAXd S;//Eigen::MatrixXd S;
        // Matrix of normalized sigmaC gravity field model coefficients
        static MatnMAXxnMAXd sigmaC;//Eigen::MatrixXd sigmaC;
        // Matrix of normalized sigmaS gravity field model coefficients
        static MatnMAXxnMAXd sigmaS;//Eigen::MatrixXd sigmaS;
        // Matrix of normalized associated Legendre functions
        static MatnMAXxnMAXd P;
        // Matrix of first derivative of normalized associated Legendre functions
        static MatnMAXxnMAXd Pd;
        
        static MatnMAXxnMAXd A;
        
        static MatnMAXxnMAXd B;
        
        static MatnMAXxnMAXd F;
        
        static MatnMAXxnMAXd CS;
        
        static MatnMAXxnMAXd CS1;
        
        static VectornMAXd Sigma1;
        
        static VectornMAXd sin_mlambda;
        
        static VectornMAXd cos_mlambda;
        // Gravitational coefficient Earth (from gravity model file)
        static double mu;
        // Earth's equatorial radius (from gravity model file)
        static double R;
        
        public:
        // Epoch at which the gravity field is evaluated. This variable is used only when gravity field models containing time variable parameters are read (e.g. EIGEN-6S)
        double grav_epoch;
        };
   }; // End of namespace magnetic
   
                              

#endif
