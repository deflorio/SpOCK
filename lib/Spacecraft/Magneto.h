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

#ifndef MAGNETO_H
#define MAGNETO_H

#include <random>

#include <Subsystem.h>
#include <MagneticField.h>
#include <Transformations.h>
#include <VarTypes.h>
// External libraries
#include <Eigen/Core>

using namespace math;
using namespace std;
using namespace subsystem;
using namespace magnetic;

namespace magneto
   {
    //------------------------------------------------------------------------------
    //! Class MAGNETO
    //------------------------------------------------------------------------------
    /*!
       Class derived from SUBSYS which implements the magnetometer and magnetotorquer
       models
     */
    //------------------------------------------------------------------------------ 
    class MAGNETO : public SUBSYS
        {
        public:
        //! Constructor.
        /*!
            Using class SUBSYS constructor
          */
        using SUBSYS::SUBSYS;
        //! Destructor.
        ~MAGNETO();
        //------------------------------------------------------------------------------
        //
        // Class methods specification
        //
        //------------------------------------------------------------------------------
        // Evaluate the operability status (subsystem on or off)
        VectorXd Output(double epoch, const Ref<const VectorXd>& currentstate, const Ref<const VectorXd>& auxstate);
        // Store in clas members constant parameters contained in base class' member vector ConstPrm
        void Init();
        
        public:
        /** Vector of magnetic torquer ontimes commanded by the attitude controller.*/
        Vec3d ontime;
        
        private:
        // Evaluate the operability status (subsystem on or off)
        void status(double epoch, const Ref<const VectorXd>& currentstate, const Ref<const VectorXd>& auxstate) ;
        
        private:
        /** Object of class MAGFIELD type.*/
        MAGFIELD MagneticField;
        /** Magnetic field model path.*/
        string magneticfield;
        /** Conversion factor from ontime to dipole produced by a magnetic torquer ( m = ontime2dipole*ontime ).*/
        double ontime2dipole;
        };

   }; // End of namespace magnetic

#endif
