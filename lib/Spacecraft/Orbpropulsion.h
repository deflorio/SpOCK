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

#ifndef ORBPROPULSION_H
#define ORBPROPULSION_H

#include <random>

#include <Subsystem.h>
#include <Transformations.h>
#include <VarTypes.h>
// External libraries
#include <Eigen/Core>

using namespace math;
using namespace std;
using namespace subsystem;

namespace orbpropulsion
   {
    //------------------------------------------------------------------------------
    //! Class ORBPROPULSION
    //------------------------------------------------------------------------------
    /*!
       Class derived from SUBSYS which implements the magnetometer and magnetotorquer
       models
     */
    //------------------------------------------------------------------------------ 
    class ORBPROPULSION : public SUBSYS
        {
        public:
        //! Constructor.
        /*!
            Using class SUBSYS constructor
          */
        using SUBSYS::SUBSYS;
        //! Destructor.
        ~ORBPROPULSION();
        //------------------------------------------------------------------------------
        //
        // Class methods specification
        //
        //------------------------------------------------------------------------------
        // Evaluate the operability status (subsystem on or off)
        VectorXd Output(double epoch, const Ref<const VectorXd>& currentstate, const Ref<const VectorXd>& auxstate);
        // Store in clas members constant parameters contained in base class' member vector ConstPrm
        void Init();
        // Convert thrust parameters of  the propulsion system to dv parameters given a spacecraft mass
        void thrust2dv(double SC_mass);
        // Convert dvs to ontimes given the propulsion system parameters
        Vec3d dv2ontimes(Vec3d& dv);
        // Convert an orbit impulsive maneuvers vector to a continuous maneuvers vector given the propulsion system parameters
        void impman2contman(vector<maneuver> &impman, int SIM_STEP);
        
        public:
        /** dv vector commanded to the propulsion system*/
        Vec3d dv_CMD;
        /** Reference system of the commanded maneuver*/
        string man_sys;
        
        private:
        // Evaluate the operability status (subsystem on or off)
        void status(double epoch, const Ref<const VectorXd>& currentstate, const Ref<const VectorXd>& auxstate) ;
        
        private:
        /** Propulsion system nominal thrust*/
        double thrust;
        /** Propulsion system nominal dv when thrusting for 1 s*/
        double dv_thr;
        /** Propulsion system thrust resolution*/
        double thrust_res;
        /** Propulsion system dv resolution*/
        double dv_res;
        /** Maximum thrust generated when firing continuously after the exaustion of the propellent*/
        double thrust_MAX;
        /** Propulsion system dv budget*/
        double dv_MAX;
        /** Normal distribution for positioning dv on the conference representing all the possible real dvs given the attitude accuracy.*/
        normal_distribution<double> theta_distr;
        /** Random numbers generator.*/
        mt19937 theta_generator;
        /** Propulsion system dv parameters comuted from thrust parameters*/
        };

   }; // End of namespace magnetic

#endif
