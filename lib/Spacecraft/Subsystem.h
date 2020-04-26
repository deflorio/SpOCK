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

#ifndef SUBSYS_H
#define SUBSYS_H

#include <map>
#include <random>

#include <VarTypes.h>
// External libraries
#include <Eigen/Core>

using namespace math;
using namespace std;
using namespace Eigen;
using namespace SC; // SYS_params in VarTypes

namespace subsystem
   {
    //------------------------------------------------------------------------------
    //! Class SUBSYS
    //------------------------------------------------------------------------------
    /*!
       Base class for spacecraft subsystems
     */
    //------------------------------------------------------------------------------ 
    class SUBSYS
        {
        public:
        //! Constructors.
        SUBSYS();
        SUBSYS(SYS_params& param);
        SUBSYS(SYS_params& param, const Ref<const VectorXd>& Telecommands);
        //! Destructor.
        virtual ~SUBSYS();
        //------------------------------------------------------------------------------
        //
        // Class methods specification
        //
        //------------------------------------------------------------------------------
        // Parameters allocation
        void Setup(SYS_params& param);
        //------------------------------------------------------------------------------
        // Abstract method VectorXd TM(double epoch, const Ref<const VectorXd>& currentstate, const Ref<const VectorXd>& auxstate)
        //------------------------------------------------------------------------------
        /**
          * Compute the subsystem output (TM data) given the input data
          *
          * @param epoch          Epoch of the input states
          * @param currentstate   State vector (e.g. attitude or orbit dynamic model state vector)
          * @param auxstate       Auxiliary state vector (e.g. attitude or orbit dynamic model state vector)
          *
          * @return Sensor readings
          */
        //------------------------------------------------------------------------------  
        virtual VectorXd Output(double epoch, const Ref<const VectorXd>& currentstate, const Ref<const VectorXd>& auxstate) = 0;
        //------------------------------------------------------------------------------
        // Abstract method bool status(double epoch, const Ref<const VectorXd>& currentstate, const Ref<const VectorXd>& auxstate)
        //------------------------------------------------------------------------------
        /**
          * Evaluate the operability status (subsystem on or off)
          *
          * @param epoch          Epoch of the input states
          * @param currentstate   State vector (e.g. attitude or orbit dynamic model state vector)
          * @param auxstate       Auxiliary state vector (e.g. attitude or orbit dynamic model state vector)
          */
        //------------------------------------------------------------------------------  
        virtual void status(double epoch, const Ref<const VectorXd>& currentstate, const Ref<const VectorXd>& auxstate) = 0;
        //------------------------------------------------------------------------------
        // Abstract method void Init(const Ref<const VectorXd>& ConstPrm)
        //------------------------------------------------------------------------------
        /**
           * Store subsystem constant parameters contained in SYS_Parameters members
           * in variables specific to the derived class.
           * @note This method is called inside method Setup but is abstract in order to
           * let storing the values of vector SYS_Parameters members in specific members
           * of a derived class
           *
           * @param ConstPrm    Subsystem constant parameters
           * @see SYS_Parameters in method Setup
          */
        //------------------------------------------------------------------------------      
        virtual void Init() = 0;
        
        public:
        /** Boolean variable for subsystem switched on or off by TC.*/
        bool On;
        /** Subsystem on (on = true) or off (on = false) due to operational conditions(e.g. a Sun sensor cannot work when the spacecraft is in eclipse).*/
        bool subsystem_on;
        /** Mode in which a sensor can be used.*/
        string mode;
        
        protected:
        /** Subsystem parameters (@see VarTypes.h).*/
        SYS_params SYS_Parameters;
        /** Subsystem telecommands.*/
        VectorXd TCs;
        /** Subsystem name.*/
        string Name;
        /** Position of subsystem geometric center with respect to the spacecraft body-fixed frame center.*/
        Vec3d Position;
        /** Subsystem constant parameters.*/
        VectorXd ConstPrm;
        /** Subsystem constant parameters.*/
        VectorXd AuxPrm;
        /** First main axis of the subsystem defined in the SC body-fixed frame.*/
        Vec3d X;
        /** Second main axis of the subsystem defined in the SC body-fixed frame.*/
        Vec3d Y;
        /** Third main axis of the subsystem defined in the SC body-fixed frame.*/
        Vec3d Z;
        /** Transformation matrix from SC body-fixed frame subsystem frame.*/
        Mat3x3d SC2SYS;
        /** Operability limits (e.g. field of view limitations, etc.).*/
        VectorXd OPS_limits;
        /** Rate of the output [Hz].*/
        int MaxUpdateRate;
        /** Accuracy (e.g. estimation, actuation, etc.).*/
        VectorXd Accuracy;
        /** Operability range.*/
        double Range;
        /** Space environment models data files ('see VarTypes.h).*/
        EnvModels SpaceEnv;
        /** Normal distribution error of measurement or actuation output.*/
        normal_distribution<double> Error1;
        /** Normal distribution error of measurement or actuation output.*/
        normal_distribution<double> Error2;
        /** Normal distribution error of measurement or actuation output.*/
        normal_distribution<double> Error3;
        /** Random numbers generator.*/
        mt19937 generator1;//(random_device()());
        /** Random numbers generator.*/
        mt19937 generator2;//(random_device()());
        /** Random numbers generator.*/
        mt19937 generator3;//(random_device()());
        ///** Normal distribution error of measurement or actuation output.*/
        //normal_distribution<double> Error[10];
        ///** Normal distribution error of measurement or actuation output.*/
        //mt19937 generator[10];//(random_device()());
        };
    }; // End of namespace subsystem

#endif
