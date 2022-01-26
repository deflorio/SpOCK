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
#include <string>
#include <math.h>
#include <random>

#include <Subsystem.h>
//#include <Constants.h>

using namespace std;
//using namespace math;
//using namespace constants;

namespace subsystem
    {
    //------------------------------------------------------------------------------
    // SUBSYS implementation
    //------------------------------------------------------------------------------
    // Initialization of static members
    //string Name = " ";
    //Vec3d SUBSYS::Position = Vec3d::Zero();
    //Mat3x3d SUBSYS::SC2SYS = Mat3x3d::Zero();
    //Vec3d SUBSYS::MainAxis1 = Vec3d::Zero();
    //Vec3d SUBSYS::MainAxis2 = Vec3d::Zero();
    //Vec3d SUBSYS::MainAxis3 = Vec3d::Zero();
    
    //------------------------------------------------------------------------------
    /**
     * Default constructor
     */
    //------------------------------------------------------------------------------         
    SUBSYS::SUBSYS()
        {
        //SYS_Parameters = {};
        SYS_Parameters.SC2SYS = Mat3x3d::Zero();
        SYS_Parameters.Position = Vec3d::Zero();
        SYS_Parameters.Name = "";
        SYS_Parameters.Range = 0.0;
        SYS_Parameters.MaxUpdateRate = 0;
        SYS_Parameters.on_off = false;
        
        On = true; // This value will be passe to member On by TC in a later version of the software
        subsystem_on = true;
        mode = " ";
        X = Vec3d::Zero();
        Y = Vec3d::Zero();
        Z = Vec3d::Zero();
        };
    //------------------------------------------------------------------------------
    /**
     * Constructor
     *
     * @param param          Spacecraft parameters (@see VarTypes.h)
     */
    //------------------------------------------------------------------------------         
    SUBSYS::SUBSYS(SYS_params& param)
      {
      Setup(param);
      
      subsystem_on = true;
      mode = " ";
      };
    //------------------------------------------------------------------------------
    /**
     * Constructor
     *
     * @param param          Spacecraft parameters (@see VarTypes.h)
     * @param Telecommands   Subsystem telecommands
     */
    //------------------------------------------------------------------------------         
    SUBSYS::SUBSYS(SYS_params& param,
                   const Ref<const VectorXd>& Telecommands)
                  {
                  Setup(param);
                  
                  subsystem_on = true;
                  mode = " ";
                  
                  TCs = Telecommands;
                  };
    //------------------------------------------------------------------------------
    // Destructor
    //------------------------------------------------------------------------------  
    SUBSYS::~SUBSYS() {};
    //------------------------------------------------------------------------------
    // Method void Setup(SYS_params& param)
    //------------------------------------------------------------------------------
    /**
     * Store in class spaceraft parameters. Function used by constructors
     *
     * @param param    Spacecraft parameters (@see VarTypes.h)
     */
    //------------------------------------------------------------------------------      
    void SUBSYS::Setup(SYS_params& param)
        {
        // Store subsystem's parameters
        SYS_Parameters = param;
        
        On = SYS_Parameters.on_off;
        Name = SYS_Parameters.Name;
        Position = SYS_Parameters.Position;
        ConstPrm = SYS_Parameters.ConstPrm;
        AuxPrm = SYS_Parameters.AuxPrm;
        SC2SYS = SYS_Parameters.SC2SYS;
        OPS_limits = SYS_Parameters.OPS_limits;
        MaxUpdateRate = SYS_Parameters.MaxUpdateRate;
        Accuracy = SYS_Parameters.Accuracy;
        Range = SYS_Parameters.Range;
        SpaceEnv = SYS_Parameters.SpaceEnv;
        
        X = SC2SYS.row(0);
        Y = SC2SYS.row(1);
        Z = SC2SYS.row(2);
        
        switch(Accuracy.size())
            {
            case 2: {
                    double mean1 = Accuracy(0);
                    double stddev1 = Accuracy(1);
        
                    normal_distribution<double> d1(mean1,stddev1);
                    Error1 = d1;
                    
                    break;
                    }
                    
            case 4: {
                    double mean1 = Accuracy(0);
                    double stddev1 = Accuracy(1);
                    double mean2 = Accuracy(2);
                    double stddev2 = Accuracy(3);
                    
                    normal_distribution<double> d1(mean1,stddev1);
                    Error1 = d1;
                    
                    normal_distribution<double> d2(mean2,stddev2);
                    Error2 = d2;
                    
                    break;
                    }
                    
            case 6: {
                    double mean1 = Accuracy(0);
                    double stddev1 = Accuracy(1);
                    double mean2 = Accuracy(2);
                    double stddev2 = Accuracy(3);
                    double mean3 = Accuracy(4);
                    double stddev3 = Accuracy(5);
                    
                    normal_distribution<double> d1(mean1,stddev1);
                    Error1 = d1;
                    
                    normal_distribution<double> d2(mean2,stddev2);
                    Error2 = d2;
                    
                    normal_distribution<double> d3(mean3,stddev3);
                    Error3 = d3;
                    
                    break;
                    }
            }
        
        };
    

}; // End of namespace spaceenvironment



