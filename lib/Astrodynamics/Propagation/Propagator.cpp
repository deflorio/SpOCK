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

#include <Propagator.h>
#include <Constants.h>
#include <Transformations.h>

using namespace std;
using namespace math;
using namespace boost::numeric::odeint;
using namespace constants;
using namespace Eigen;

namespace propagator
    {
    //------------------------------------------------------------------------------
    // Class PROP implementation
    //------------------------------------------------------------------------------
    // Initialization of static members
    Mat3x3d PROP::MoI = Mat3x3d::Zero();
    Mat3x3d PROP::invMoI = Mat3x3d::Zero();
    //------------------------------------------------------------------------------
    /**
     * Default constructor
     */
    //------------------------------------------------------------------------------         
    PROP::PROP()
       {
        SC_Parameters = {};
        Models = {};
        
        ggrad_on = false;
        mag_on  = false;
        drag_on = false;
        srp_on = false;
		nMAX = 0;
        sunmoon_on = false;
        
        MoI = Mat3x3d::Zero();
        detMoI = 0.0;
        invMoI = Mat3x3d::Zero();
        //Faces = Vec3d::Zero();
        Mdip = Vec3d::Zero();
        datapath = "";
        planetephem = "";
        magneticfield = "";
        atmosphere = "";
        gravityfield =  "";
        inittime = 0.0;
        CD = 0.0;
        C_SRP = 0.0;
        Area_D = 0.0;
        Area_R = 0.0;
        integ_first_step = true;
        //ECItoBody = Mat3x3d::Zero();
       };
    //------------------------------------------------------------------------------
    /**
     * Constructor
     *
     * @param param    Spacecraft parameters (@see VarTypes.h)
     */
    //------------------------------------------------------------------------------         
    PROP::PROP(SC_params& param)
      {
      Setup(param);
      
      ggrad_on = false;
      mag_on  = false;
      drag_on = false;
      srp_on = false;
      nMAX = 0;
      sunmoon_on = false;
        
      inittime = 0.0;
      integ_first_step = true;
      };
    //------------------------------------------------------------------------------
    /**
     * Constructor
     *
     * @param param    Spacecraft parameters (@see VarTypes.h)
     * @param models   Environment models paths (@see VarTypes.h)
     */
    //------------------------------------------------------------------------------    
    PROP::PROP(SC_params& param,
               EnvModels& models)
      {
      Setup(param, models);
      
      ggrad_on = false;
      mag_on  = false;
      drag_on = false;
      srp_on = false;
      nMAX = 0;
      sunmoon_on = false;
        
      inittime = 0.0;
      integ_first_step = true;
      };  
    //------------------------------------------------------------------------------
    // Destructor
    //------------------------------------------------------------------------------  
    PROP::~PROP() {};
    //------------------------------------------------------------------------------
    // Method void Setup(SC_params& param)
    //------------------------------------------------------------------------------
    /**
     * Store in class spaceraft parameters. Function used by constructors
     *
     * @param param    Spacecraft parameters (@see VarTypes.h)
     */
    //------------------------------------------------------------------------------      
    void PROP::Setup(SC_params& param)
        {
        // Store propagation parameters
        SC_Parameters = param;
        
        //Faces = SC_Parameters.Faces;
        SC_Faces = SC_Parameters.Segment;
        
        SC_mass = SC_Parameters.SC_mass;
        Mdip = SC_Parameters.Mdip;
        CD = SC_Parameters.CD;
        C_SRP = SC_Parameters.C_SRP;
        Area_D = SC_Parameters.Area_D;
        Area_R = SC_Parameters.Area_R;
        MoI = SC_Parameters.MoI;
        
        detMoI = MoI(0,0)*( MoI(1,1)*MoI(2,2) - MoI(2,1)*MoI(1,2) ) - MoI(0,1)*( MoI(1,0)*MoI(2,2) - MoI(2,0)*MoI(1,2) ) + MoI(0,2)*( MoI(1,0)*MoI(2,1) - MoI(2,0)*MoI(1,1) );
        
        invMoI(0,0) = ( MoI(1,1)*MoI(2,2) - MoI(1,2)*MoI(2,1) )/detMoI;
        invMoI(0,1) = ( MoI(0,2)*MoI(2,1) - MoI(0,1)*MoI(2,2) )/detMoI;
        invMoI(0,2) = ( MoI(0,1)*MoI(1,2) - MoI(0,2)*MoI(1,1) )/detMoI;
        
        invMoI(1,0) = ( MoI(1,2)*MoI(2,0) - MoI(1,0)*MoI(2,2) )/detMoI;
        invMoI(1,1) = ( MoI(0,0)*MoI(2,2) - MoI(0,2)*MoI(2,0) )/detMoI;
        invMoI(1,2) = ( MoI(0,2)*MoI(1,0) - MoI(0,0)*MoI(1,2) )/detMoI;
        
        invMoI(2,0) = ( MoI(1,0)*MoI(2,1) - MoI(1,1)*MoI(2,0) )/detMoI;
        invMoI(2,1) = ( MoI(0,1)*MoI(2,0) - MoI(0,0)*MoI(2,1) )/detMoI;
        invMoI(2,2) = ( MoI(0,0)*MoI(1,1) - MoI(0,1)*MoI(1,0) )/detMoI;
       };
    //------------------------------------------------------------------------------
    // Method Setup(SC_params& param, EnvModels& models)
    //------------------------------------------------------------------------------
    /**
     * Store in class spaceraft parameters and environment models paths.
     * Function used by constructors
     *
     * @param param    Spacecraft parameters (@see VarTypes.h)
     * @param models   Environment models paths (@see VarTypes.h)
     */
    //------------------------------------------------------------------------------        
    void PROP::Setup(SC_params& param, EnvModels& models)
        {
        Setup(param);
        
        Models = models;
        
        datapath = Models.datapath;
        planetephem = Models.planetephem;
        magneticfield = Models.magneticfield;
        atmosphere = Models.atmosphere;
        gravityfield = Models.gravityfield;
        };   
    //------------------------------------------------------------------------------
    // Method void Init(double init_time, const Ref<const VectorXd>& init_state)
    //------------------------------------------------------------------------------
    /**
     * Initialization with propagation start epoch and initial state
     *
     * @param init_time    Propagation start epoch
     * @param init_state   Dynamic model initial state vector
     */
    //------------------------------------------------------------------------------  
    void PROP::Init(
                  double init_time,
                  const Ref<const VectorXd>& init_state)
                  {
                  inittime = init_time;
                  initstate = init_state;
                  state = initstate;
                  
                  integ_first_step = true;
                  
                  //ComputeAction(init_time, init_state, init_state);
                  };
    //------------------------------------------------------------------------------
    // Method void Init(double init_time, const Ref<const VectorXd>& init_state, const Ref<const VectorXd>& orb_state)
    //------------------------------------------------------------------------------
    /**
     * Initialization with propagation start epoch, initial state and initial orbital state
     *
     * @param init_time    Propagation start epoch
     * @param init_state   Dynamic model initial state vector
     * @param orb_state    Orbital state vector
     */
    //------------------------------------------------------------------------------  
    void PROP::Init(
                  double init_time,
                  const Ref<const VectorXd>& init_state,
                  const Ref<const VectorXd>& orb_state)
                  {
                  inittime = init_time;
                  initstate = init_state;
                  state = initstate;
                  
                  integ_first_step = true;
                  
                  orbstate = orb_state;
                  
                  //ComputeAction(init_time, init_state, orb_state);                  
                  };
    //------------------------------------------------------------------------------
    // Method void StepperSetup(double eps_abs, double eps_rel, double factor_x, double factor_dxdt)
    //------------------------------------------------------------------------------
    /**
     * Setup numerical integrator parameters
     *
     * @param eps_abs       Absolute tolerance level
     * @param eps_rel       Relative tolerance level
     * @param factor_x      Factor for the weight of the derivative
     * @param factor_dxdt   Factor for the weight of the state
     */
    //------------------------------------------------------------------------------   
    void PROP::StepperSetup(double eps_abs,
                            double eps_rel,
                            double factor_x,
                            double factor_dxdt)
                            {
                            bulirsch_stoer<state_type> setup_stepper(eps_abs, eps_rel, factor_x, factor_dxdt);
                            bulirsch_stoer_stepper = setup_stepper;
                            };
    
}; // End of namespace propagator



