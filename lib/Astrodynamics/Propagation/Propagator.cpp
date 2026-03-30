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

#include <boost/math/interpolators/cubic_hermite.hpp>
#include <boost/math/interpolators/quintic_hermite.hpp>

using namespace std;
using namespace math;
using namespace boost::numeric::odeint;
using namespace constants;
using namespace Eigen;
using namespace boost::math::interpolators;

namespace propagator
    {
    //------------------------------------------------------------------------------
    // Class PROP implementation
    //------------------------------------------------------------------------------
    // Initialization of static members
    Mat3x3d PROP::MoI = Mat3x3d::Zero();
    Mat3x3d PROP::invMoI = Mat3x3d::Zero();
    
    Vec3d PROP::eop = Vec3d::Zero();
    double PROP::tai_utc = 0.0;
    //------------------------------------------------------------------------------
    /**
     * Default constructor
     */
    //------------------------------------------------------------------------------         
    PROP::PROP()
       {
        SC_Parameters.MoI = Mat3x3d::Zero();
        SC_Parameters.Mdip = Vec3d::Zero();
        SC_Parameters.SC_mass = 0.0;
        SC_Parameters.CD = 0.0;
        SC_Parameters.C_SRP = 0.0;
        SC_Parameters.Area_D = 0.0;
        SC_Parameters.Area_R = 0.0;
        
        ODEINT_stepper = "RK4";
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
        sunmoon = "";
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
        sunmoon = Models.sunmoon;
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
                  state = Vec4d::Zero();
                  orbstate = init_state;
                  
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
     * 
     */
    //------------------------------------------------------------------------------   
    //void PROP::StepperSetup(double eps_abs,
    //                        double eps_rel,
    //                        double factor_x,
    //                        double factor_dxdt)
    //                        {
    //                        bulirsch_stoer<state_type> setup_stepper(eps_abs, eps_rel, factor_x, factor_dxdt);
    //                        bulirsch_stoer_stepper = setup_stepper;
    //                        };
    ////------------------------------------------------------------------------------
    //// void QH_Interpolation(const int InterpStep, const Ref<const MatrixXd>& timeposvelacc, Ref<MatrixXd> orbstate_interpolated)
    ////------------------------------------------------------------------------------
    ///**
    // * Perform quintic Hermite interpolation of state vector (x, y, z, dx/dt, dy/dt, dz/dt) over a time interval and with a specific interpolation step
    // *
    // * @param InterpStep            Step of interpolated points                 
    // * @param timeposvelacc         Matrix of points used for the computation of the interpolation polynomial
    // *                              (each row contains time, x, y, z, vx, vy, vz, ax, ay, az)
    // * 
    // * @return Matrix of interpolated points (each row contains time, x, y, z, vx, vy, vz)
    // */
    ////------------------------------------------------------------------------------   
    //void PROP::QH_Interpolation(const int InterpStep,
    //                            const Ref<const MatrixXd>& timeposvelacc,
    //                            Ref<MatrixXd> orbstate_interpolated)
    //                            {
    //                            const int InterPoints = timeposvelacc.rows();
    //                            double t_intpl;
    //                            int ind = 1;
    //                            
    //                            vector<double> v_null(InterPoints);
    //                            for(int i = 0; i < InterPoints; i++) v_null[i] = 0.0;
    //                            
    //                            vector<double> time_interpol(InterPoints);
    //                            
    //                            vector<double> x_interpol(InterPoints);
    //                            vector<double> dxdt_interpol(InterPoints);
    //                            vector<double> dx2dt2_interpol(InterPoints);
    //                            
    //                            int lastrow = orbstate_interpolated.rows() - 1;
    //                            
    //                            orbstate_interpolated.row(0) = timeposvelacc.row(0).segment(0,7);
    //                            
    //                            for(int i = 0; i < InterPoints; i++)
    //                                {
    //                                time_interpol[i] = timeposvelacc(i,0);
    //                                x_interpol[i] = timeposvelacc(i,1);
    //                                dxdt_interpol[i] = timeposvelacc(i,4);
    //                                dx2dt2_interpol[i] = timeposvelacc(i,7);
    //                                }
    //                            
    //                            auto spline_x = quintic_hermite(std::move(time_interpol), std::move(x_interpol), std::move(dxdt_interpol), std::move(dx2dt2_interpol));
    //                            
    //                            time_interpol.clear();   time_interpol = v_null;
    //                            x_interpol.clear();      x_interpol = v_null;
    //                            dxdt_interpol.clear();   dxdt_interpol = v_null;
    //                            dx2dt2_interpol.clear(); dx2dt2_interpol = v_null;
    //                            
    //                            for(int i = 0; i < InterPoints; i++)
    //                                {
    //                                time_interpol[i] = timeposvelacc(i,0);
    //                                x_interpol[i] = timeposvelacc(i,2);
    //                                dxdt_interpol[i] = timeposvelacc(i,5);
    //                                dx2dt2_interpol[i] = timeposvelacc(i,8);
    //                                }
    //                            
    //                            auto spline_y = quintic_hermite(std::move(time_interpol), std::move(x_interpol), std::move(dxdt_interpol), std::move(dx2dt2_interpol));
    //                            
    //                            time_interpol.clear();   time_interpol = v_null;
    //                            x_interpol.clear();      x_interpol = v_null;
    //                            dxdt_interpol.clear();   dxdt_interpol = v_null;
    //                            dx2dt2_interpol.clear(); dx2dt2_interpol = v_null;
    //                            
    //                            for(int i = 0; i < InterPoints; i++)
    //                                {
    //                                time_interpol[i] = timeposvelacc(i,0);
    //                                x_interpol[i] = timeposvelacc(i,3);
    //                                dxdt_interpol[i] = timeposvelacc(i,6);
    //                                dx2dt2_interpol[i] = timeposvelacc(i,9);
    //                                }
    //                            
    //                            auto spline_z = quintic_hermite(std::move(time_interpol), std::move(x_interpol), std::move(dxdt_interpol), std::move(dx2dt2_interpol));
    //                            
    //                            time_interpol.clear();   time_interpol = v_null;
    //                            x_interpol.clear();      x_interpol = v_null;
    //                            dxdt_interpol.clear();   dxdt_interpol = v_null;
    //                            dx2dt2_interpol.clear(); dx2dt2_interpol = v_null;
    //                            
    //                            
    //                            t_intpl = orbstate_interpolated(0,0) + InterpStep;
    //                            
    //                            for(int i = 1; i < lastrow; i++)
    //                                {
    //                                if( t_intpl == timeposvelacc(ind,0) )
    //                                    {
    //                                     orbstate_interpolated.row(i) = timeposvelacc.row(ind).segment(0,7);
    //                                     
    //                                     t_intpl += InterpStep;
    //                                     ind++;
    //                                     //continue;
    //                                    }
    //                                else
    //                                    {
    //                                    orbstate_interpolated(i,0) = t_intpl;
    //                                    orbstate_interpolated(i,1) = spline_x(t_intpl);
    //                                    orbstate_interpolated(i,2) = spline_y(t_intpl);
    //                                    orbstate_interpolated(i,3) = spline_z(t_intpl);
    //                                    orbstate_interpolated(i,4) = spline_x.prime(t_intpl);
    //                                    orbstate_interpolated(i,5) = spline_y.prime(t_intpl);
    //                                    orbstate_interpolated(i,6) = spline_z.prime(t_intpl);
    //                                    
    //                                    t_intpl += InterpStep;
    //                                    }
    //                                }
    //                                
    //                            orbstate_interpolated.row(lastrow) = timeposvelacc.row(InterPoints-1).segment(0,7);
    //                            };
    ////------------------------------------------------------------------------------
    //// void QH_Interpolation(double t_intpl, const Ref<const MatrixXd>& timeposvelacc, VectorNd<9>& orbstate_interpolated)
    ////------------------------------------------------------------------------------
    ///**
    // * Perform quintic Hermite interpolation of state vector and acceleration (x, y, z, dx/dt, dy/dt, dz/dt, dx2/dt2, dy2/dt2, dz2/dt2) at a specific time
    // *
    // * @param t_intpl               Interpolation time
    // * @param timeposvelacc         Matrix of points used for the computation of the interpolation polynomial
    // *                              (each row contains time, x, y, z, vx, vy, vz, ax, ay, az)
    // * 
    // * @return Interpolated state vector and acceleration (x, y, z, dx/dt, dy/dt, dz/dt, dx2/dt2, dy2/dt2, dz2/dt2)
    // */
    ////------------------------------------------------------------------------------   
    //void PROP::QH_Interpolation(double t_intpl,
    //                            const Ref<const MatrixXd>& timeposvelacc,
    //                            VectorNd<9>& orbstate_interpolated)
    //                            {
    //                            const int InterPoints = timeposvelacc.rows();
    //                            
    //                            for(int i = 0; i < InterPoints; i++)
    //                                {
    //                                if(t_intpl == timeposvelacc(i,0)) // No need of interpolation
    //                                    {
    //                                    orbstate_interpolated = timeposvelacc.row(i).segment(1,9);
    //                                    return;
    //                                    }
    //                                }
    //                            
    //                            vector<double> v_null(InterPoints);
    //                            for(int i = 0; i < InterPoints; i++) v_null[i] = 0.0;
    //                            
    //                            vector<double> time_interpol(InterPoints);
    //                            vector<double> x_interpol(InterPoints);
    //                            vector<double> dxdt_interpol(InterPoints);
    //                            vector<double> dx2dt2_interpol(InterPoints);
    //                            
    //                            for(int i = 0; i < InterPoints; i++)
    //                                {
    //                                time_interpol[i] = timeposvelacc(i,0);
    //                                x_interpol[i] = timeposvelacc(i,1);
    //                                dxdt_interpol[i] = timeposvelacc(i,4);
    //                                dx2dt2_interpol[i] = timeposvelacc(i,7);
    //                                }
    //                            
    //                            auto spline_x = quintic_hermite(std::move(time_interpol), std::move(x_interpol), std::move(dxdt_interpol), std::move(dx2dt2_interpol));
    //                            
    //                            time_interpol.clear();   time_interpol = v_null;
    //                            x_interpol.clear();      x_interpol = v_null;
    //                            dxdt_interpol.clear();   dxdt_interpol = v_null;
    //                            dx2dt2_interpol.clear(); dx2dt2_interpol = v_null;
    //                            
    //                            for(int i = 0; i < InterPoints; i++)
    //                                {
    //                                time_interpol[i] = timeposvelacc(i,0);
    //                                x_interpol[i] = timeposvelacc(i,2);
    //                                dxdt_interpol[i] = timeposvelacc(i,5);
    //                                dx2dt2_interpol[i] = timeposvelacc(i,8);
    //                                }
    //                            
    //                            auto spline_y = quintic_hermite(std::move(time_interpol), std::move(x_interpol), std::move(dxdt_interpol), std::move(dx2dt2_interpol));
    //                            
    //                            time_interpol.clear();   time_interpol = v_null;
    //                            x_interpol.clear();      x_interpol = v_null;
    //                            dxdt_interpol.clear();   dxdt_interpol = v_null;
    //                            dx2dt2_interpol.clear(); dx2dt2_interpol = v_null;
    //                            
    //                            for(int i = 0; i < InterPoints; i++)
    //                                {
    //                                time_interpol[i] = timeposvelacc(i,0);
    //                                x_interpol[i] = timeposvelacc(i,3);
    //                                dxdt_interpol[i] = timeposvelacc(i,6);
    //                                dx2dt2_interpol[i] = timeposvelacc(i,9);
    //                                }
    //                            
    //                            auto spline_z = quintic_hermite(std::move(time_interpol), std::move(x_interpol), std::move(dxdt_interpol), std::move(dx2dt2_interpol));
    //                            
    //                            time_interpol.clear();   time_interpol = v_null;
    //                            x_interpol.clear();      x_interpol = v_null;
    //                            dxdt_interpol.clear();   dxdt_interpol = v_null;
    //                            dx2dt2_interpol.clear(); dx2dt2_interpol = v_null;
    //                            
    //                            // Interpolation                                
    //                            orbstate_interpolated(0) = spline_x(t_intpl);
    //                            orbstate_interpolated(1) = spline_y(t_intpl);
    //                            orbstate_interpolated(2) = spline_z(t_intpl);
    //                            orbstate_interpolated(3) = spline_x.prime(t_intpl);
    //                            orbstate_interpolated(4) = spline_y.prime(t_intpl);
    //                            orbstate_interpolated(5) = spline_z.prime(t_intpl);
    //                            orbstate_interpolated(6) = spline_x.double_prime(t_intpl);
    //                            orbstate_interpolated(7) = spline_y.double_prime(t_intpl);
    //                            orbstate_interpolated(8) = spline_z.double_prime(t_intpl);
    //                            };
    ////------------------------------------------------------------------------------
    //// void CH_Interpolation(const int InterPoints, const int InterpStep, const Ref<const MatrixXd>& timeposvel, Ref<MatrixXd> orbstate_interpolated)
    ////------------------------------------------------------------------------------
    ///**
    // * Perform quintic Hermite interpolation of state vector (x, y, z, dx/dt, dy/dt, dz/dt)
    // *
    // * @param InterPoints           Number of points used for the computation of the interpolation polynomial
    // * @param InterpStep            Step of interpolated points                 
    // * @param timeposvel            Matrix of points used for the computation of the interpolation polynomial
    // *                              (each row contains time, x, y, z, vx, vy, vz)
    // * 
    // * @return Matrix of interpolated points (each row contains time, x, y, z, vx, vy, vz)  
    // */
    ////------------------------------------------------------------------------------   
    //void PROP::CH_Interpolation(const int InterPoints,
    //                            const int InterpStep,
    //                            const Ref<const MatrixXd>& timeposvel,
    //                            Ref<MatrixXd> orbstate_interpolated)
    //                            {
    //                            double t_intpl;
    //                            int ind = 1;
    //                            
    //                            vector<double> v_null(InterPoints);
    //                            for(int i = 0; i < InterPoints; i++) v_null[i] = 0.0;
    //                            
    //                            vector<double> time_interpol(InterPoints);
    //                            
    //                            vector<double> x_interpol(InterPoints);
    //                            vector<double> dxdt_interpol(InterPoints);
    //                            
    //                            int lastrow = orbstate_interpolated.rows() - 1;
    //                            
    //                            orbstate_interpolated.row(0) = timeposvel.row(0).segment(0,7);
    //                            
    //                            //cout << orbstate_interpolated(0,0) << "   " << orbstate_interpolated(0,1) << "   " << orbstate_interpolated(0,2) << "   " << orbstate_interpolated(0,3) << "   " << orbstate_interpolated(0,4) << "   " << orbstate_interpolated(0,5) << "   " << orbstate_interpolated(0,6) << endl;
    //                            
    //                            for(int i = 0; i < InterPoints; i++)
    //                                {
    //                                time_interpol[i] = timeposvel(i,0);
    //                                x_interpol[i] = timeposvel(i,1);
    //                                dxdt_interpol[i] = timeposvel(i,4);
    //                                }
    //                            
    //                            auto spline_x = cubic_hermite(std::move(time_interpol), std::move(x_interpol), std::move(dxdt_interpol));
    //                            
    //                            time_interpol.clear();   time_interpol = v_null;
    //                            x_interpol.clear();      x_interpol = v_null;
    //                            dxdt_interpol.clear();   dxdt_interpol = v_null;
    //                            
    //                            for(int i = 0; i < InterPoints; i++)
    //                                {
    //                                time_interpol[i] = timeposvel(i,0);
    //                                x_interpol[i] = timeposvel(i,2);
    //                                dxdt_interpol[i] = timeposvel(i,5);
    //                                }
    //                            
    //                            auto spline_y = cubic_hermite(std::move(time_interpol), std::move(x_interpol), std::move(dxdt_interpol));
    //                            
    //                            time_interpol.clear();   time_interpol = v_null;
    //                            x_interpol.clear();      x_interpol = v_null;
    //                            dxdt_interpol.clear();   dxdt_interpol = v_null;
    //                            
    //                            for(int i = 0; i < InterPoints; i++)
    //                                {
    //                                time_interpol[i] = timeposvel(i,0);
    //                                x_interpol[i] = timeposvel(i,3);
    //                                dxdt_interpol[i] = timeposvel(i,6);
    //                                }
    //                            
    //                            auto spline_z = cubic_hermite(std::move(time_interpol), std::move(x_interpol), std::move(dxdt_interpol));
    //                            
    //                            time_interpol.clear();   time_interpol = v_null;
    //                            x_interpol.clear();      x_interpol = v_null;
    //                            dxdt_interpol.clear();   dxdt_interpol = v_null;
    //                            
    //                            
    //                            t_intpl = orbstate_interpolated(0,0) + InterpStep;
    //                            
    //                            for(int i = 1; i < lastrow; i++)
    //                                {
    //                                //cout << "t_intpl = " << t_intpl << endl;
    //                                //cout << "timeposvel(ind,0) = " << timeposvel(ind,0) << endl;
    //                                if( t_intpl == timeposvel(ind,0) )
    //                                    {
    //                                     orbstate_interpolated.row(i) = timeposvel.row(ind).segment(0,7);
    //                                     
    //                                     t_intpl += InterpStep;
    //                                     ind++;
    //                                     //continue;
    //                                    }
    //                                else
    //                                    {
    //                                    orbstate_interpolated(i,0) = t_intpl;
    //                                    
    //                                    orbstate_interpolated(i,1) = spline_x(t_intpl);
    //                                    orbstate_interpolated(i,4) = spline_x.prime(t_intpl);
    //                                    
    //                                    orbstate_interpolated(i,2) = spline_y(t_intpl);
    //                                    orbstate_interpolated(i,5) = spline_y.prime(t_intpl);
    //                                    
    //                                    orbstate_interpolated(i,3) = spline_z(t_intpl);
    //                                    orbstate_interpolated(i,6) = spline_z.prime(t_intpl);
    //                                    
    //                                    //t_intpl++;
    //                                    t_intpl += InterpStep;
    //                                    }
    //                                
    //                                //cout << orbstate_interpolated(i,0) << "   " << orbstate_interpolated(i,1) << "   " << orbstate_interpolated(i,2) << "   " << orbstate_interpolated(i,3) << "   " << orbstate_interpolated(i,4) << "   " << orbstate_interpolated(i,5) << "   " << orbstate_interpolated(i,6) << endl;
    //                                }
    //                                
    //                            orbstate_interpolated.row(lastrow) = timeposvel.row(InterPoints-1).segment(0,7);
    //                                
    //                            //cout << orbstate_interpolated(lastrow,0) << "   " << orbstate_interpolated(lastrow,1) << "   " << orbstate_interpolated(lastrow,2) << "   " << orbstate_interpolated(lastrow,3) << "   " << orbstate_interpolated(lastrow,4) << "   " << orbstate_interpolated(lastrow,5) << "   " << orbstate_interpolated(lastrow,6) << endl;   
    //                            };
    ////------------------------------------------------------------------------------
    //// void CH_Interpolation(double t_intpl, const Ref<const MatrixXd>& timeposvel, VectorNd<7>& orbstate_interpolated)
    ////------------------------------------------------------------------------------
    ///**
    // * Perform cubic Hermite interpolation of state vector (x, y, z, dx/dt, dy/dt, dz/dt) at a specific time
    // *
    // * @param t_intpl               Interpolation time
    // * @param timeposvel            Matrix of points used for the computation of the interpolation polynomial
    // *                              (each row contains time, x, y, z, vx, vy, vz)
    // * 
    // * @return Interpolated state vector (x, y, z, vx, vy, vz)
    // */
    ////------------------------------------------------------------------------------   
    //void PROP::CH_Interpolation(double t_intpl,
    //                            const Ref<const MatrixXd>& timeposvel,
    //                            VectorNd<6>& orbstate_interpolated)
    //                            {
    //                            const int InterPoints = timeposvel.rows();
    //                            
    //                            for(int i = 0; i < InterPoints; i++)
    //                                {
    //                                if(t_intpl == timeposvel(i,0)) // No need of interpolation
    //                                    {
    //                                    orbstate_interpolated = timeposvel.row(i).segment(1,6);
    //                                    return;
    //                                    }
    //                                }
    //                            
    //                            vector<double> v_null(InterPoints);
    //                            for(int i = 0; i < InterPoints; i++) v_null[i] = 0.0;
    //                            
    //                            vector<double> time_interpol(InterPoints);
    //                            vector<double> x_interpol(InterPoints);
    //                            vector<double> dxdt_interpol(InterPoints);
    //                            
    //                            for(int i = 0; i < InterPoints; i++)
    //                                {
    //                                time_interpol[i] = timeposvel(i,0);
    //                                x_interpol[i] = timeposvel(i,1);
    //                                dxdt_interpol[i] = timeposvel(i,4);
    //                                }
    //                            
    //                            auto spline_x = cubic_hermite(std::move(time_interpol), std::move(x_interpol), std::move(dxdt_interpol));
    //                            
    //                            time_interpol.clear();   time_interpol = v_null;
    //                            x_interpol.clear();      x_interpol = v_null;
    //                            dxdt_interpol.clear();   dxdt_interpol = v_null;
    //                            
    //                            for(int i = 0; i < InterPoints; i++)
    //                                {
    //                                time_interpol[i] = timeposvel(i,0);
    //                                x_interpol[i] = timeposvel(i,2);
    //                                dxdt_interpol[i] = timeposvel(i,5);
    //                                }
    //                            
    //                            auto spline_y = cubic_hermite(std::move(time_interpol), std::move(x_interpol), std::move(dxdt_interpol));
    //                            
    //                            time_interpol.clear();   time_interpol = v_null;
    //                            x_interpol.clear();      x_interpol = v_null;
    //                            dxdt_interpol.clear();   dxdt_interpol = v_null;
    //                            
    //                            for(int i = 0; i < InterPoints; i++)
    //                                {
    //                                time_interpol[i] = timeposvel(i,0);
    //                                x_interpol[i] = timeposvel(i,3);
    //                                dxdt_interpol[i] = timeposvel(i,6);
    //                                }
    //                            
    //                            auto spline_z = cubic_hermite(std::move(time_interpol), std::move(x_interpol), std::move(dxdt_interpol));
    //                            
    //                            time_interpol.clear();   time_interpol = v_null;
    //                            x_interpol.clear();      x_interpol = v_null;
    //                            dxdt_interpol.clear();   dxdt_interpol = v_null;
    //                            
    //                            // Interpolation                                
    //                            orbstate_interpolated(0) = spline_x(t_intpl);
    //                            orbstate_interpolated(1) = spline_y(t_intpl);
    //                            orbstate_interpolated(2) = spline_z(t_intpl);
    //                            orbstate_interpolated(3) = spline_x.prime(t_intpl);
    //                            orbstate_interpolated(4) = spline_y.prime(t_intpl);
    //                            orbstate_interpolated(5) = spline_z.prime(t_intpl);
    //                            };                                
                                
                                
                                
                                
                                
                                
                                
                                
    
}; // End of namespace propagator



