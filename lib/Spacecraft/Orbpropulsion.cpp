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

#include <Orbpropulsion.h>
#include <Constants.h>

//#include <boost/math/distributions/normal.hpp>

using namespace std;
using namespace math;
using namespace constants;
using namespace mathconst;
//using boost::math::normal;
//using namespace spaceenvironment;


namespace orbpropulsion
    {
    //------------------------------------------------------------------------------
    // ORBPROPULSION implementation
    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------
    // Destructor
    //------------------------------------------------------------------------------
    ORBPROPULSION::~ORBPROPULSION() {};
    //------------------------------------------------------------------------------
    // Abstract method void Init(const Ref<const VectorXd>& ConstPrm)
    //------------------------------------------------------------------------------
    /**
     * Store in class members constant parameters contained in base class' members of
     * struct SYS_Parameters
     *
     * @param ConstPrm    Subsystem constant parameters
     * @see SYS_Parameters.ConstPrm in base class method Setup
     */
    //------------------------------------------------------------------------------      
    void ORBPROPULSION::Init()
                        {
                        thrust = ConstPrm(0);
                        thrust_res = ConstPrm(1);
                        thrust_MAX = OPS_limits(0);
                        normal_distribution<double> dt(0.0,PI2);
                        theta_distr = dt;
                        
                        //thrust2dv(SC_mass);
                        };
    //------------------------------------------------------------------------------
    // Method bool status(double epoch, const Ref<const VectorXd>& currentstate, const Ref<const VectorXd>& auxstate)
    //------------------------------------------------------------------------------
    /**
     * Evaluate the operability status (subsystem on or off) and put its value in the class member subsystem_on
     *
     * @param epoch          GPS epoch (seconds) of the input state
     * @param currentstate   Spacecraft attitude state (quaternion)
     * @param auxstate       Spacecraft position vector (ECI)
     */
    //------------------------------------------------------------------------------  
    void ORBPROPULSION::status(double epoch,
                         const Ref<const VectorXd>& currentstate,
                         const Ref<const VectorXd>& auxstate)
                        {
                        subsystem_on = true;
                        };
    //------------------------------------------------------------------------------
    // Method VectorXd Output(double epoch, const Ref<const VectorXd>& currentstate, const Ref<const VectorXd>& auxstate)
    //------------------------------------------------------------------------------
    /**
      * Compute the subsystem output (Output data) given the input data
      *
      * @param epoch          GPS epoch (seconds) of the input state
      * @param currentstate   Spacecraft attitude state (quaternion)
      * @param auxstate       Spacecraft state vector (ECI)
      *
      * @return     dv_CMD maneuver vector in ECI and with the inclusion of the accuracy parameters of the propulsion system in use
      */
    //------------------------------------------------------------------------------  
    VectorXd ORBPROPULSION::Output(double epoch,
                             const Ref<const VectorXd>& currentstate,
                             const Ref<const VectorXd>& auxstate)
                            {
                            Vec3d output = Vec3d::Zero();
                            // Update dv budget
                            dv_MAX = dv_MAX - dv_CMD.norm();
                            if(dv_MAX < 0)
                                {
                                cerr << "Not enough dv budget to execute the orbit maneuver!" << endl;
                                return(output);
                                }
                                
                            Vec3d dv_REAL = Vec3d::Zero();
                            dv_REAL = dv_CMD;
                            
                            // Spacecraft attitude quaternion
                            Vec4d q_attitude;
                            q_attitude = currentstate;
                            // Spacecraft ECI state vector
                            Vector6d state_ECI = auxstate;
                            
                            double dv_acc, dv_norm, dv_newnorm, att_acc, alpha, theta;
                            
                            // Insert propulsion performance error
                            dv_norm = dv_CMD.norm();
                            dv_acc = Error1(generator1);
                            
                            dv_newnorm = dv_norm + dv_norm*dv_acc/100.0; // dv_acc is given in percentage
                            
                            if(dv_acc != 0.0) dv_REAL = vectcorr_norm(dv_REAL, dv_newnorm);
                            
                            // Insert attitude error
                            att_acc = Error2(generator2);
                            // Total accuracy
                            alpha = att_acc;
                            // Spatial circle parameter
                            theta = theta_distr(theta_generator);
                            
                            //cout << "alpha: " << alpha*RAD2DEG << "\n" << endl;
                            //cout << "dv_CMD.norm(): " << dv_CMD.norm() << "\n" << endl;
                            
                            if(alpha != 0.0) dv_REAL = vectcorr_cone(dv_REAL, alpha, theta);
                            
                            //cout << "dv_REAL.norm(): " << dv_REAL.norm() << "\n" << endl;
                            //cout << dv_REAL << "\n" << endl;
                            
                            if(man_sys.compare("RTN") == 0)
                                {
                                Mat3x3d ECI2RTN = Mat3x3d::Zero();
                                ECI2RTN = ECI2RTN_Matrix(state_ECI);
                                //output = (ECI2RTN.transpose())*dv_REAL;
                                output = RTN2ECI(dv_REAL, state_ECI);
                                }
                            else if(man_sys.compare("Body") == 0)
                                {
                                Vec4d q_currentstate_inv = Vec4d::Zero();    
                                q_currentstate_inv = q_inv(q_attitude);
                                output = TransbyQ(dv_REAL,q_currentstate_inv);
                                }
                                
                            return(output);
                            };
    //------------------------------------------------------------------------------
    // Method Vec3d thrust2dv(double SC_mass)
    //------------------------------------------------------------------------------
    /**
     * Convert thrust parameters of propulsion system to dv parameters given a spacecraft mass
     *
     * @param SC_mass    Spacecraft mass [kg]
     */
    //------------------------------------------------------------------------------  
    void ORBPROPULSION::thrust2dv(double SC_mass)
                        {
                        double dt, dt_MAX;
                        
                        dt = 1.0; // [s]
                        
                        dv_thr = thrust*dt/SC_mass; // dv obtained by thrusting for 1 s
                        
                        dv_res = dv_thr*(thrust_res/thrust);
                        
                        dt_MAX = thrust_MAX/thrust;
                        
                        dv_MAX = dv_thr*dt_MAX;
                        
                        //dv_prms << dv_thr, dv_res, dv_MAX;
                        //
                        //return(dv_prms);
                        };
    //------------------------------------------------------------------------------
    // Method void impman2contman(vector<maneuver> &impman)
    //------------------------------------------------------------------------------
    /**
     * Convert an orbit impulsive maneuvers vector to a continuous maneuvers vector given the propulsion system parameters
     * and the spacecraft mass and the integration step mass.
     * The dvs in impman are procesed in a way that if a dv is larger than the one dv_STEP that can be obtained in one integration
     * step (e.g. 1 s, 10 s, etc.) it is splitted in slots of dv_STEP which will be executed at each integration step as a continuous maneuver
     * (e.g. if dv = 1 m/s, dv_STEP = 1e-2 m/s and SIM_STEP = 10 s, dv is splitted in 10 dv = 1e-1,
     * if dv = 1 m/s, dv_STEP = 1e-2 m/s and SIM_STEP = 1 s, dv is splitted in 100 dv = 1e-2, etc.)
     *
     * @param impman     Vector of impulsive maneuvers [m/s] (see VarTypes.h)
     * @param SC_mass    Spacecraft mass [kg]
     * @param SIM_STEP   Integration step size of numerical propagation [s]
     */
    //------------------------------------------------------------------------------  
    void ORBPROPULSION::impman2contman(vector<maneuver> &impman,
                                       int SIM_STEP)
                                        {
                                        unsigned int man_ind;
                                        vector<maneuver>::size_type man_size;
                                        man_size = impman.size();
                                        double m_init_time;
                                        Vec3d dv = Vec3d::Zero();
                                        //Vec3d dv_prms = Vec3d::Zero();
                                        double dv_unit, minus_dv_unit;
                                        Vec3d unit_dv = Vec3d::Zero();
                                        maneuver unit_maneuver;
                                        
                                        //dv_prms = thrust2dv(SC_mass);
                                        //dv_unit = SIM_STEP*dv_prms(0); // dv obtained by thrusting for SIM_STEP s. A maneuver with dv < dv_unit is considered as impulsive. The most accurate simuation is
                                        //                                // given by choosing SIM_STEP = 1 s, since in this case dv_unit = dv_thr (see function thrust2dv)
                                        //dv_res = dv_prms(1);
                                        //dv_MAX = dv_prms(2);
                                        
                                        dv_unit = SIM_STEP*dv_thr; // dv obtained by thrusting for SIM_STEP s. A maneuver with dv < dv_unit is considered as impulsive. The most accurate simuation is
                                                                        // given by choosing SIM_STEP = 1 s, since in this case dv_unit = dv_thr (see function thrust2dv)
                                        minus_dv_unit = -dv_unit;
                                                                        
                                        
                                        for( man_ind = 0; man_ind < man_size; man_ind++ )
                                            {
                                            if( man_size > 1 && (impman[man_ind].init_time <= impman[man_ind - 1].init_time) )
                                            {
                                            cerr << "ERROR: initial time of a continuous maneuver < final time of previous continuous maneuver" << endl;
                                            return;
                                            }
                                            
                                            m_init_time = impman[man_ind].init_time;
                                            dv = impman[man_ind].ManVec;
                                            
                                            if(dv_res != 0.0) // Take into account the dv resolution if this parameter is given in the propulsion system characterization
                                                {
                                                if( dv.norm() < dv_res ) impman[man_ind].ManVec = Vec3d::Zero();
                                                //for(int i = 0 ; i < 3; i++)
                                                //    {
                                                //    if(fabs(dv(i)) < dv_res) dv(i) = 0.0; //dv(i) = dv_res*round(dv(i)/dv_res);
                                                //    }
                                                }
                                            
                                            //if( (dv(0) > dv_unit) || (dv(1) > dv_unit) || (dv(2) > dv_unit) )
                                            if( dv.norm() >= dv_unit )
                                                {
                                                unit_maneuver = impman[man_ind];
                                                
                                                m_init_time = m_init_time + SIM_STEP;
                                                
                                                for(int i = 0 ; i < 3; i++)
                                                    {
                                                    //m_init_time = m_init_time + SIM_STEP;
                                                       
                                                    if( dv(i) > 0.0 && dv(i) > dv_unit )
                                                      {
                                                      dv(i) = dv(i) - dv_unit;
                                                      unit_dv(i) = dv_unit;
                                                      }
                                                    else if( dv(i) < 0.0 && fabs(dv(i)) > dv_unit )
                                                      {
                                                      dv(i) = dv(i) - minus_dv_unit;
                                                      unit_dv(i) = minus_dv_unit;
                                                      }
                                                    else
                                                      {
                                                      unit_dv(i) = dv(i);
                                                      }
                                                    }
                                                
                                                if( ( dv.norm() >= dv_unit ) || ( dv.norm() < dv_unit && dv.norm() > dv_res ) )
                                                    {    
                                                    // Make place for the insertion of a new element
                                                    man_size = man_size + 1;
                                                    impman.resize(man_size);
                                                    // Translate up the element over position man_ind
                                                    for( unsigned int i = man_size-1; i > man_ind; i-- ) impman[i] = impman[i-1];
                                                    // Change the element that was in position man_ind
                                                    impman[man_ind+1].init_time = m_init_time;
                                                    impman[man_ind+1].ManVec = dv;
                                                    }
                                                // Insert new element in position man_ind
                                                unit_maneuver.ManVec = unit_dv;
                                                impman[man_ind] = unit_maneuver; 
                                                }
                                            }
                                        };
    //------------------------------------------------------------------------------
    

}; // End of namespace orbpropulsion



