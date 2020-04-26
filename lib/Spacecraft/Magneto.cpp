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

#include <Magneto.h>
#include <Constants.h>

//#include <boost/math/distributions/normal.hpp>

using namespace std;
using namespace math;
using namespace constants;
using namespace mathconst;

namespace magneto
    {
    //------------------------------------------------------------------------------
    // MAGNETO implementation
    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------
    // Destructor
    //------------------------------------------------------------------------------
    MAGNETO::~MAGNETO() {};
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
    void MAGNETO::Init()
                        {
                        magneticfield = SpaceEnv.magneticfield;
                        MagneticField.setfilespath(magneticfield);
                        ontime2dipole = ConstPrm(0);    
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
    void MAGNETO::status(double epoch,
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
      * @param auxstate       Spacecraft position vector (ECEF)
      *
      * @return Readings of magnetometer or magnetic torquers' torque vector depending on the value of class member 'Name'.
      *         Output of the magnetometer is the estimated Earth magnetic field [nanotesla] vector (magnetometer frame).
      *         Output of the magnetic torquers is the torque vector [Nm] produced by the actuator (SC body-fixed frame).
      */
    //------------------------------------------------------------------------------  
    VectorXd MAGNETO::Output(double epoch,
                             const Ref<const VectorXd>& currentstate,
                             const Ref<const VectorXd>& auxstate)
                            {
                            Vec3d output = Vec3d::Zero();
                            //normal_distribution<double> Error1 = Error[0];
                            //mt19937 generator1 = generator[0];
                            
                            Vec4d q_currentstate = currentstate;
                            Vec3d pos_ECEF = auxstate;
                            
                            Vec3d B_ECI = Vec3d::Zero();
                            Vec3d B_Body = Vec3d::Zero();
                            
                            MagneticField.SetReferenceFrame("ECI");
                            try{ B_ECI = MagneticField.field_vec(epoch, pos_ECEF); }
                            catch(const string errmsg)
                            {
                            cerr << "In MAGNETO::Output " + errmsg << endl;
                            exit(EXIT_FAILURE);
                            }
                            
                            B_Body = TransbyQ(B_ECI,q_currentstate);
                            
                            //B_Body = B_Body*1E-9; // Conversion from nanotesla to tesla
                            
                            if(On)
                              {
                              double err;
                                
                              if(Name.compare("Magnetometer") == 0)
                                {
                                output = SC2SYS*B_Body;
                                
                                for(int i = 0; i < 3; i++)
                                    {
                                    err = Error1(generator1); // Measurement error
                                    output(i) = output(i) + err;
                                    }
                                }
                              
                              if(Name.compare("Magnetorquer") == 0)
                                {
                                Vec3d m, mSC, Torque;
                                B_Body = B_Body*1E-9; // Conversion from nanotesla to tesla
                                //cout << "ontime: " << ontime << endl;
                                m = ontime2dipole*ontime; // Magnetic dipole vector produced by the magnetic torquers (magnetic torquers frame)
                                
                                for(int i = 0; i < 3; i++)
                                    {
                                    err = Error1(generator1); // Measurement error
                                    if(m(i) != 0.0) m(i) = m(i) + err;
                                    }
                                
                                mSC = (SC2SYS.transpose())*m; // Transform vector m from magnetic torquers frame to SC body-fixed frame ( inverse of SC2SYS = transpose of SC2SYS )
                                Torque = mSC.cross(B_Body);
                                
                                // ontime = Vec3d::Zero();
                                
                                output = Torque;
                                }
                              }
                                
                            return(output);     
                            };
    

}; // End of namespace spaceenvironment



