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
#include <omp.h>

#include <Gravity.h>
#include <Constants.h>
#include <IO_utils.h>

using namespace std;
using namespace math;
using namespace constants;

namespace gravity
    {
    //------------------------------------------------------------------------------
    // GRAV implementation
    //------------------------------------------------------------------------------
    MatnMAXxnMAXd GRAV::C = MatnMAXxnMAXd::Zero();
    MatnMAXxnMAXd GRAV::S = MatnMAXxnMAXd::Zero();
    MatnMAXxnMAXd GRAV::sigmaC = MatnMAXxnMAXd::Zero();
    MatnMAXxnMAXd GRAV::sigmaS = MatnMAXxnMAXd::Zero();
    MatnMAXxnMAXd GRAV::P = MatnMAXxnMAXd::Zero();
    MatnMAXxnMAXd GRAV::Pd = MatnMAXxnMAXd::Zero();
    MatnMAXxnMAXd GRAV::A = MatnMAXxnMAXd::Zero();
    MatnMAXxnMAXd GRAV::B = MatnMAXxnMAXd::Zero();
    MatnMAXxnMAXd GRAV::F = MatnMAXxnMAXd::Zero();
    MatnMAXxnMAXd GRAV::CS = MatnMAXxnMAXd::Zero();
    MatnMAXxnMAXd GRAV::CS1 = MatnMAXxnMAXd::Zero();
    VectornMAXd GRAV::Sigma1 = VectornMAXd::Zero();
    VectornMAXd GRAV::sin_mlambda = VectornMAXd::Zero();
    VectornMAXd GRAV::cos_mlambda = VectornMAXd::Zero();
    double GRAV::mu = 0.0;
    double GRAV::R = 0.0;
    //------------------------------------------------------------------------------
    // Destructor
    //------------------------------------------------------------------------------
    GRAV::~GRAV() {};
    //------------------------------------------------------------------------------
    // Method double Legendre_norm(int n, int m, double x)
    //------------------------------------------------------------------------------
    /**
     * Normalized associated Legendre functions
     *
     * @param n	Degree
     * @param m	Order
     * @param x	Function variable
     *
     * @return		Return the value of the normalized associated Legendre function of degree n and order m corresponding to x
     *
     */
    //------------------------------------------------------------------------------
    double GRAV::Legendre_norm(unsigned int n,
                                unsigned int m,
                                double x)
                            {
                            double delta0m, Pnm, Pnm_norm;
                            
                            if(m > n) Pnm_norm = 0.0;
                            else
                                            {
                                            // Note: For applications in geodesy and related fields it is important to use the definition of the Legendre associated functions
                                            // WITHOUT the factor (-l)^m in order to be consistent with published geopotential coefficients. The term pow(-1.0, -m) is used because
                                            // assoc_legendre uses the definition WITH the factor (-l)^m and (-l)^m * (-l)^-m = 1
                                            //Pnm = pow(-1.0, -m)*boost::math::legendre_p(n, m, x);
                                            Pnm = pow(-1.0, -m)*assoc_legendre(n, m, x);
                                            
                                            delta0m = (m == 0) ? 1.0 : 0;
                                            
                                            //cout << "n = " << n << " m = " << m <<" n - m = " << n - m << " n + m = " << n + m << endl;
                                            
                                            Pnm_norm = Pnm*sqrt( (2.0 - delta0m)*(2.0*n + 1.0)*boost::math::factorial<double>(n - m)/boost::math::factorial<double>(n + m) );
                                            }
                            
                            return(Pnm_norm);
                            };
    //------------------------------------------------------------------------------
    // Method double d_Legendre_norm(int n, int m, double x)
    //------------------------------------------------------------------------------
    /**
     * First derivative of normalized associated Legendre functions
     *
     * @param n	Degree
     * @param m	Order
     * @param x	Function variable
     *
     * @return		Return the value of the first derivative of the normalized associated Legendre function of degree n and order m corresponding to x
     *
     */
    //------------------------------------------------------------------------------
    double GRAV::d_Legendre_norm(unsigned int n,
                                unsigned int m,
                                double x)
                                {
                                double delta0m, d_Pnm_norm;
                                                                                                        
                                delta0m = (m == 0) ? 1.0 : 0;
                                
                                d_Pnm_norm = sqrt( (n + m + 1.0)*(n - m)/( (1.0 + delta0m)*(1.0 - x*x) ) )*Legendre_norm(n,m+1,x) - ( m*x/(1.0 - x*x) )*Legendre_norm(n,m,x);
                                
                                //cout << "d_Pnm_norm: " << d_Pnm_norm << endl;
                                
                                return(d_Pnm_norm);
                                };
    //------------------------------------------------------------------------------
    // Method void Legendre_norm_FC(double u, double w)
    //------------------------------------------------------------------------------
    /**
     * Normalized associated Legendre functions and first derivatives (forward columns recursive method)
     *
     * @param u		cos(phi) = sin(theta) where theta + phi  = pi/2.0 and phi = geocentric latitude
     * @param u		sin(phi) = cos(theta)
     *
     * @return		Return the matrices whose elements are the value of the normalized associated Legendre function of degree n and order m and its first derivative 
     * 			corresponding to u = cos(theta). The computation of the Legendre functions is based on a forward columns (FC) recursive method
     */
    //------------------------------------------------------------------------------
    void GRAV::Legendre_norm_FC(double u,
                                double w)
                                {
                                double a = w/u;
                                
                                P(0,0) = 1.0;
                                P(1,1) = sqrt(3.0)*u;
                                // Compute P(1,0) with the forward raw recursive relation
                                P(1,0) = w*sqrt(3)*P(0,0);
                                
                                //double prod_pi;
                                //Compute first column (m = 0) and diagonal elements of P
                                for(int n = 2; n <= n_max; n++)
                                    {
                                    P(n,0) = A(n,0)*w*P(n-1,0) - B(n,0)*P(n-2,0);
                                    Pd(n,0) = ( n*w*P(n,0) - F(n,0)*P(n-1,0) )/u;
                                    
                                    //Diagonal elements
                                    P(n,n) = u*sqrt( (2.0*n + 1.0)/(2.0*n) )*P(n-1,n-1);
                                    //prod_pi = 1.0;
                                    //for(int i = 2; i <= n; i++) prod_pi = prod_pi*sqrt( (2.0*i + 1.0)/(2.0*i) );
                                    //P(n,n) = pow(u,n)*sqrt(3.0)*prod_pi;
                                    Pd(n,n) = n*a*P(n,n);
                                    }
                                
                                //Compute off-diagonal terms
                                #pragma omp parallel for collapse(1)
                                for(int m = 1; m <= n_max; m++)
                                    {
                                    for(int n = m + 1; n <= n_max; n++)
                                        {
                                        P(n,m) = A(n,m)*w*P(n-1,m) - B(n,m)*P(n-2,m);
                                        Pd(n,m) = ( n*w*P(n,m) - F(n,m)*P(n-1,m) )/u;
                                        }
                                    }
                                
                                };
    //------------------------------------------------------------------------------
    // Method Vec3d field_vec(double time, const Ref<const VectorXd>& orbstate)
    //------------------------------------------------------------------------------
    /**
     * Implementation of the abstract method of class SPACEENV
     * Compute the acceleration vector due to the Earth's gravity field
     *
     * @param time        GPS epoch (seconds) of the input state
     * @param orbstate    Spacecraft position vector (ECEF)
     *
     * @return 3-dimensional gravity acceleration vector (ECI)
     *
     */
    //------------------------------------------------------------------------------
    Vec3d GRAV::field_vec(double time,
                          const Ref<const VectorXd>& orbstate)
                            {
                            Vec3d acc, acc_ECEF;
                            
                            Vec3d SC_posECEF = orbstate;
                            
                            double x = SC_posECEF(0);
                            double y = SC_posECEF(1);
                            double z = SC_posECEF(2);
                            
                            //cout << "time: " << time << endl;
                            //lonlath = ECEF2lonlath(SC_posECEF);
                            //lambda = lonlath(0);
                            //phi = lonlath(1);
                            double lambda = atan2(y,x);
                            double phi = atan2(z,sqrt(x*x + y*y));
                            double r = sqrt(x*x + y*y + z*z);
																												
							double sinphi = sin(phi);
                            double cosphi = cos(phi);
                            double sinlambda = sin(lambda);
                            double coslambda = cos(lambda);
                            
                            double u = cosphi; //sin(theta)
                            double w = sinphi; //cos(theta)
                            
                            double dUdr, dUdphi, dUdlambda, dUdx, dUdy, dUdz;
                            
                            double drdx = cosphi*coslambda;
                            double drdy = cosphi*sinlambda;
                            double drdz = sinphi;
                            
                            double dphidx = -(1/r)*sinphi*coslambda;
                            double dphidy = -(1/r)*sinphi*sinlambda;
                            double dphidz = (1/r)*cosphi;
                            
                            double dlambdadx = -sinlambda/(r*cosphi);
                            double dlambdady = coslambda/(r*cosphi);
                            double dlambdadz = 0.0;
                            
                            double q = 0.0;
                            double n1 = 0.0;
                            
                            Legendre_norm_FC(u,w);
                            
                            // Tabularization of sin(m*lambda) and cos(m*lambda)
                            #pragma omp parallel for collapse(1)
                            for(int m = 0; m <= n_max; m++)
                                {
                                sin_mlambda(m) = sin(m*lambda);
                                cos_mlambda(m) = cos(m*lambda);
                                }
                            // Tabularization of q*( C(n,m)*cos_mlambda(m) + S(n,m)*sin_mlambda(m) ) and q*( -C(n,m)*sin_mlambda(m) + S(n,m)*cos_mlambda(m) ) with internal loop unrolling
                            #pragma omp parallel for collapse(1)
                            for(int n = 2; n <= n_max; n++)
                                {
                                q = pow(R/r, n);
                                
                                for(int m = 0; m <= n; m++)
                                    {
                                    CS(n,m) = q*( C(n,m)*cos_mlambda(m) + S(n,m)*sin_mlambda(m) );
                                    CS1(n,m) = q*( -C(n,m)*sin_mlambda(m) + S(n,m)*cos_mlambda(m) );
                                    }
                                }
                            
                            // Partial derivative of potential U with respect to r
                            #pragma omp parallel for collapse(1)
                            for(int n = 2; n <= n_max; n++)
                                {
                                VectorXd Sigma2 = VectorXd::Zero(n_max+1);
                                
                                n1 = n + 1.0;
                                
                                for(int m = 0; m <= n; m++)
                                    {
                                    ////double Pnm = Legendre_norm(n,m,sinphi);//P(n,m);//
                                    //Sigma2(m) = (n + 1.0)*q*P(n,m)*( C(n,m)*cos(m*lambda) + S(n,m)*sin(m*lambda) );
                                    Sigma2(m) = n1*P(n,m)*CS(n,m);
                                    }
                                    
                                Sigma1(n) = Sigma2.sum();
                                }
                            
                            dUdr = -( mu/(r*r) )*( 1.0 + Sigma1.sum());
                            
                            Sigma1 = VectornMAXd::Zero();
                            // Partial derivative of potential U with respect to phi
                            #pragma omp parallel for collapse(1)
                            for(int n = 2; n <= n_max; n++)
                                {
                                VectorXd Sigma2 = VectorXd::Zero(n_max+1);
                                
                                for(int m = 0; m <= n; m++)
                                    {
                                    ////double dPnmdphi = cosphi*d_Legendre_norm(n,m,sinphi);//-Pd(n,m);//cosphi*Ld(n,m);//
                                    //double dPnmdphi = -Pd(n,m);
                                    //Sigma2(m) = q*(-Pd(n,m))*( C(n,m)*cos(m*lambda) + S(n,m)*sin(m*lambda) );
                                    Sigma2(m) = -Pd(n,m)*CS(n,m);
                                    }
                                    
                                Sigma1(n) = Sigma2.sum();
                                }
                            
                            dUdphi = (mu/r)*Sigma1.sum();
                            //dUdphi = 0.0;
                            
                            Sigma1 = VectornMAXd::Zero();
                            // Partial derivative of potential U with respect to lambda
                            #pragma omp parallel for collapse(1)
                            for(int n = 2; n <= n_max; n++)
                                {
                                VectorXd Sigma2 = VectorXd::Zero(n_max+1);
                                
                                for(int m = 0; m <= n; m++)
                                    {
                                    ////double Pnm = Legendre_norm(n,m,sinphi);//P(n,m);//
                                    //Sigma2(m) = m*q*P(n,m)*( -C(n,m)*sin(m*lambda) + S(n,m)*cos(m*lambda) );
                                    Sigma2(m) = m*P(n,m)*CS1(n,m);
                                    }
                                    
                                Sigma1(n) = Sigma2.sum();
                                }
                            
                            dUdlambda = (mu/r)*Sigma1.sum();
                            
                            //cout << "dUdr: " << dUdr << "  dUdphi: " << dUdphi << "  dUdlambda: " << dUdlambda << endl;
                                    
                            // Compute accelerations
                            
                            dUdx = dUdr*drdx + dUdphi*dphidx + dUdlambda*dlambdadx;
                            
                            dUdy = dUdr*drdy + dUdphi*dphidy + dUdlambda*dlambdady;
                            
                            dUdz = dUdr*drdz + dUdphi*dphidz + dUdlambda*dlambdadz;
                            
                            acc_ECEF(0) = dUdx;
                            acc_ECEF(1) = dUdy;    
                            acc_ECEF(2) = dUdz;
                            
                            //Vec3d r_vec;
                            //r_vec << x, y, z;
                            //cout << "phi_acc: " << acc_ECEF.dot(r_vec)/(acc_ECEF.norm()*r_vec.norm()) << endl;
                            //cout << "acc_ECEF(0): " << acc_ECEF(0) << "  acc_ECEF(1): " << acc_ECEF(1) << "  acc_ECEF(2): " << acc_ECEF(2) << endl;
                            //cout << "acc_ECEF: " << acc_ECEF.norm() << endl;
            
                            acc = v3D_transform(time, acc_ECEF, "ITRF93", "J2000"); // Acceleration vector in ECI frame
                            //if(refsys.compare("ECEF") == 0) acc = acc_ECEF;
                            //cout << "acc(0): " << acc(0) << "  acc(1): " << acc(1) << "  acc(2): " << acc(2) << endl;
                            
                            return(acc);
                            };
    //------------------------------------------------------------------------------
    // Method getmodel_coeff()
    //------------------------------------------------------------------------------
    /**
     * Get surface optical coefficients
     *
     * @return 3-dimensional vector containing in order the specular reflectivity,
     * diffuse reflectivity and transmitted portions of incoming photons coefficients
     * currently in use
     */
    //------------------------------------------------------------------------------
    void GRAV::getmodel_coeff()
                            {
                            Matrix3D gravmodel_coeff;
                            gravmodel_coeff.resize(boost::extents[n_max+1][n_max+1][4]);
                            string modelfile = modelfilepath + "/gravityfield/" + modelname + ".gfc";
                            
                            gravmodel_coeff = read_gfc(modelfile.c_str(), n_max, grav_epoch, mu, R); // Function from IO_utils.h
                            
                            int n, m;
                            Matrix3D_index N = 0, M = 0;
                            
                            #pragma omp parallel for collapse(2)
                            for(n = 0; n <= n_max; n++)
                                for(m = 0; m <= n_max; m++)
                                    {
                                    N = n; M = m;
                                    C(n,m) = gravmodel_coeff[N][M][0];
                                    S(n,m) = gravmodel_coeff[N][M][1];
                                    sigmaC(n,m) = gravmodel_coeff[N][M][2];
                                    sigmaS(n,m) = gravmodel_coeff[N][M][3];
                                    
                                    A(n,m) = sqrt( (2.0*n - 1.0)*(2.0*n + 1.0)/( (n - m)*(n + m) ) );
                                    B(n,m) = sqrt( (2.0*n + 1.0)*(n + m - 1.0)*(n - m - 1.0)/( (n - m)*(n + m)*(2.0*n - 3.0) ) );
                                    F(n,m) = sqrt( (n*n - m*m)*(2.0*n + 1.0)/(2.0*n - 1.0) );
                                    }
                            };

}; // End of namespace gravity