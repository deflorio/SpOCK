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
#include <fstream>
#include <string>
#include <math.h>

#include <Events.h>
#include <Transformations.h>
#include <AnalyticalModels.h>
#include <Constants.h>
#include <Eigen/Core>
#include <Eigen/Geometry> 
#include <VarTypes.h>

extern "C"
      {
      #include <SpiceUsr.h>
      }

using namespace std;
using namespace math;
using namespace Eigen;
using namespace constants;
using namespace mathconst;
using namespace astro;


//------------------------------------------------------------------------------
// Matrix<double,Dynamic,6> GSsContacts(int SC_num, int PL_num, Vec4d& orbel, const MatrixXd& orbstateECEF, const MatrixXd& targets, string output_path, double FOV, double FOVdir, double duration, int simstep);
//------------------------------------------------------------------------------
/**
 * Compute sensor contacts of one or more spacecraft with one or more targets
 *
 * @param SC_num          Number of spacecraft considered
 * @param PL_num          Number of orbital plane considered
 * @param orbel           4-D vector containing nominal mean orbital elements a, ex, ey and i
 * @param orbstateECEF    7-D vector containing epoch (0) and ECEF orbital state (1-6) at epoch 
 * @param TGs_grid_lons   Vector containing all longitudes of target locations of the targets grid [deg]
 * @param TGs_grid_lats   Vector containing all latitudes of target locations of the targets grid [deg]
 * @param TGs_list        Array of type ground::TG containing name, longitude, latitude and altitude of targets
 * @param output_path     Location of output files containing contacts information    
 * @param FOV             Field of view half angle (perpendicolar to the direction of motion) [deg]
 * @param FOVdir          Field of view half angle (along the direction of motion) [deg]
 * @param comp_duration        Time span of analysis [s]
 * @param simstep         Step of input ephemerides [s]
 *
 * @return .csv file containing contact information for all spacecraft (AllContacts_file)
 * @return .csv file containing contact information for the spacecraft considered (written inside the function) 
 */
//------------------------------------------------------------------------------                      
void TGsContacts( int SC_num,
                  int PL_num,
                  Vec4d& orbel,
                  const MatrixXd& orbstateECEF,
                  VectorXd& TGs_grid_lons,
                  VectorXd& TGs_grid_lats,
                  ground::TG* TGs_list,
                  string output_path,
                  double FOV,
                  double FOVdir,
                  double comp_duration,
                  int simstep,
                  bool TGs_on,
                  bool TGs_grid_on,
                  ofstream& AllContacts_file)
                  {
                  int last_element = comp_duration/simstep;
                  
                  bool maxel_computed = false;
                  
                  //////////////////////////// Put in spacecraft structs //////////////////////////
                  
                  double a_nom, e, inc, incdeg, h_sat;
                  double sensFOV, sensFOVdir, eta, alpha, Ftrans, F, Fdir;
                  double aF, bF;
                  double n, dOMdt, domdt, nd, ki;
                  double i_app, delta_lon, delta_lat;
                  
                  double max_elevation, time_maxel;
                  //double SATlon, SATlat, deltaL, lambda;
                  //double sinrho, sin_lambda, cos_lambda, epsilon_first, epsilon_last;
                  double Az, Az_in, Az_out, Az_maxel, El, El_first, El_last, El_in, El_out;
                  Vec3d AzElAlt;
                  Vector6d stateECEF;
                  VectorNd<7> ephem_row;
                  Vec3d posECEF, lonlath, lonlatrow;
                  
                  // Input orbital elements
                  a_nom = orbel(0);
                  e = sqrt(orbel(1)*orbel(1) + orbel(2)*orbel(2));
                  inc = orbel(3);
                  incdeg = inc*RAD2DEG;
          
                  h_sat = a_nom - R_EARTH;
                  
                  // Ground swath (lateral)
                  sensFOV = FOV*DEG2RAD;
                  eta = (R_EARTH + h_sat)/R_EARTH;
                  alpha = -sensFOV + asin(eta*sin(sensFOV));
                  Ftrans = alpha;
                  F = R_EARTH*alpha;
                  // Ground swath (along the motion)
                  sensFOVdir = FOVdir*DEG2RAD;
                  eta = (R_EARTH + h_sat)/R_EARTH;
                  alpha = -sensFOVdir + asin(eta*sin(sensFOVdir));
                  Fdir = alpha;
                  
                  // Elliptical footprint parameters
                  aF = Ftrans*RAD2DEG; // Semi-major axis of payload footprint
                  bF = Fdir*RAD2DEG; // Semi-minor axis of payload footprint
                  
                  // Apparent inclination
                  Vector6d J4mod_output = J4model(a_nom,e,inc);
                  
                  n = J4mod_output(0);
                  dOMdt = J4mod_output(2);
                  domdt = J4mod_output(3);
                  
                  nd = n + domdt;
                  ki = nd/(OMEGA_EARTH - dOMdt); // Daily recurrence frequency of the orbit
                  
                  i_app = atan( n*sin(inc)/( n*cos(inc) - (OMEGA_EARTH - dOMdt) ) );
                  i_app = mod(i_app,PI);
          
                  delta_lon = fabs(F*sin(i_app));
                  delta_lat = fabs(F*cos(i_app));
          
                  delta_lon = (delta_lon/R_EARTH)*RAD2DEG;
                  delta_lat = (delta_lat/R_EARTH)*RAD2DEG;
                  
                  // Longitude and latitude
                  int matrows = orbstateECEF.rows();
                  //int len = lonlat.rows();
                  int len = matrows;
                  
                  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
                  //////////////////////////////////////////// CONTACTS WITH TARGETS //////////////////////////////////////////
                  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
                  
                  string TG_filename = output_path + "/S" + to_string(SC_num) + "-P" + to_string(PL_num) + "_TG_contacts.csv";
                  //cout << TG_filename << endl;
                  ofstream Contacts_file;
                  Contacts_file.open(TG_filename);
                  
                  Contacts_file << "TG,Epoch UTC in,GPS time [s] in,Duration [s],El. in,El. out,Epoch UTC Max. El.,Max. El. [deg],Az. in [deg],Az. out [deg],Az. Max El. [deg],lon [deg],lat [deg]" << endl;
          
                  int grid_rows = TGs_grid_lons.size();
                  int grid_cols = TGs_grid_lats.size();
                  
                  SpiceChar time_in_string[35];
                  SpiceChar time_maxel_string[35];
                  //SpiceChar time_out_string[35];
                  ConstSpiceChar* format;
                  //format = "DD-MM-YYYY HR:MN:SC";
                  format = "YYYY-MM-DD HR:MN:SC";
          
                  if(TGs_grid_on)
                    {
                    for(int row = 0 ; row < grid_rows; row++)
                      {
                      for(int col = 0 ; col < grid_cols; col++)
                        {
                        //f(h) = 0;
                        
                        double TGlon = TGs_grid_lons(row);
                        double TGlat = TGs_grid_lats(col);
                        
                        string TGname = "TG_" + to_string(row) + "-" + to_string(col);
                        
                        double time, time_in, time_out, sec_J2000, x0, y0, xT, yT, duration;
                        VectorNd<2> xyT = VectorNd<2>::Zero();
                        
                        int i = 0;
                        int* i_ptr = &i;
                        int ind;
                        
                        while(i < last_element)
                            {
                            ephem_row = orbstateECEF.row(i);
                            time = ephem_row(0);
                            stateECEF = ephem_row.segment(1,6);
                            
                            posECEF = ephem_row.segment(1,3);                                        
                            lonlath = ECEF2lonlath(posECEF);
                            
                            x0 = mod(lonlath(0),PI2)*RAD2DEG;
                            y0 = lonlath(1)*RAD2DEG;
                            
                            AzElAlt = ECEF2AzElAlt(stateECEF, TGlon*DEG2RAD, TGlat*DEG2RAD);
                            
                            Az = AzElAlt(0);
                            El = AzElAlt(1);
                            
                            Az = Az*RAD2DEG;
                            El = El*RAD2DEG;
                            
                            xyT  = lonlat2satcart(x0, y0, TGlon, TGlat, incdeg, ki);
                            
                            xT = xyT(0);
                            yT = xyT(1);
                            
                            double theta = mod(atan2(yT,xT)*RAD2DEG,PI2);
                            if(isnan(theta)) cout << "Theta is a NaN" << endl;
                            
                            if( (xT*xT)*(bF*bF) + (yT*yT)*(aF*aF) <= (aF*aF)*(bF*bF) ) // Condition of inclusion of target-point into the elliptical footprint of the payload
                              {
                              time_in = time;
                              sec_J2000 = GPS2ET(time_in);
                              
                              timout_c( sec_J2000, format, 35, time_in_string );
                              
                              El_first = El;
                              El_in = El;
                              Az_in = Az;
          
                              ind = i + 1;
          
                              do
                                {
                                ephem_row = orbstateECEF.row(ind);
                                time = ephem_row(0);
                                stateECEF = ephem_row.segment(1,6);
                                
                                posECEF = ephem_row.segment(1,3);                                        
                                lonlath = ECEF2lonlath(posECEF);
                                
                                AzElAlt = ECEF2AzElAlt(stateECEF, TGlon*DEG2RAD, TGlat*DEG2RAD);
                                
                                Az = AzElAlt(0);
                                El = AzElAlt(1);
                                
                                Az = Az*RAD2DEG;
                                El = El*RAD2DEG;
                                
                                El_last = El;
                                
                                x0 = mod(lonlath(0),PI2)*RAD2DEG;
                                y0 = lonlath(1)*RAD2DEG;
        
                                xyT  = lonlat2satcart(x0, y0, TGlon, TGlat, incdeg, ki);
                            
                                xT = xyT(0);
                                yT = xyT(1);
                                
                                if( (El_last <= El_first) && !maxel_computed)
                                    {
                                    Az_maxel = Az;
                                    max_elevation = El_first;
                                    
                                    time_maxel = time;
                                    sec_J2000 = GPS2ET(time_maxel);
                                    timout_c( sec_J2000, format, 35, time_maxel_string );
                                    
                                    maxel_computed = true;
                                    }
                              
                                El_first = El_last;
                                
                                ind++;
                                
                                }while( ( (xT*xT)*(bF*bF) + (yT*yT)*(aF*aF) <= (aF*aF)*(bF*bF) ) && ind < len );//while( ( (xT/aF)*(xT/aF) + (yT/bF)*(yT/bF) <= 1.0 ) && ind < len );
                        
                              time_out = time;
                              duration = time_out - time_in;
                              maxel_computed = false;
                              
                              Az_out = Az;
                              El_out = El;
                              
                              Contacts_file << TGname << "," << string(time_in_string) << "," << time_in << "," << duration << "," << El_in << "," << El_out << "," << string(time_maxel_string) << "," << max_elevation << "," << Az_in << "," << Az_out << "," << Az_maxel << "," << TGlon << "," << TGlat << "," << endl;
                              
                              AllContacts_file << TGname << "," << string(time_in_string) << "," << time_in << "," << duration << "," << El_in << "," << El_out << "," << string(time_maxel_string) << "," << max_elevation << "," << Az_in << "," << Az_out << "," << Az_maxel << "," << TGlon << "," << TGlat << "," << PL_num << "," << SC_num << endl;
                              
                              *i_ptr = ind;
                              }
                              
                              
                            i++;
                            }
                        }
                      }
                    }
                    
                  if(TGs_on)
                    {
                    unsigned int TGind = 0;
                    unsigned int listsize = 1000;
                    
                    while( !TGs_list[TGind].name.empty() && TGind < listsize )
                        {
                        string TGname = TGs_list[TGind].name;//"TG" + to_string(TGind);
                    
                        double TGlon = TGs_list[TGind].lon;//targets(TGind,0);
                        double TGlat = TGs_list[TGind].lat;//targets(TGind,1);
                      
                        double time, time_in, time_out, sec_J2000, x0, y0, xT, yT, duration;
                        VectorNd<2> xyT = VectorNd<2>::Zero();
                      
                        int i = 0;
                        int* i_ptr = &i;
                        int ind;
                      
                        while(i < last_element)
                            {
                            ephem_row = orbstateECEF.row(i);
                            time = ephem_row(0);
                            stateECEF = ephem_row.segment(1,6);
                            
                            posECEF = ephem_row.segment(1,3);                                        
                            lonlath = ECEF2lonlath(posECEF);
                            
                            x0 = mod(lonlath(0),PI2)*RAD2DEG;
                            y0 = lonlath(1)*RAD2DEG;
                            
                            xyT  = lonlat2satcart(x0, y0, TGlon, TGlat, incdeg, ki);
                            
                            xT = xyT(0);
                            yT = xyT(1);
                            
                            double theta = mod(atan2(yT,xT)*RAD2DEG,PI2);
                            if(isnan(theta)) cout << "Theta is a NaN" << endl;
                            
                            AzElAlt = ECEF2AzElAlt(stateECEF, TGlon*DEG2RAD, TGlat*DEG2RAD);
                            
                            Az = AzElAlt(0);
                            El = AzElAlt(1);
                            
                            Az = Az*RAD2DEG;
                            El = El*RAD2DEG;
                            
                            if( (xT*xT)*(bF*bF) + (yT*yT)*(aF*aF) <= (aF*aF)*(bF*bF) ) // Condition of inclusion of target-point into the elliptical footprint of the payload
                              {
                              time_in = time;
                              sec_J2000 = GPS2ET(time_in);
                              
                              timout_c( sec_J2000, format, 35, time_in_string );
                              
                              El_first = El;
                              El_in = El;
                              Az_in = Az;
                              
                              ind = i + 1;
          
                              do
                                {
                                ephem_row = orbstateECEF.row(ind);
                                time = ephem_row(0);
                                stateECEF = ephem_row.segment(1,6);
                                
                                posECEF = ephem_row.segment(1,3);                                        
                                lonlath = ECEF2lonlath(posECEF);
                                
                                AzElAlt = ECEF2AzElAlt(stateECEF, TGlon*DEG2RAD, TGlat*DEG2RAD);
                                
                                Az = AzElAlt(0);
                                El = AzElAlt(1);
                                
                                Az = Az*RAD2DEG;
                                El = El*RAD2DEG;
                                
                                El_last = El;
                                
                                x0 = mod(lonlath(0),PI2)*RAD2DEG;
                                y0 = lonlath(1)*RAD2DEG;
        
                                xyT  = lonlat2satcart(x0, y0, TGlon, TGlat, incdeg, ki);
                            
                                xT = xyT(0);
                                yT = xyT(1);
                                
                                if( (El_last <= El_first) && !maxel_computed)
                                    {
                                    Az_maxel = Az;
                                    max_elevation = El_first;
                                    
                                    time_maxel = time;
                                    sec_J2000 = GPS2ET(time_maxel);
                                    timout_c( sec_J2000, format, 35, time_maxel_string );
                                    
                                    maxel_computed = true;
                                    }
                              
                                El_first = El_last;
                                
                                ind++;//ind = ind + 1;
                                }while( ( (xT*xT)*(bF*bF) + (yT*yT)*(aF*aF) <= (aF*aF)*(bF*bF) ) && ind < len );//while( ( (xT/aF)*(xT/aF) + (yT/bF)*(yT/bF) <= 1.0 ) && ind < len );
                              
                              time_out = time;
                              duration = time_out - time_in;
                              maxel_computed = false;
                              
                              Az_out = Az;
                              El_out = El;
                              
                              Contacts_file << TGname << "," << string(time_in_string) << "," << time_in << "," << duration << "," << El_in << "," << El_out << "," << string(time_maxel_string) << "," << max_elevation << "," << Az_in << "," << Az_out << "," << Az_maxel << "," << TGlon << "," << TGlat << "," << endl;
                              
                              AllContacts_file << TGname << "," << string(time_in_string) << "," << time_in << "," << duration << "," << El_in << "," << El_out << "," << string(time_maxel_string) << "," << max_elevation << "," << Az_in << "," << Az_out << "," << Az_maxel << "," << TGlon << "," << TGlat << "," << PL_num << "," << SC_num << endl;
                              
                              *i_ptr = ind;
                              }
                              
                            i++;
                            }
                        TGind++;
                        }  
                    }
          
                  Contacts_file.close();
                  };                  
//------------------------------------------------------------------------------
// GSsContacts(int SC_num, int PL_num, Vec4d& orbel, const MatrixXd& orbstateECEF, const MatrixXd& groundstations, string output_path, double min_elevation, double duration, int simstep);
//------------------------------------------------------------------------------
/**
 * Compute sensor contacts of one or more spacecraft with one or more targets
 *
 * @param SC_num          Number of spacecraft considered
 * @param PL_num          Number of orbital plane considered
 * @param orbel           4-D vector containing nominal mean orbital elements a, ex, ey and i
 * @param orbstateECEF    7-D vector containing epoch (0) and ECEF orbital state (1-6) at epoch 
 * @param GSs_list        Array of type ground::GS containing name, longitude, latitude, altitude and minimum spacecraft elevation of ground stations
 * @param output_path     Location of output files containing contacts information    
 * @param comp_duration        Time span of analysis [s]
 * @param simstep         Step of input ephemerides [s]
 *
 * @return .csv file containing contact information for all spacecraft (AllContacts_file)
 * @return .csv file containing contact information for the spacecraft considered (written inside the function)        
 */
//------------------------------------------------------------------------------                      
void GSsContacts( int SC_num,
                  int PL_num,
                  Vec4d& orbel,
                  const MatrixXd& orbstateECEF,
                  ground::GS* GSs_list,
                  string output_path,
                  double comp_duration,
                  int simstep,
                  ofstream& AllContacts_file)
                  {
                  int last_element = comp_duration/simstep;
                  
                  bool maxel_computed = false;
                  
                  int matrows = orbstateECEF.rows();
                  MatrixXd lonlat = orbstateECEF.block(0,0,matrows,3);
                  
                  int len = matrows;
          
                  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
                  /////////////////////////////////////// CONTACTS WITH GROUND STATIONS ///////////////////////////////////////
                  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
                  
                  string GS_filename = output_path + "/S" + to_string(SC_num) + "-P" + to_string(PL_num) + "_GS_contacts.csv";
                  //cout << GS_filename << endl;
                  ofstream Contacts_file;
                  Contacts_file.open(GS_filename);
                  
                  Contacts_file << "GS,AOS UTC,LOS UTC,AOS [GPS secs],Duration [m],Epoch UTC Max. El.,Max Elevation [deg],Az. AOS [deg],Az. LOS [deg],Az. Max El. [deg], lon [deg],lat [deg]" << endl;
                  
                  SpiceChar time_in_string[35];
                  SpiceChar time_out_string[35];
                  SpiceChar time_maxel_string[35];
                  ConstSpiceChar* format;
                  //format = "DD-MM-YYYY HR:MN:SC";
                  format = "YYYY-MM-DDTHR:MN:SCZ";
                  
                  string GSname;
                  double GSlon, GSlat, min_elevation, max_elevation;
                  double time, time_in, time_out, time_maxel, sec_J2000, duration;
                  double Az, Az_in, Az_out, Az_maxel, El, El_first, El_last;
                  Vec3d AzElAlt;
                  Vector6d stateECEF;
                  VectorNd<7> ephem_row;
                  
                  unsigned int GSind = 0;
                  unsigned int listsize = 1000;
                  
                  while( !GSs_list[GSind].name.empty() && GSind < listsize )
                        {
                        GSname = GSs_list[GSind].name;//"GS" + to_string(GSind);
                        
                        GSlon = GSs_list[GSind].lon;//targets(GSind,0);
                        GSlat = GSs_list[GSind].lat;//targets(GSind,1);
                        min_elevation = GSs_list[GSind].minelev;
                        
                        int i = 0;
                        int* i_ptr = &i;
                        int ind;
                        
                        while(i < last_element)
                                {
                                ephem_row = orbstateECEF.row(i);
                                time = ephem_row(0);
                                stateECEF = ephem_row.segment(1,6);
                                
                                AzElAlt = ECEF2AzElAlt(stateECEF, GSlon*DEG2RAD, GSlat*DEG2RAD);
                                
                                Az = AzElAlt(0);
                                El = AzElAlt(1);
                                
                                Az = Az*RAD2DEG;
                                El = El*RAD2DEG;
                            
                                time_out = 0.0;
                                  
                                if( El >= min_elevation )
                                  {
                                  time_in = time;
                                  sec_J2000 = GPS2ET(time_in);
                                  
                                  timout_c( sec_J2000, format, 35, time_in_string );
                                  
                                  El_first = El;
                                  Az_in = Az;
                                  
                                  ind = i + 1;
              
                                  do
                                    {
                                    ephem_row = orbstateECEF.row(ind);
                                    time = ephem_row(0);
                                    stateECEF = ephem_row.segment(1,6);
                                    
                                    AzElAlt = ECEF2AzElAlt(stateECEF, GSlon*DEG2RAD, GSlat*DEG2RAD);
                                    
                                    Az = AzElAlt(0);
                                    El = AzElAlt(1);
                                    
                                    Az = Az*RAD2DEG;
                                    El = El*RAD2DEG;
                                    
                                    El_last = El;
                                
                                    if( (El_last <= El_first) && !maxel_computed)
                                          {
                                          Az_maxel = Az;
                                          max_elevation = El_first;
                                          
                                          time_maxel = time;
                                          sec_J2000 = GPS2ET(time_maxel);
                                          timout_c( sec_J2000, format, 35, time_maxel_string );
                                          
                                          maxel_computed = true;
                                          }
                                    
                                    if(El_first >= min_elevation && El_last <= min_elevation && min_elevation != 0.0)
                                          {
                                          time_out = time;
                                          Az_out = Az;      
                                          }
                                    
                                    El_first = El_last;
                                    
                                    ind++;
                                    
                                    }while( ( El >= min_elevation ) && ind < len );
                                  
                                  if(time_out == 0.0)
                                    {
                                    time_out = time;
                                    Az_out = Az;
                                    }
                                  
                                  sec_J2000 = GPS2ET(time_out);
                                  timout_c( sec_J2000, format, 35, time_out_string );
                                  
                                  duration = time_out - time_in;
                                  maxel_computed = false;
                                  
                                  Contacts_file << GSname << "," << string(time_in_string) << "," << string(time_out_string) << "," << fixed << time_in << "," << duration/60.0 << "," << string(time_maxel_string) << "," << max_elevation << "," << Az_in << "," << Az_out << "," << Az_maxel << "," << GSlon << "," << GSlat << "," << endl;
                                  
                                  AllContacts_file << GSname << "," << string(time_in_string) << "," << string(time_out_string) << "," << fixed << time_in << "," << duration/60.0 << "," << string(time_maxel_string) << "," << max_elevation << "," << Az_in << "," << Az_out << "," << Az_maxel << "," << GSlon << "," << GSlat << "," << PL_num << "," << SC_num << endl;
                                  
                                  *i_ptr = ind;
                                  }
                                  
                                i++;
                                }
                        GSind++;
                        }
          
                  Contacts_file.close();
                  };           
//------------------------------------------------------------------------------
// void Umbras(int SC_num, int PL_num, const MatrixXd& orbpos, string output_path, double comp_duration, int simstep, ofstream& Umbras_file);
//------------------------------------------------------------------------------
/**
 * Compute umbra and penumbra entry and exit times of one or more spacecraft
 *
 * @param SC_num          Number of spacecraft considered
 * @param PL_num          Number of orbital plane considered
 * @param orbpos          4-D vector containing epoch (0) and ECI position (1-3) at epoch
 * @param output_path     Location of output files containing contacts information    
 * @param comp_duration   Time span of analysis [s]
 * @param simstep         Step of input ephemerides [s]
 *
 * @return .csv file containing umbra and penumbra information for all spacecraft (Umbras_file)
 * @return .csv file containing umbra and penumbra information for the spacecraft considered (written inside the function)       
 */
//------------------------------------------------------------------------------                      
void Umbras(int SC_num,
            int PL_num,
            SOLSYS Solar,
            const MatrixXd& orbpos,
            string output_path,
            double comp_duration,
            int simstep,
            ofstream& AllUmbras_file)
                  {
                  int last_element = comp_duration/simstep;
                  
                  int matrows = orbpos.rows();
                  MatrixXd lonlat = orbpos.block(0,0,matrows,3);
                  
                  int len = matrows;
          
                  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
                  /////////////////////////////////////// CONTACTS WITH GROUND STATIONS ///////////////////////////////////////
                  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
                  
                  string umbra_filename = output_path + "/S" + to_string(SC_num) + "-P" + to_string(PL_num) + "_Eclipse_times.csv";
                  //cout << umbra_filename << endl;
                  ofstream Umbras_file;
                  Umbras_file.open(umbra_filename);
                  
                  Umbras_file << "Penumbra start UTC,Umbra start UTC,Umbra end UTC,Penumbra end UTC,Umbra duration [m],Penumbra duration [m]" << endl;
                  
                  SpiceChar time_umbrain_string[35];
                  SpiceChar time_umbraout_string[35];
                  SpiceChar time_penumbrain_string[35];
                  SpiceChar time_penumbraout_string[35];
                  
                  ConstSpiceChar* format;
                  format = "YYYY-MM-DDTHR:MN:SCZ";
                  
                  bool umbra, penumbra;
                  
                  double time, sec_J2000, time_umbrain, time_umbraout, time_penumbrain, time_penumbraout, umbra_duration, penumbra_duration;
                  Vec3d posECI = Vec3d::Zero();
                  Vec4d ephem_row = Vec4d::Zero();
                  
                  int i = 0;
                  int* i_ptr = &i;
                  int ind;
                  
                  while(i < last_element)
                        {
                        ephem_row = orbpos.row(i);
                        time = ephem_row(0);
                        posECI = ephem_row.segment(1,3);
				
				time_umbrain = 0.0;
                        time_umbraout = 0.0;
                        time_penumbrain = 0.0;
                        time_penumbraout = 0.0;
                        
                        Solar.eclipse(time, posECI, umbra, penumbra);
                          
                        if(penumbra)
                              {
                              time_penumbrain = time;
                              sec_J2000 = GPS2ET(time_penumbrain);
                              
                              timout_c( sec_J2000, format, 35, time_penumbrain_string );
                              
                              ind = i + 1;
          
                              do
                                {
                                ephem_row = orbpos.row(ind);
                                time = ephem_row(0);
                                posECI = ephem_row.segment(1,3);
                                
                                Solar.eclipse(time, posECI, umbra, penumbra);
                            
                                if(umbra)
                                        {
                                        time_umbrain = time;
                                        sec_J2000 = GPS2ET(time_umbrain);
                                        timout_c( sec_J2000, format, 35, time_umbrain_string );
                                        
                                        ind++;
                                        
                                        while(umbra)
                                          {
                                          time_umbraout = time;
                                          
                                          ephem_row = orbpos.row(ind);
                                          time = ephem_row(0);
                                          posECI = ephem_row.segment(1,3);
                                        
                                          Solar.eclipse(time, posECI, umbra, penumbra);
                                        
                                          if(umbra) ind++;
                                          }
                                          
                                        //time_umbraout = time;
                                        sec_J2000 = GPS2ET(time_umbraout);
                                        timout_c( sec_J2000, format, 35, time_umbraout_string );
                                        
                                        time_penumbraout = time;
                                        }
                                
                                ind++;
                                
                                }while( penumbra && ind < len );
                              
                              time_penumbraout = time;
                              
                              sec_J2000 = GPS2ET(time_penumbraout);
                              timout_c( sec_J2000, format, 35, time_penumbraout_string );
                              
                              umbra_duration = time_umbraout - time_umbrain;
                              penumbra_duration = time_penumbraout - time_penumbrain;
					
					Umbras_file << string(time_penumbrain_string) << "," << string(time_umbrain_string) << "," << string(time_umbraout_string) << "," << string(time_penumbraout_string) << "," << umbra_duration/60.0 << "," << penumbra_duration/60.0 << "," << endl;
                              
                              AllUmbras_file << string(time_penumbrain_string) << "," << string(time_umbrain_string) << "," << string(time_umbraout_string) << "," << string(time_penumbraout_string) << "," << umbra_duration/60.0 << "," << penumbra_duration/60.0 << "," << PL_num << "," << SC_num << endl;
                                 
                              *i_ptr = ind;
                              }
                            
                          i++;
                          }
          
                  Umbras_file.close();
                  };