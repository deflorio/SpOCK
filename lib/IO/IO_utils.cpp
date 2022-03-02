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

#include "IO_utils.h"
#include <Transformations.h>
#include <Constants.h>

#include <simparam_schema-pimpl.h>

using namespace std;
using namespace SC;
using namespace constants;

//-------------------------------------------------------------------------------------
// Eigen::MatrixXd read_csvfile(const char* filename, int cols, bool header = false)
//-------------------------------------------------------------------------------------
/**
 * Load and read a csv file
 *
 * @param filename              Complete path to csv file (e.g. dir1/dir2/dir3/filenae.csv)
 * @param cols                  Number of columns to be read in csv file
 * @param header (optional)     If header == true all the lines not containing only numbers
 *                              are skipped. Default value is header = false.
 *
 * @return Eigen matrix containing the data read in the csv file
 */
//------------------------------------------------------------------------------------- 
Eigen::MatrixXd read_csvfile(const char* filename,
                             int cols,
                             bool header)
                              {
                              Eigen::MatrixXd loaded_file;
                              
                              ifstream file(filename);
                              string line;
                              
                              // Count number of lines and resize matrix loaded_file
                              int rows = 0;
                              if(file.is_open())
                                 {
                                 while( getline(file,line) ) rows++;
                                 file.close();
                                 }
                                
                              loaded_file.resize(rows,cols);
                              
                              typedef boost::tokenizer<boost::escaped_list_separator<char>> ephem_tokenizer;
                              
                              ifstream file1(filename);
                              string line1;
                              
                              if(file1.is_open())
                                 {
                                 int ind = 0;
                                 int linenum = 0;
                                 while( getline(file1,line1) )
                                        {
                                        int k = 0;
                                        bool isnum;
                                        ephem_tokenizer tok(line1);
                                        
                                        if(header)
                                            {
                                            isnum = isNumber(line1);
                                            if(!isnum)
                                                {
                                                cout << "Not considered line N. " << linenum << ": " << line1 << " of file " << filename << " because containing not number characters" << endl;
                                                rows--;
                                                loaded_file.resize(rows,cols);
                                                linenum++;
                                                continue;
                                                }
                                            else
                                                {
                                                for(ephem_tokenizer::iterator beg = tok.begin(); beg != tok.end(); ++beg)
                                                    {
                                                    //cout << stod(*beg) << "\n";
                                                    loaded_file(ind,k) = stod(*beg);
                                                    k++;
                                                    }
                                                ind ++;
                                                linenum++;
                                                }
                                            }
                                        else
                                            {
                                            for(ephem_tokenizer::iterator beg = tok.begin(); beg != tok.end(); ++beg)
                                                {
                                                //cout << stod(*beg) << "\n";
                                                loaded_file(ind,k) = stod(*beg);
                                                k++;
                                                }
                                            ind ++;
                                            }
                                        }
                                          
                                 file1.close();
                                 }
                              else
                                 {
                                 string filename_str(filename);
                                 cerr << "csv file " + filename_str + " cannot be opened" << endl;
                                 exit(EXIT_FAILURE);
                                 }
                                  
                              return loaded_file;    
                              };
//-------------------------------------------------------------------------------------
// bool isNumber(const string& str)
//-------------------------------------------------------------------------------------
/**
 * Determine if a string is a number. NOTE: letters eE has not be considered since a number
 * can be written also in scientific notation (e.g. 1e-3)
 *
 * @param String to be analyzed
 *
 * @return 1 if the string considered contains only numbers
 */
//-------------------------------------------------------------------------------------                               
bool isNumber(const string& str)                             
                  {
                  return str.find_first_of("AaBbCcDdFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXxYyZz") == string::npos;
                  }   
//-------------------------------------------------------------------------------------
// Matrix3D read_gfc(const char* filename, int maxdeg, double epoch, double& mu, double& Re)
//-------------------------------------------------------------------------------------
/**
 * Load and read a geopotential model file
 *
 * @param filename    Complete path to gfc gravity model file
 * @param maxdeg      Maximum degree and order of coefficients to be read and stored in coefficients matrices
 * @param epoch       Epoch in GPS seconds. This input is used only when gravity field models containing time variable parameters are read (e.g. EIGEN-6S)
 * @param mu          Value of Earth's gravitational coefficient read in model's file
 * @param Re          Value of Earth's equatorial radius read in model's file
 *
 * @return 3D matrix containing all gravity model coefficients read in the file
 */
//------------------------------------------------------------------------------------- 
Matrix3D read_gfc(const char* filename,
                  int maxdeg,
                  double epoch,
                  double& mu,
                  double& Re)
                        {
                        Matrix3D loaded_file;
                        loaded_file.resize(boost::extents[maxdeg+1][maxdeg+1][4]);
                        
                        ifstream file(filename);
                        string line, readtoken;
                        string filenamestring(filename);
                        
                        int mod_maxdeg = 0;
                        //int rows = 0;
                        
                        Matrix3D_index L = 0, M = 0;
                        double C, S, sigmaC, sigmaS;
                        
                        typedef boost::tokenizer< boost::char_separator<char> > grav_tokenizer;
                        boost::char_separator<char> sep("  ");
                        
                        if(file.is_open())
                            {
                            while(getline(file,line))
                                 {
                                 if(line.empty()) continue;
                                 if(isspace(line[0])) continue;
                                 
                                 grav_tokenizer tokens(line, sep);
                                 
                                 grav_tokenizer::iterator tok_iter = tokens.begin();
                                 readtoken = *tok_iter;
                                 //cout << readtoken << "\n";
                                 
                                 if(readtoken.compare("max_degree") == 0)
                                   {
                                   tok_iter++;
                                   mod_maxdeg = stoi(*tok_iter);
                                   //cout << mod_maxdeg << "\n";
                                   if(maxdeg > mod_maxdeg)
                                     {
                                     cerr << "The maximum degree and order of the selected model is" << mod_maxdeg << "x" << mod_maxdeg << endl;
                                     exit(EXIT_FAILURE);
                                     }
                                   }
                                   
                                 if(readtoken.compare("earth_gravity_constant") == 0)
                                   {
                                   tok_iter++;
                                   mu = stod(*tok_iter);
                                   //cout << mu << "\n";
                                   }
                                   
                                 if(readtoken.compare("radius") == 0)
                                   {
                                   tok_iter++;
                                   Re = stod(*tok_iter);
                                   //cout << Re << "\n";
                                   }
                                 
                                 if(readtoken.compare("gfc") == 0)
                                   {
                                   tok_iter = tokens.begin();
                                   
                                   tok_iter++;
                                   L = stoi(*tok_iter);
                                   tok_iter++;
                                   M = stoi(*tok_iter);
                                   
                                   if( (L > maxdeg) || (M > maxdeg) ) continue;
                                   if( (L == maxdeg + 1) && (M == maxdeg) ) break;
                                   
                                   tok_iter++;
                                   C = stod(*tok_iter);
                                   tok_iter++;
                                   S = stod(*tok_iter);
                                   tok_iter++;
                                   sigmaC = stod(*tok_iter);
                                   tok_iter++;
                                   sigmaS = stod(*tok_iter);
                                   
                                   //cout << L << "\n";
                                   //cout << M << "\n";
                                   
                                   loaded_file[L][M][0] = C;
                                   loaded_file[L][M][1] = S;
                                   loaded_file[L][M][2] = sigmaC;
                                   loaded_file[L][M][3] = sigmaS;
                                   
                                   //cout.precision(20);
                                   //cout << loaded_file[L][M][0] << "\n";
                                   //cout << loaded_file[L][M][1] << "\n";
                                   //cout << loaded_file[L][M][2] << "\n";
                                   //cout << loaded_file[L][M][3] << "\n";
                                   }
                                   
                                 double Ct, St, sigmaCt, sigmaSt, year, month, day, t0_GPS, delta_t;
                                 Vec4d trnd, asin1, acos1, asin2, acos2;
                                 string t0;
                                 Vector6d UTCdate;
                                   
                                 if(readtoken.compare("gfct") == 0 && (filenamestring.find("EIGEN")!=std::string::npos))
                                   {  
                                   tok_iter = tokens.begin();
                                   
                                   tok_iter++;
                                   L = stoi(*tok_iter);
                                   tok_iter++;
                                   M = stoi(*tok_iter);
                                   
                                   //cout << "\n" << L << "\n";
                                   //cout << "\n" << M << "\n";
                                   
                                   if( (L > maxdeg) || (M > maxdeg) ) continue;
                                   if( (L == maxdeg + 1) && (M == maxdeg) ) break;
                                   
                                   tok_iter++;
                                   Ct = stod(*tok_iter);
                                   tok_iter++;
                                   St = stod(*tok_iter);
                                   tok_iter++;
                                   sigmaCt = stod(*tok_iter);
                                   tok_iter++;
                                   sigmaSt = stod(*tok_iter);
                                   tok_iter++;
                                   t0 = *tok_iter;
                                   
                                   year = stod(t0.substr(0,4));
                                   month = stod(t0.substr(4,2));
                                   day = stod(t0.substr(6,2));
                                   
                                   UTCdate << year, month, day, 0.0, 0.0, 0.0;
                                   
                                   t0_GPS = UTCdate2GPSsecs(UTCdate);
                                   
                                   delta_t = (epoch - t0_GPS)/timescales::JULIAN_YEAR;//31557600.0; // Number of astronomical years (())31557600.0 = number of seconds in 365.25 days)
                                   
                                   // trnd line
                                   getline(file,line);
                                   tok_iter = tokens.begin();
                                   //cout << line << "\n";
                                   tok_iter++;
                                   //L = stoi(*tok_iter);
                                   tok_iter++;
                                   //M = stoi(*tok_iter);
                                   tok_iter++;
                                   trnd(0) = stod(*tok_iter);
                                   tok_iter++;
                                   trnd(1) = stod(*tok_iter);
                                   tok_iter++;
                                   trnd(2) = stod(*tok_iter);
                                   tok_iter++;
                                   trnd(3) = stod(*tok_iter);
                                   
                                   // acos1 line
                                   getline(file,line);
                                   tok_iter = tokens.begin();
                                   //cout << line << "\n";
                                   tok_iter++;
                                   //L = stoi(*tok_iter);
                                   tok_iter++;
                                   //M = stoi(*tok_iter);
                                   tok_iter++;
                                   acos1(0) = stod(*tok_iter);
                                   tok_iter++;
                                   acos1(1) = stod(*tok_iter);
                                   tok_iter++;
                                   acos1(2) = stod(*tok_iter);
                                   tok_iter++;
                                   acos1(3) = stod(*tok_iter);
                                   
                                   // asin1 line
                                   getline(file,line);
                                   tok_iter = tokens.begin();
                                   //cout << line << "\n";
                                   tok_iter++;
                                   //L = stoi(*tok_iter);
                                   tok_iter++;
                                   //M = stoi(*tok_iter);
                                   tok_iter++;
                                   asin1(0) = stod(*tok_iter);
                                   tok_iter++;
                                   asin1(1) = stod(*tok_iter);
                                   tok_iter++;
                                   asin1(2) = stod(*tok_iter);
                                   tok_iter++;
                                   asin1(3) = stod(*tok_iter);
                                   
                                   // acos2 line
                                   getline(file,line);
                                   tok_iter = tokens.begin();
                                   //cout << line << "\n";
                                   tok_iter++;
                                   //L = stoi(*tok_iter);
                                   tok_iter++;
                                   //M = stoi(*tok_iter);
                                   tok_iter++;
                                   acos2(0) = stod(*tok_iter);
                                   tok_iter++;
                                   acos2(1) = stod(*tok_iter);
                                   tok_iter++;
                                   acos2(2) = stod(*tok_iter);
                                   tok_iter++;
                                   acos2(3) = stod(*tok_iter);
                                   
                                   // asin2 line
                                   getline(file,line);
                                   tok_iter = tokens.begin();
                                   //cout << line << "\n";
                                   tok_iter++;
                                   //L = stoi(*tok_iter);
                                   tok_iter++;
                                   //M = stoi(*tok_iter);
                                   tok_iter++;
                                   asin2(0) = stod(*tok_iter);
                                   tok_iter++;
                                   asin2(1) = stod(*tok_iter);
                                   tok_iter++;
                                   asin2(2) = stod(*tok_iter);
                                   tok_iter++;
                                   asin2(3) = stod(*tok_iter);
                                   
                                   // Computation of coefficients
                                   C = Ct + trnd(0)*delta_t + asin1(0)*sin((mathconst::PI2/1.0)*delta_t)+acos1(0)*cos((mathconst::PI2/1.0)*delta_t) + asin2(0)*sin((mathconst::PI2/0.5)*delta_t)+acos2(0)*cos((mathconst::PI2/0.5)*delta_t);
                                   
                                   S = St + trnd(1)*delta_t + asin1(1)*sin((mathconst::PI2/1.0)*delta_t)+acos1(1)*cos((mathconst::PI2/1.0)*delta_t) + asin2(1)*sin((mathconst::PI2/0.5)*delta_t)+acos2(1)*cos((mathconst::PI2/0.5)*delta_t);
                                 
                                   sigmaC = sigmaCt + trnd(2)*delta_t + asin1(2)*sin((mathconst::PI2/1.0)*delta_t)+acos1(2)*cos((mathconst::PI2/1.0)*delta_t) + asin2(2)*sin((mathconst::PI2/0.5)*delta_t)+acos2(2)*cos((mathconst::PI2/0.5)*delta_t);
                                   
                                   sigmaS = sigmaSt + trnd(3)*delta_t + asin1(3)*sin((mathconst::PI2/1.0)*delta_t)+acos1(3)*cos((mathconst::PI2/1.0)*delta_t) + asin2(3)*sin((mathconst::PI2/0.5)*delta_t)+acos2(3)*cos((mathconst::PI2/0.5)*delta_t);
                                   
                                   loaded_file[L][M][0] = C;
                                   loaded_file[L][M][1] = S;
                                   loaded_file[L][M][2] = sigmaC;
                                   loaded_file[L][M][3] = sigmaS;
                                   }
                                   
                                  if(readtoken.compare("gfct") == 0 && (filenamestring.find("GGM")!=std::string::npos))
                                   {  
                                   tok_iter = tokens.begin();
                                   
                                   tok_iter++;
                                   L = stoi(*tok_iter);
                                   tok_iter++;
                                   M = stoi(*tok_iter);
                                   
                                   //cout << "\n" << L << "\n";
                                   //cout << "\n" << M << "\n";
                                   
                                   if( (L > maxdeg) || (M > maxdeg) ) continue;
                                   if( (L == maxdeg + 1) && (M == maxdeg) ) break;
                                   
                                   tok_iter++;
                                   Ct = stod(*tok_iter);
                                   tok_iter++;
                                   St = stod(*tok_iter);
                                   tok_iter++;
                                   sigmaCt = stod(*tok_iter);
                                   tok_iter++;
                                   sigmaSt = stod(*tok_iter);
                                   tok_iter++;
                                   t0 = "20000101";
                                   
                                   year = stod(t0.substr(0,4));
                                   month = stod(t0.substr(4,2));
                                   day = stod(t0.substr(6,2));
                                   
                                   UTCdate << year, month, day, 0.0, 0.0, 0.0;
                                   
                                   t0_GPS = UTCdate2GPSsecs(UTCdate);
                                   
                                   delta_t = (epoch - t0_GPS)/timescales::JULIAN_YEAR;//31557600.0; // Number of astronomical years (())31557600.0 = number of seconds in 365.25 days)
                                   
                                   // dot line
                                   getline(file,line);
                                   tok_iter = tokens.begin();
                                   //cout << line << "\n";
                                   tok_iter++;
                                   //L = stoi(*tok_iter);
                                   tok_iter++;
                                   //M = stoi(*tok_iter);
                                   tok_iter++;
                                   trnd(0) = stod(*tok_iter);
                                   tok_iter++;
                                   trnd(1) = stod(*tok_iter);
                                   tok_iter++;
                                   trnd(2) = stod(*tok_iter);
                                   tok_iter++;
                                   trnd(3) = stod(*tok_iter);
                                   
                                   // Computation of coefficients
                                   C = Ct + trnd(0)*delta_t;
                                   
                                   S = St + trnd(1)*delta_t;
                                 
                                   sigmaC = sigmaCt + trnd(2)*delta_t;
                                   
                                   sigmaS = sigmaSt + trnd(3)*delta_t;
                                   
                                   loaded_file[L][M][0] = C;
                                   loaded_file[L][M][1] = S;
                                   loaded_file[L][M][2] = sigmaC;
                                   loaded_file[L][M][3] = sigmaS;
                                   }
                                   //cout << readtoken << endl;
                                 }
                            file.close();
                            }
                        else
                            {
                            string filename_str(filename);
                            cerr << "Gravity model file " + filename_str + " cannot be opened" << endl;
                            exit(EXIT_FAILURE);
                            }  
                           
                        return loaded_file;    
                        };          
//-------------------------------------------------------------------------------------
// Eigen::MatrixXf read_SpaceWeather(const char* filename, double start_epoch, int sim_duration)
//-------------------------------------------------------------------------------------
/**
 * Load and read a space weather file
 * @note This function is based on the free source C++ library GeographicLib
 * @see http://geographiclib.sourceforge.net
 *
 * @param filename      Complete path to gfc gravity model file
 * @param start_epoch   Epoch from which to read the data [GPS seconds]
 * @param sim_duration   Time span of interest [s]
 *
 * @return Eigen matrix containing all space weather data in read file
 */
//------------------------------------------------------------------------------------- 
Eigen::MatrixXf read_SpaceWeather(const char* filename,
                                  double start_epoch,
                                  int sim_duration)
                                {
                                Eigen::MatrixXf loaded_file;
                                Vector6d UTCdate;
                                float year, month, day, t0_GPS, doy;
                                float end_epoch = start_epoch + sim_duration;
                                int rows = 0;
                                int cols = 32;
                                
                                int ind = 0;
                                
                                ifstream file(filename);
                                string line, readtoken, doystr, iydstr;
                                //char doychar[10];
                                string filenamestring(filename);
                                
                                loaded_file.resize(rows,cols);
                                //cout << loaded_file << endl;
                                typedef boost::tokenizer< boost::char_separator<char> > weather_tokenizer;
                                boost::char_separator<char> sep("  ");
                                
                                if(file.is_open())
                                    {
                                    while(getline(file,line))
                                        {
                                        if(line.empty()) continue;
                                        
                                        weather_tokenizer tokens(line, sep);
                                        
                                        weather_tokenizer::iterator tok_iter = tokens.begin();
                                        readtoken = *tok_iter;
                                        //cout << readtoken << "\n";
                                        if( isdigit(readtoken.c_str()[0]) == 0 ) continue;
                                        
                                        //cout << line << endl;
                                        
                                        year = stod(line.substr(0,4));
                                        month = stod(line.substr(5,2));
                                        day = stod(line.substr(8,2));
                                          
                                        UTCdate << year, month, day, 0.0, 0.0, 0.0;
                            
                                        t0_GPS = UTCdate2GPSsecs(UTCdate);
                                        
                                        if( (start_epoch - t0_GPS) > 4.0*timescales::JULIAN_DAY ) continue; // We want to store data starting from 3 days before the given start epoch
                                        if( (t0_GPS - end_epoch) > 1.0*timescales::JULIAN_DAY) break; // Stop reading the file
                                        
                                        // Build YYDOY
                                        doy = GPS2DOY(t0_GPS);
                                        //sprintf(doychar, "%03d", (int)doy);
                                        //doystr = doychar;
                                        ////iydstr = line.substr(2,2) + doystr; // YYDOY
                                        //iydstr = line.substr(0,4) + doystr; // YYYYDOY
                                        
                                        rows++;
                                        loaded_file.conservativeResize(rows,cols);
                                        loaded_file.row(ind).setZero();
                                        //cout << line << endl;
                                        
                                        loaded_file(ind,0) = t0_GPS;
                                        loaded_file(ind,1) = doy;//stod(iydstr);//
                                        if(!isspace(line.substr(11,4).c_str()[3])) loaded_file(ind,2) = stod(line.substr(11,4)); // Bartels Solar Rotation Number
                                        if(!isspace(line.substr(16,2).c_str()[1])) loaded_file(ind,3) = stod(line.substr(16,2)); // Number of Day within the Bartels 27-day cycle (01-27)
                                        if(!isspace(line.substr(19,2).c_str()[1])) loaded_file(ind,4) = stod(line.substr(19,2)); // Planetary 3-hour Range Index Kp for 0000-0300 UT
                                        if(!isspace(line.substr(22,2).c_str()[1])) loaded_file(ind,5) = stod(line.substr(22,2)); // Planetary 3-hour Range Index Kp for 0300-0600 UT
                                        if(!isspace(line.substr(25,2).c_str()[1])) loaded_file(ind,6) = stod(line.substr(25,2)); // Planetary 3-hour Range Index Kp for 0600-0900 UT
                                        if(!isspace(line.substr(28,2).c_str()[1])) loaded_file(ind,7) = stod(line.substr(28,2)); // Planetary 3-hour Range Index Kp for 0900-1200 UT
                                        if(!isspace(line.substr(31,2).c_str()[1])) loaded_file(ind,8) = stod(line.substr(31,2)); // Planetary 3-hour Range Index Kp for 1200-1500 UT
                                        if(!isspace(line.substr(34,2).c_str()[1])) loaded_file(ind,9) = stod(line.substr(34,2)); // Planetary 3-hour Range Index Kp for 1500-1800 UT
                                        if(!isspace(line.substr(37,2).c_str()[1])) loaded_file(ind,10) = stod(line.substr(37,2)); // Planetary 3-hour Range Index Kp for 1800-2100 UT
                                        if(!isspace(line.substr(40,2).c_str()[1])) loaded_file(ind,11) = stod(line.substr(40,2)); // Planetary 3-hour Range Index Kp for 2100-0000 UT
                                        if(!isspace(line.substr(43,3).c_str()[2])) loaded_file(ind,12) = stod(line.substr(43,3)); // Sum of the 8 K_p indices for the day expressed to the nearest third of a unit
                                        if(!isspace(line.substr(47,3).c_str()[2])) loaded_file(ind,13) = stod(line.substr(47,3)); // Planetary Equivalent Amplitude Ap for 0000-0300 UT
                                        if(!isspace(line.substr(51,3).c_str()[2])) loaded_file(ind,14) = stod(line.substr(51,3)); // Planetary Equivalent Amplitude Ap for 0300-0600 UT
                                        if(!isspace(line.substr(55,3).c_str()[2])) loaded_file(ind,15) = stod(line.substr(55,3)); // Planetary Equivalent Amplitude Ap for 0600-0900 UT
                                        if(!isspace(line.substr(59,3).c_str()[2])) loaded_file(ind,16) = stod(line.substr(59,3)); // Planetary Equivalent Amplitude Ap for 0900-1200 UT
                                        if(!isspace(line.substr(63,3).c_str()[2])) loaded_file(ind,17) = stod(line.substr(63,3)); // Planetary Equivalent Amplitude Ap for 1200-1500 UT
                                        if(!isspace(line.substr(67,3).c_str()[2])) loaded_file(ind,18) = stod(line.substr(67,3)); // Planetary Equivalent Amplitude Ap for 1500-1800 UT
                                        if(!isspace(line.substr(71,3).c_str()[2])) loaded_file(ind,19) = stod(line.substr(71,3)); // Planetary Equivalent Amplitude Ap for 1800-2100 UT
                                        if(!isspace(line.substr(75,3).c_str()[2])) loaded_file(ind,20) = stod(line.substr(75,3)); // Planetary Equivalent Amplitude Ap for 2100-0000 UT
                                        if(!isspace(line.substr(79,3).c_str()[2])) loaded_file(ind,21) = stod(line.substr(79,3)); // Arithmetic average of the 8 Ap indices for the day
                                        if(!isspace(line.substr(83,3).c_str()[2])) loaded_file(ind,22) = stod(line.substr(83,3)); // C_p or Planetary Daily Character Figure
                                        if(!isspace(line.substr(87,1).c_str()[0])) loaded_file(ind,23) = stod(line.substr(87,1)); // C_9. A conversion of the 0-to-2.5 range of the C_p
                                        if(!isspace(line.substr(89,3).c_str()[2])) loaded_file(ind,24) = stod(line.substr(89,3)); // International Sunspot Number
                                        if(!isspace(line.substr(93,5).c_str()[4])) loaded_file(ind,25) = stod(line.substr(93,5)); // 10.7-cm Solar Radio Flux F10.7 Adjusted to 1 AU
                                        if(!isspace(line.substr(99,1).c_str()[0])) loaded_file(ind,26) = stod(line.substr(99,1)); // Flux Qualifier
                                        if(!isspace(line.substr(101,5).c_str()[4])) loaded_file(ind,27) = stod(line.substr(101,5)); // Centered 81-day arithmetic average of F10.7 (adjusted)
                                        if(!isspace(line.substr(107,5).c_str()[4])) loaded_file(ind,28) = stod(line.substr(107,5)); // Last 81-day arithmetic average of F10.7 (adjusted)
                                        if(!isspace(line.substr(113,5).c_str()[4])) loaded_file(ind,29) = stod(line.substr(113,5)); // Observed (unadjusted) value of F10.7
                                        if(!isspace(line.substr(119,5).c_str()[4])) loaded_file(ind,30) = stod(line.substr(119,5)); // Centered 81-day arithmetic average of F10.7 (observed)
                                        if(!isspace(line.substr(125,5).c_str()[4])) loaded_file(ind,31) = stod(line.substr(125,5)); // Last 81-day arithmetic average of F10.7 (observed)
                                        
                                        ind++;
                                        }
                                    file.close();
                                    }
                                else
                                    {
                                    string filename_str(filename);
                                    cerr << "Space weather indices file " + filename_str + " cannot be opened" << endl;
                                    exit(EXIT_FAILURE);
                                    }
                                   
                                return loaded_file;    
                                };
//-------------------------------------------------------------------------------------
// void write_matfile(const Eigen::VectorXd& ephem, int cols, const char* matvarname, const char* filename)
//-------------------------------------------------------------------------------------
/**
 * Write data to a Matlab mat file
 *
 * @param ephem      Vector containing data to be written
 * @param cols       Number of columns to split data in mat file
 * @param matvarname Name of variable in Matlab environment
 * @param filename   Name of output mat file
 *
 */
//------------------------------------------------------------------------------------- 
//void write_matfile(const Eigen::VectorXd& ephem,
//                   int cols,
//                   const char* matvarname,
//                   const char* filename)
//                  {
//                  // Open .mat file
//                  mat_t *matfp;
//                  matfp = Mat_CreateVer(filename,NULL,MAT_FT_MAT5);
//                  if(NULL == matfp)
//                    {
//                     fprintf(stderr,"Error creating MAT file ephem.mat\n");
//                     //return EXIT_FAILURE;
//                    }
//                  
//                  matvar_t *matvar;
//                  
//                  int arr_len = ephem.size();
//                  //cout << arr_len << endl;
//                  double arr_ephem[arr_len] = {0};
//                  for(int k = 0; k < arr_len; k++) arr_ephem[k] = ephem(k);
//                  //cout << "Sono qui" << endl;
//                  
//                  arr_len = arr_len/cols;
//                  size_t dims[2] = {arr_len,cols};
//                  
//                  matvar = Mat_VarCreate(matvarname,MAT_C_DOUBLE,MAT_T_DOUBLE,2,dims,arr_ephem,0);
//                  if(NULL == matvar) fprintf(stderr,"Error creating matlab variable stateEME_MAIN\n");
//                  
//                  Mat_VarWrite(matfp,matvar,MAT_COMPRESSION_NONE);
//                  Mat_VarFree(matvar);
//                  Mat_Close(matfp);
//                  };
//#ifdef USE_MATLAB
//
//int Vec2matfile(const char *filename,
//                 const Eigen::VectorXd& vec_in)
//                {
//                int arr_len = vec_in.size();
//                double arr_vec_in[arr_len] = {0};
//                for(int k = 0; k < arr_len; k++) arr_vec_in[k] = vec_in(k);
//	
//                int status;
//                MATFile *pmat;
//                mxArray *mat_vec_in;
//                
//                pmat = matOpen(filename, "w");
//                if(pmat == NULL)
//	{
//	printf("Error creating file %s\n", filename);
//	printf("(Do you have write permission in this directory?)\n");
//	return(EXIT_FAILURE);
//	}
//	
//                mat_vec_in = mxCreateDoubleMatrix(1,arr_len,mxREAL);
//                
//                if(mat_vec_in == NULL)
//	{
//	printf("%s : Out of memory on line %d\n", __FILE__, __LINE__);
//	printf("Unable to create mxArray.\n");
//	return(EXIT_FAILURE);
//	}
//                memcpy((void *)(mxGetPr(mat_vec_in)), (void *)arr_vec_in, sizeof(arr_vec_in));
//                
//                status = matPutVariableAsGlobal(pmat, "GlobalDouble", mat_vec_in);
//                if(status != 0)
//	{
//	printf("Error using matPutVariableAsGlobal\n");
//	return(EXIT_FAILURE);
//	}
//                
//                mxDestroyArray(mat_vec_in);
//                
//                if(matClose(pmat) != 0)
//	{
//	printf("Error closing file %s\n",filename);
//	return(EXIT_FAILURE);
//	}
//	
//                return(0);
//                };
//                
//#endif
//-------------------------------------------------------------------------------------
// int XML_parser(...)
//-------------------------------------------------------------------------------------
/**
 * Load and read simulation parameters XML file and validate it against an XML schema file
 *
 * @param XML_simparam_file      Complete path to XML file (e.g. dir1/dir2/dir3/filenae.xml)
 * @param Orbit_ephemeris etc    Variables in which read simulation parameters are put
 *
 * @return exception against XML schema
 *
 * @note The implementation of this parser is based on CodeSynthesis XSD binding compiler and the xerces library
 * @see http://www.codesynthesis.com/products/xsd/
 * @see https://xerces.apache.org/xerces-c/
 */
//------------------------------------------------------------------------------------- 
int XML_parser(const string XML_simparam_file,
               string& Orbit_ephemeris,
               string& Attitude_ephemeris,
               string& Data_path,
               string& planetephemeris,
               string& eop,
               string& pck_data,
               string& leapsecond,
               string& magneticfield,
               string& gravityfield,
               string& atmosphere,
               string& orbfile_name,
               string& attfile_name,
               string& sensors_filename,
               string& csv_torques_name,
               string& csv_accelerations_name,
               int& SIM_STEP,
               int& SIM_DURATION,
               Vector6d& init_orbtime,
               Vector6d& init_orbstate,
               double& phi,
               double& theta,
               double& psi,
               double& om_x,
               double& om_y,
               double& om_z,
               bool& initstate_in_RTN,
               bool& realtime,
               double& realtime_wait,
               bool& ggrad_on,
               bool& mag_on,
               bool& drag_on,
               bool& srp_on,
               int& nMAX,
               bool& sunmoon_on,
               string& Drag_Model,
               string& SRP_Model,
               string& AttitudeType,
               bool& attctrl_on,
               string& AttCtrlType,
               bool& orbctrl_on,
               string& OrbCtrlType,
               double& SC_mass,
               Mat3x3d& MoI,
               Vec3d& CoG,
               double& SC_Cd,
               double& SC_Cr,
               double& SC_Area_D,
               double& SC_Area_R,
               Vec3d& Mdip,
               Face& F_Xplus,
               Face& F_Xminus,
               Face& F_Yplus,
               Face& F_Yminus,
               Face& F_Zplus,
               Face& F_Zminus,
               SYS_params& Sensor_prm_SUN,
               SYS_params& Sensor_prm_EARTH,
               SYS_params& Sensor_prm_CSS1,
               SYS_params& Sensor_prm_CSS2,
               SYS_params& Sensor_prm_CSS3,
               SYS_params& Sensor_prm_CSS4,
               SYS_params& Sensor_prm_CSS5,
               SYS_params& Sensor_prm_CSS6,
               SYS_params& Sensor_prm_MAG,
               SYS_params& Sensor_prm_MAGstowed,
               SYS_params& Sensor_prm_RS,
               SYS_params& Sensor_prm_MAGTRQ,
               SYS_params& Sensor_prm_WHEEL1,
               SYS_params& Sensor_prm_WHEEL2,
               SYS_params& Sensor_prm_WHEEL3,
               SYS_params& Solarpan1_prm,
               SYS_params& Solarpan2_prm,
               SYS_params& Solarpan3_prm,
               SYS_params& OrbitPropulsion1_prm,
               SYS_params& OrbitPropulsion2_prm,
               vector<maneuver>& all_maneuvers)
               {
               // Instantiate individual parsers.
               //
               //::simparam_pimpl simparam_p;
               //::fileheader_pimpl fileheader_p;
               //::xml_schema::string_pimpl string_p;
               //::license_pimpl license_p;
               //::xml_schema::uri_pimpl uri_p;
               //::xml_schema::date_pimpl date_p;
               //::reference_pimpl reference_p;
               //::xml_schema::any_simple_type_pimpl any_simple_type_p;
               //::SC_Faces_pimpl SC_Faces_p;
               //::Face_pimpl Face_p;
               //::Area_pimpl Area_p;
               //::AreaType_pimpl AreaType_p;
               //::Versor_pimpl Versor_p;
               //::xml_schema::double_pimpl double_p;
               //::posV_pimpl posV_p;
               //::unit_pimpl unit_p;
               //::SC_properties_pimpl SC_properties_p;
               //::InertiaMatrix_pimpl InertiaMatrix_p;
               //::MoI_pimpl MoI_p;
               //::InertiaType_pimpl InertiaType_p;
               //::CoG_pimpl CoG_p;
               //::Mass_pimpl Mass_p;
               //::MassType_pimpl MassType_p;
               //::Coefficients_pimpl Coefficients_p;
               //::PositiveNumber_pimpl PositiveNumber_p;
               //::InputFiles_pimpl InputFiles_p;
               //::OutputFiles_pimpl OutputFiles_p;
               //::SimParameters_pimpl SimParameters_p;
               //::durstep_pimpl durstep_p;
               //::initstate_pimpl initstate_p;
               //::simoptions_pimpl simoptions_p;
               //::xml_schema::boolean_pimpl boolean_p;
               //::SensorsActuators_pimpl SensorsActuators_p;
               //::RotationMatrix_3x3_pimpl RotationMatrix_3x3_p;
               
              ::simparam_pimpl simparam_p;
              ::fileheader_pimpl fileheader_p;
              ::xml_schema::string_pimpl string_p;
              ::license_pimpl license_p;
              ::xml_schema::uri_pimpl uri_p;
              ::xml_schema::date_pimpl date_p;
              ::reference_pimpl reference_p;
              ::SC_Faces_pimpl SC_Faces_p;
              ::Face_pimpl Face_p;
              ::Area_pimpl Area_p;
              ::AreaType_pimpl AreaType_p;
              ::Versor_pimpl Versor_p;
              ::xml_schema::double_pimpl double_p;
              ::posV_pimpl posV_p;
              ::SC_properties_pimpl SC_properties_p;
              ::InertiaMatrix_pimpl InertiaMatrix_p;
              ::MoI_pimpl MoI_p;
              ::InertiaType_pimpl InertiaType_p;
              ::CoG_pimpl CoG_p;
              ::Mass_pimpl Mass_p;
              ::MassType_pimpl MassType_p;
              ::Coefficients_pimpl Coefficients_p;
              ::PositiveNumber_pimpl PositiveNumber_p;
              ::Areas_pimpl Areas_p;
              ::InputFiles_pimpl InputFiles_p;
              ::OutputFiles_pimpl OutputFiles_p;
              ::SimParameters_pimpl SimParameters_p;
              ::durstep_pimpl durstep_p;
              ::xml_schema::duration_pimpl duration_p;
              ::ORB_initstate_pimpl ORB_initstate_p;
              ::xml_schema::date_time_pimpl date_time_p;
              ::velV_pimpl velV_p;
              ::ATT_initstate_pimpl ATT_initstate_p;
              ::Angle_pimpl Angle_p;
              ::AngleType_pimpl AngleType_p;
              ::simoptions_pimpl simoptions_p;
              ::xml_schema::boolean_pimpl boolean_p;
              ::Dimensioned_pimpl Dimensioned_p;
              ::nMAX_pimpl nMAX_p;
              ::SensorsActuators_pimpl SensorsActuators_p;
              ::constparam_pimpl constparam_p;
              ::auxparam_pimpl auxparam_p;
              ::opslimit_pimpl opslimit_p;
              ::accuracy_pimpl accuracy_p;
              ::RotationMatrix_3x3_pimpl RotationMatrix_3x3_p;
              ::Maneuvers_pimpl Maneuvers_p;
              ::Man_pimpl Man_p;
              ::Vector_pimpl Vector_p;
              ::name_pimpl name_p;
           
              // Connect the parsers together.
              //
              simparam_p.parsers (fileheader_p,
                        SC_Faces_p,
                        SC_properties_p,
                        InputFiles_p,
                        OutputFiles_p,
                        SimParameters_p,
                        SensorsActuators_p,
                        Maneuvers_p,
                        string_p);

              fileheader_p.parsers (string_p,
                                    string_p,
                                    string_p,
                                    license_p,
                                    string_p,
                                    date_p,
                                    string_p,
                                    string_p,
                                    string_p,
                                    string_p,
                                    reference_p);
          
              license_p.parsers (string_p,
                                 uri_p);
          
              reference_p.parsers (string_p,
                                  string_p,
                                  string_p,
                                  string_p);
          
              SC_Faces_p.parsers (Face_p);
          
              Face_p.parsers (Area_p,
                              Versor_p,
                              string_p,
                              posV_p,
                              posV_p,
                              string_p);
          
              Area_p.parsers (AreaType_p);
          
              Versor_p.parsers (double_p,
                                double_p,
                                double_p,
                                string_p);
          
              posV_p.parsers (double_p,
                              double_p,
                              double_p,
                              string_p,
                              string_p);
          
              SC_properties_p.parsers (InertiaMatrix_p,
                                       CoG_p,
                                       Coefficients_p,
                                       Areas_p,
                                       posV_p);
          
              InertiaMatrix_p.parsers (MoI_p,
                                       MoI_p,
                                       MoI_p,
                                       double_p,
                                       double_p,
                                       double_p,
                                       InertiaType_p);
          
              MoI_p.parsers (InertiaType_p);
          
              CoG_p.parsers (Mass_p,
                             posV_p);
          
              Mass_p.parsers (MassType_p);
          
              Coefficients_p.parsers (PositiveNumber_p,
                                      PositiveNumber_p);
              
              Areas_p.parsers (PositiveNumber_p,
                               PositiveNumber_p);
          
              InputFiles_p.parsers (string_p,
                                    string_p,
                                    string_p,
                                    string_p,
                                    string_p,
                                    string_p,
                                    string_p,
                                    string_p,
                                    string_p,
                                    string_p,
                                    string_p);

              OutputFiles_p.parsers (string_p,
                                     string_p,
                                     string_p,
                                     string_p,
                                     string_p,
                                     string_p);
          
              SimParameters_p.parsers (durstep_p,
                                       ORB_initstate_p,
                                       ATT_initstate_p,
                                       simoptions_p);
          
              durstep_p.parsers (duration_p,
                                 duration_p);
          
              ORB_initstate_p.parsers (date_time_p,
                             posV_p,
                             velV_p);

              velV_p.parsers (double_p,
                              double_p,
                              double_p,
                              string_p,
                              string_p);
          
              ATT_initstate_p.parsers (Angle_p,
                                       Angle_p,
                                       Angle_p,
                                       Angle_p,
                                       Angle_p,
                                       Angle_p);
          
              Angle_p.parsers (AngleType_p);
          
              simoptions_p.parsers (boolean_p,
                                    boolean_p,
                                    Dimensioned_p,
                                    boolean_p,
                                    boolean_p,
                                    boolean_p,
                                    boolean_p,
                                    nMAX_p,
                                    boolean_p,
                                    string_p,
                                    string_p,
                                    string_p,
                                    boolean_p,
                                    string_p,
                                    boolean_p,
                                    string_p);

              Dimensioned_p.parsers (string_p);
          
              SensorsActuators_p.parsers (boolean_p,
                                          constparam_p,
                                          auxparam_p,
                                          opslimit_p,
                                          accuracy_p,
                                          RotationMatrix_3x3_p,
                                          string_p);
          
              constparam_p.parsers (string_p,
                                    string_p);
          
              auxparam_p.parsers (string_p,
                                  string_p);
          
              opslimit_p.parsers (string_p,
                                  string_p);
          
              accuracy_p.parsers (string_p,
                                  string_p);
          
              RotationMatrix_3x3_p.parsers (double_p,
                                            double_p,
                                            double_p,
                                            double_p,
                                            double_p,
                                            double_p,
                                            double_p,
                                            double_p,
                                            double_p,
                                            string_p);
          
              Maneuvers_p.parsers (Man_p);
          
              Man_p.parsers (boolean_p,
                                double_p,
                                double_p,
                                Vector_p,
                                name_p);
          
              Vector_p.parsers (double_p,
                                double_p,
                                double_p,
                                string_p,
                                string_p);
               
               
              try
                {
                // Parse the XML document.
                //
                ::xml_schema::document doc_p(simparam_p, "simparam");
            
                simparam_p.pre ();
                doc_p.parse(XML_simparam_file);
                simparam_p.post_simparam ();
                }
              catch (const ::xml_schema::exception& e)
                {
                std::cerr << e << std::endl;
                return 1;
                }
                
              ///////////////////////////////////////////////////////////////
              ///////////////////////// FILES PATHS /////////////////////////
              ///////////////////////////////////////////////////////////////
              // Input files paths
              Orbit_ephemeris = InputFiles_p.Orbit_ephemeris_in;
              Attitude_ephemeris = InputFiles_p.Attitude_ephemeris_in;
              Data_path = InputFiles_p.Data_path_in;
              planetephemeris = InputFiles_p.planetephemeris_in;
              pck_data = InputFiles_p.pck_data_in;
              eop = InputFiles_p.eop_in;
              leapsecond = InputFiles_p.leapsecond_in;
              magneticfield = InputFiles_p.magn_model_in;
              gravityfield = InputFiles_p.gravityfield_in;
              atmosphere = InputFiles_p.atmosphere_in;
              // Output files paths
              orbfile_name = OutputFiles_p.orbfile_name_in;
              attfile_name = OutputFiles_p.attfile_name_in;
              sensors_filename = OutputFiles_p.sensors_filename_in;
              csv_torques_name = OutputFiles_p.csv_torques_name_in;
              csv_accelerations_name = OutputFiles_p.csv_accelerations_name_in;
              /////////////////////////////////////////////////////////////////////////
              ///////////////////////// SIMULATION PARAMETERS /////////////////////////
              /////////////////////////////////////////////////////////////////////////
              // Simulation step
              SIM_STEP = durstep_p.sim_step;
              // Simulation duration
              SIM_DURATION = durstep_p.sim_duration;
              // Initial orbit UTC date and time
              init_orbtime = ORB_initstate_p.UTCdate;
              // Initial orbit state
              Vec3d Pos_vec = ORB_initstate_p.Pos_vec;
              Vec3d Vel_vec = ORB_initstate_p.Vel_vec;
              
              init_orbstate.segment(0,3) = Pos_vec;
              init_orbstate.segment(3,3) = Vel_vec;
              // Initial attitude state
              phi = ATT_initstate_p.phi_in;
              theta = ATT_initstate_p.theta_in;
              psi = ATT_initstate_p.psi_in;
              om_x = ATT_initstate_p.om_x_in;
              om_y = ATT_initstate_p.om_y_in;
              om_z = ATT_initstate_p.om_z_in;
              // Initial state in RTN frame on/off
              initstate_in_RTN = simoptions_p.initstate_in_RTN_in;
              // Real time simulation on/off
              realtime = simoptions_p.realtime_in;
              // Real time simulation on/off
              realtime_wait = simoptions_p.realtime_wait_in;
              // Gravity gradient torque on/off
              ggrad_on = simoptions_p.ggrad_on_in;
              // Magnetic torque on/off
              mag_on = simoptions_p.mag_on_in;
              // Atmospheric drag torque on/off
              drag_on = simoptions_p.drag_on_in;
              // Solar radiation pressure torque on/off
              srp_on = simoptions_p.srp_on_in;
              // Maximum order and degree of gravitational field model used for the orbit propagation
              nMAX = simoptions_p.nMAX_in;
              // Third body perturbation on/off
              sunmoon_on = simoptions_p.sunmoon_on_in;
              // Atmospheric drag model used
              Drag_Model = simoptions_p.Drag_Model_in;
              // Solar radiation pressure model used
              SRP_Model = simoptions_p.SRP_Model_in;
              // Attitude type during orbit propagation
              AttitudeType = simoptions_p.AttitudeType_in;
              // Attitude control on/off
              attctrl_on = simoptions_p.attctrl_on_in;
              // Attitude control type
              AttCtrlType = simoptions_p.AttCtrlType_in;
              // Orbit control on/off
              orbctrl_on = simoptions_p.orbctrl_on_in;
              // Orbit control type
              OrbCtrlType = simoptions_p.OrbCtrlType_in;
              
              ////////////////////////////////////////////////////////////////////////
              ///////////////////////// SPACECRAFT PROPERTIES/////////////////////////
              ////////////////////////////////////////////////////////////////////////
              // Spacecraft mass
              SC_mass = CoG_p.SC_mass_in;
              // Spacecraft center of mmass position vector in body-fixed coordinates
              CoG(0) = CoG_p.CoG_pos_vec(0);
              CoG(1) = CoG_p.CoG_pos_vec(1);
              CoG(2) = CoG_p.CoG_pos_vec(2);
              // Moments of inertia matrix. Moment of inertia taken at the center of mass and aligned with the body-fixed frame [kg*m^2]
              MoI(0,0) = InertiaMatrix_p.Ixx_in;
              MoI(1,1) = InertiaMatrix_p.Iyy_in;
              MoI(2,2) = InertiaMatrix_p.Izz_in;
              MoI(0,1) = -InertiaMatrix_p.Ixy_in;
              MoI(0,2) = -InertiaMatrix_p.Ixz_in;
              MoI(1,2) = -InertiaMatrix_p.Iyz_in;
              
              MoI(1,0) = MoI(0,1);
              MoI(2,0) = MoI(0,2);
              MoI(2,1) = MoI(1,2);
              // Drag coefficient
              SC_Cd = Coefficients_p.Cd_in;
              // SRP coefficient
              SC_Cr = Coefficients_p.Cr_in;
              // Drag area to be used with atmospheric drag simple model
              SC_Area_D = Areas_p.Area_D_in;
              // Radiation area to be used with solar radiation pressure simple model
              SC_Area_R = Areas_p.Area_R_in;
              // Spacecraft magnetic dipole moment vector in body-fixed coordinates
              Mdip(0) = SC_properties_p.Mdip_in(0);
              Mdip(1) = SC_properties_p.Mdip_in(1);
              Mdip(2) = SC_properties_p.Mdip_in(2);
              // Spacecraft faces
              F_Xplus.Area = SC_Faces_p.SC_Face_in[0].Area;
              F_Xplus.Material = SC_Faces_p.SC_Face_in[0].Material;
              F_Xplus.cP << SC_Faces_p.SC_Face_in[0].cP(0), SC_Faces_p.SC_Face_in[0].cP(1), SC_Faces_p.SC_Face_in[0].cP(2);
              F_Xplus.cA << SC_Faces_p.SC_Face_in[0].cA(0), SC_Faces_p.SC_Face_in[0].cA(1), SC_Faces_p.SC_Face_in[0].cA(2);
              
              F_Xminus.Area = SC_Faces_p.SC_Face_in[1].Area;
              F_Xminus.Material = SC_Faces_p.SC_Face_in[1].Material;
              F_Xminus.cP << SC_Faces_p.SC_Face_in[1].cP(0), SC_Faces_p.SC_Face_in[1].cP(1), SC_Faces_p.SC_Face_in[1].cP(2);
              F_Xminus.cA << SC_Faces_p.SC_Face_in[1].cA(0), SC_Faces_p.SC_Face_in[1].cA(1), SC_Faces_p.SC_Face_in[1].cA(2);
              
              F_Yplus.Area = SC_Faces_p.SC_Face_in[2].Area;
              F_Yplus.Material = SC_Faces_p.SC_Face_in[2].Material;
              F_Yplus.cP << SC_Faces_p.SC_Face_in[2].cP(0), SC_Faces_p.SC_Face_in[2].cP(1), SC_Faces_p.SC_Face_in[2].cP(2);
              F_Yplus.cA << SC_Faces_p.SC_Face_in[2].cA(0), SC_Faces_p.SC_Face_in[2].cA(1), SC_Faces_p.SC_Face_in[2].cA(2);
              
              F_Yminus.Area = SC_Faces_p.SC_Face_in[3].Area;
              F_Yminus.Material = SC_Faces_p.SC_Face_in[3].Material;
              F_Yminus.cP << SC_Faces_p.SC_Face_in[3].cP(0), SC_Faces_p.SC_Face_in[3].cP(1), SC_Faces_p.SC_Face_in[3].cP(2);
              F_Yminus.cA << SC_Faces_p.SC_Face_in[3].cA(0), SC_Faces_p.SC_Face_in[3].cA(1), SC_Faces_p.SC_Face_in[3].cA(2);
              
              F_Zplus.Area = SC_Faces_p.SC_Face_in[4].Area;
              F_Zplus.Material = SC_Faces_p.SC_Face_in[4].Material;
              F_Zplus.cP << SC_Faces_p.SC_Face_in[4].cP(0), SC_Faces_p.SC_Face_in[4].cP(1), SC_Faces_p.SC_Face_in[4].cP(2);
              F_Zplus.cA << SC_Faces_p.SC_Face_in[4].cA(0), SC_Faces_p.SC_Face_in[4].cA(1), SC_Faces_p.SC_Face_in[4].cA(2);
              
              F_Zminus.Area = SC_Faces_p.SC_Face_in[5].Area;
              F_Zminus.Material = SC_Faces_p.SC_Face_in[5].Material;
              F_Zminus.cP << SC_Faces_p.SC_Face_in[5].cP(0), SC_Faces_p.SC_Face_in[5].cP(1), SC_Faces_p.SC_Face_in[5].cP(2);
              F_Zminus.cA << SC_Faces_p.SC_Face_in[5].cA(0), SC_Faces_p.SC_Face_in[5].cA(1), SC_Faces_p.SC_Face_in[5].cA(2);
              
              //////////////////////////////////////////////////////////////
              ///////////////// ADCS sensors and actuators /////////////////
              //////////////////////////////////////////////////////////////
              
              // Sun camera
              Sensor_prm_SUN = simparam_p.SensAct_prms_in[0];
              // Earth camera
              Sensor_prm_EARTH = simparam_p.SensAct_prms_in[1];
              // Coarse Sun sensor (CSS)
              Sensor_prm_CSS1 = simparam_p.SensAct_prms_in[2];
              // Coarse Sun sensor (CSS)
              Sensor_prm_CSS2 = simparam_p.SensAct_prms_in[3];
              // Coarse Sun sensor (CSS)
              Sensor_prm_CSS3 = simparam_p.SensAct_prms_in[4];
              // Coarse Sun sensor (CSS)
              Sensor_prm_CSS4 = simparam_p.SensAct_prms_in[5];
              // Coarse Sun sensor (CSS)
              Sensor_prm_CSS5 = simparam_p.SensAct_prms_in[6];
              // Coarse Sun sensor (CSS)
              Sensor_prm_CSS6 = simparam_p.SensAct_prms_in[7];
              // Magnetometer (deployed)
              Sensor_prm_MAG = simparam_p.SensAct_prms_in[8];
              // Magnetometer (stowed)
              Sensor_prm_MAGstowed = simparam_p.SensAct_prms_in[9];
              // Rate sensor
              Sensor_prm_RS = simparam_p.SensAct_prms_in[10];
              // Magnetic torquer
              Sensor_prm_MAGTRQ = simparam_p.SensAct_prms_in[11];
              // Reaction/momentum wheels
              Sensor_prm_WHEEL1 = simparam_p.SensAct_prms_in[12];
              Sensor_prm_WHEEL2 = simparam_p.SensAct_prms_in[13];
              Sensor_prm_WHEEL3 = simparam_p.SensAct_prms_in[14];
              
              // Solar panels
              Solarpan1_prm = simparam_p.SensAct_prms_in[15];
              //SP_epsilon = SensorsActuators_p.auxparam_in(0);
              
              Solarpan2_prm = simparam_p.SensAct_prms_in[16];
              //SP_epsilon = SensorsActuators_p.auxparam_in(0);
              
              Solarpan3_prm = simparam_p.SensAct_prms_in[17];
              //SP_epsilon = SensorsActuators_p.auxparam_in(0);
              
              // Propulsion systems
              OrbitPropulsion1_prm = simparam_p.SensAct_prms_in[18];
              OrbitPropulsion2_prm = simparam_p.SensAct_prms_in[19];
              
              //////////////////////////////////////////////////////////////
              //////////////// Commanded attitude maneuvers ////////////////
              //////////////////////////////////////////////////////////////
              
              all_maneuvers = Maneuvers_p.all_maneuvers;
              
              return(0);
              };
              
              
              
              
              
              
              
              
              
              
              
//-------------------------------------------------------------------------------------
// int SGP4_XML_parser(...)
//-------------------------------------------------------------------------------------
/**
 * Load and read simulation parameters XML file and validate it against an XML schema file.
 * This is a reduced version of function XML_parser for the SGP4 propagator which requires
 * a small subset of the inputs required by the precise orbit propagation.
 *
 * @param XML_simparam_file      Complete path to XML file (e.g. dir1/dir2/dir3/filenae.xml)
 *
 * @return exception against XML schema
 *
 * @note The implementation of this parser is based on CodeSynthesis XSD binding compiler and the xerces library
 * @see http://www.codesynthesis.com/products/xsd/
 * @see https://xerces.apache.org/xerces-c/
 */
//------------------------------------------------------------------------------------- 
int SGP4_XML_parser(const string XML_simparam_file,
                string& TLE_file,
                string& Data_path,
                string& eop,
                string& pck_data,
                string& leapsecond,
                string& orbfile_name,
                int& SIM_STEP,
                int& SIM_DURATION)
                {
                ::simparam_pimpl simparam_p;
                ::fileheader_pimpl fileheader_p;
                ::xml_schema::string_pimpl string_p;
                ::license_pimpl license_p;
                ::xml_schema::uri_pimpl uri_p;
                ::xml_schema::date_pimpl date_p;
                ::reference_pimpl reference_p;
                ::SC_Faces_pimpl SC_Faces_p;
                ::SC_properties_pimpl SC_properties_p;
                ::InputFiles_pimpl InputFiles_p;
                ::OutputFiles_pimpl OutputFiles_p;
                ::SimParameters_pimpl SimParameters_p;
                ::durstep_pimpl durstep_p;
                ::xml_schema::duration_pimpl duration_p;
                ::ORB_initstate_pimpl ORB_initstate_p;
                ::ATT_initstate_pimpl ATT_initstate_p;
                ::simoptions_pimpl simoptions_p;
                ::SensorsActuators_pimpl SensorsActuators_p;
                ::Maneuvers_pimpl Maneuvers_p;
             
                // Connect the parsers together.
                //
                simparam_p.parsers (fileheader_p,
                          SC_Faces_p,
                          SC_properties_p,
                          InputFiles_p,
                          OutputFiles_p,
                          SimParameters_p,
                          SensorsActuators_p,
                          Maneuvers_p,
                          string_p);
  
                fileheader_p.parsers (string_p,
                                      string_p,
                                      string_p,
                                      license_p,
                                      string_p,
                                      date_p,
                                      string_p,
                                      string_p,
                                      string_p,
                                      string_p,
                                      reference_p);
            
                license_p.parsers (string_p,
                                   uri_p);
            
                reference_p.parsers (string_p,
                                    string_p,
                                    string_p,
                                    string_p);
            
                InputFiles_p.parsers (string_p,
                                      string_p,
                                      string_p,
                                      string_p,
                                      string_p,
                                      string_p,
                                      string_p,
                                      string_p,
                                      string_p,
                                      string_p,
                                      string_p);
  
                OutputFiles_p.parsers (string_p,
                                       string_p,
                                       string_p,
                                       string_p,
                                       string_p,
                                       string_p);
            
                SimParameters_p.parsers (durstep_p,
                                         ORB_initstate_p,
                                         ATT_initstate_p,
                                         simoptions_p);
            
                durstep_p.parsers (duration_p,
                                   duration_p);
                 
                try
                  {
                  // Parse the XML document.
                  //
                  ::xml_schema::document doc_p(simparam_p, "simparam");
              
                  simparam_p.pre ();
                  doc_p.parse(XML_simparam_file);
                  simparam_p.post_simparam ();
                  }
                catch (const ::xml_schema::exception& e)
                  {
                  std::cerr << e << std::endl;
                  return 1;
                  }
                  
                ///////////////////////////////////////////////////////////////
                ///////////////////////// FILES PATHS /////////////////////////
                ///////////////////////////////////////////////////////////////
                // Input files paths
                TLE_file = InputFiles_p.Orbit_ephemeris_in;
                Data_path = InputFiles_p.Data_path_in;
                pck_data = InputFiles_p.pck_data_in;
                eop = InputFiles_p.eop_in;
                leapsecond = InputFiles_p.leapsecond_in;
                // Output files paths
                orbfile_name = OutputFiles_p.orbfile_name_in;
                /////////////////////////////////////////////////////////////////////////
                ///////////////////////// SIMULATION PARAMETERS /////////////////////////
                /////////////////////////////////////////////////////////////////////////
                // Simulation step
                SIM_STEP = durstep_p.sim_step;
                // Simulation duration
                SIM_DURATION = durstep_p.sim_duration;
                
                return(0);
              };              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
//-------------------------------------------------------------------------------------
// void ReadXMLtoTXT(...)
//-------------------------------------------------------------------------------------
/**
 * Write simulation parameters read by parser XML_parser in a text file for verification purposes
 *
 * @param txt_file               Complete path to txt file (e.g. dir1/dir2/dir3/filenae.txt)
 * @param Orbit_ephemeris etc    Simulation parameters
 *
 */
//------------------------------------------------------------------------------------- 
void ReadXMLtoTXT(const string txt_file,
                  string Orbit_ephemeris,
                  string Attitude_ephemeris,
                  string Data_path,
                  string planetephemeris,
                  string eop,
                  string pck_data,
                  string leapsecond,
                  string magneticfield,
                  string gravityfield,
                  string atmosphere,
                  string orbfile_name,
                  string attfile_name,
                  string sensors_filename,
                  string csv_torques_name,
                  string csv_accelerations_name,
                  int SIM_STEP,
                  int SIM_DURATION,
                  Vector6d init_orbtime,
                  Vector6d init_orbstate,
                  double phi,
                  double theta,
                  double psi,
                  double om_x,
                  double om_y,
                  double om_z,
                  bool initstate_in_RTN,
                  bool realtime,
                  double realtime_wait,
                  bool ggrad_on,
                  bool mag_on,
                  bool drag_on,
                  bool srp_on,
                  int nMAX,
                  bool sunmoon_on,
                  string Drag_Model,
                  string SRP_Model,
                  string AttitudeType,
                  bool attctrl_on,
                  string AttCtrlType,
                  bool orbctrl_on,
                  string OrbCtrlType,
                  double SC_mass,
                  Mat3x3d MoI,
                  Vec3d CoG,
                  double SC_Cd,
                  double SC_Cr,
                  double SC_Area_D,
                  double SC_Area_R,
                  Vec3d Mdip,
                  Face F_Xplus,
                  Face F_Xminus,
                  Face F_Yplus,
                  Face F_Yminus,
                  Face F_Zplus,
                  Face F_Zminus,
                  SYS_params Sensor_prm_SUN,
                  SYS_params Sensor_prm_EARTH,
                  SYS_params Sensor_prm_CSS1,
                  SYS_params Sensor_prm_CSS2,
                  SYS_params Sensor_prm_CSS3,
                  SYS_params Sensor_prm_CSS4,
                  SYS_params Sensor_prm_CSS5,
                  SYS_params Sensor_prm_CSS6,
                  SYS_params Sensor_prm_MAG,
                  SYS_params Sensor_prm_MAGstowed,
                  SYS_params Sensor_prm_RS,
                  SYS_params Sensor_prm_MAGTRQ,
                  SYS_params Sensor_prm_WHEEL1,
                  SYS_params Sensor_prm_WHEEL2,
                  SYS_params Sensor_prm_WHEEL3,
                  SYS_params Solarpan1_prm,
                  SYS_params Solarpan2_prm,
                  SYS_params Solarpan3_prm,
                  SYS_params OrbitPropulsion1_prm,
                  SYS_params OrbitPropulsion2_prm,
                  vector<maneuver> all_maneuvers)
                  {
                  int maneuvers_num = all_maneuvers.size();
                  
                  ofstream txtfile;
                  txtfile.open(txt_file);
                  
                  txtfile << "######################" << endl;
                  txtfile << "INPUT FILES PATHS" << endl;
                  txtfile << "######################\n" << endl;
                  txtfile << "Orbit ephemeris: " << Orbit_ephemeris << endl;
                  txtfile << "Attitude ephemeris: " << Attitude_ephemeris << endl;
                  txtfile << "Data path: " << Data_path << endl;
                  txtfile << "Planet ephemeris: " << planetephemeris << endl;
                  txtfile << "EOP: " << eop << endl;
                  txtfile << "PCK: " << pck_data << endl;
                  txtfile << "Leap second: " << leapsecond << endl;
                  txtfile << "Magnetic field: " << magneticfield << endl;
                  txtfile << "Gravity field: " << gravityfield << endl;
                  txtfile << "Atmospheric model: " << atmosphere << endl;
                  txtfile << "\n######################" << endl;
                  txtfile << "OUTPUT FILES PATHS" << endl;
                  txtfile << "######################\n" << endl;
                  txtfile << "Orbit state: " << orbfile_name << endl;
                  txtfile << "Attitude state: " << attfile_name << endl;
                  txtfile << "Sensors readings: " << sensors_filename << endl;
                  txtfile << "Torques: " << csv_torques_name << endl;
                  txtfile << "Accelerations: " << csv_accelerations_name << endl;
                  txtfile << "\n######################" << endl;
                  txtfile << "SIMULATION OPTIONS" << endl;
                  txtfile << "######################\n" << endl;
                  txtfile << "SIM_STEP: " << SIM_STEP << endl;
                  txtfile << "SIM_DURATION: " << SIM_DURATION << "\n" << endl;
                  txtfile << "Orbit initial epoch: " << init_orbtime(0) << "-" << init_orbtime(1) << "-" << init_orbtime(2) << "/" << init_orbtime(3) << ":" << init_orbtime(4) << ":" << init_orbtime(5) << "\n" << endl;
                  txtfile << "Orbit initial state: " << "X = " << init_orbstate(0) << " Y = " << init_orbstate(1) << " Z = " << init_orbstate(2) << " Vx = " << init_orbstate(3) << " Vy = " << init_orbstate(4) << " Vz = " << init_orbstate(5) << "\n" << endl;
                  txtfile << "phi: " << phi << endl;
                  txtfile << "theta: " << theta << endl;
                  txtfile << "psi: " << psi << endl;
                  txtfile << "om_x: " << om_x << endl;
                  txtfile << "om_y: " << om_y << endl;
                  txtfile << "om_z: " << om_z << "\n" << endl;
                  txtfile << "initstate_in_RTN: " << initstate_in_RTN << endl;
                  txtfile << "realtime: " << realtime << endl;
                  txtfile << "realtime_wait: " << realtime_wait << endl;
                  txtfile << "ggrad_on: " << ggrad_on << endl;
                  txtfile << "mag_on: " << mag_on << endl;
                  txtfile << "drag_on: " << drag_on << endl;
                  txtfile << "srp_on: " << srp_on << endl;
                  txtfile << "nMAX: " << nMAX << endl;
                  txtfile << "sunmoon_on: " << sunmoon_on << endl;
                  txtfile << "Drag_Model: " << Drag_Model << endl;
                  txtfile << "SRP_Model: " << SRP_Model << endl;
                  txtfile << "AttitudeType: " << AttitudeType << endl;
                  txtfile << "attctrl_on: " << attctrl_on << endl;
                  txtfile << "AttCtrlType: " << AttCtrlType << endl;
                  txtfile << "orbctrl_on: " << orbctrl_on << endl;
                  txtfile << "OrbCtrlType: " << OrbCtrlType << endl;
                  txtfile << "\n######################" << endl;
                  txtfile << "SPACECRAFT PROPERTIES" << endl;
                  txtfile << "######################\n" << endl;
                  txtfile << "Spacecraft mass: " << SC_mass << endl;
                  txtfile << "Spacecraft center of mass position vector: " << CoG(0) << " , " << CoG(1) << " , " << CoG(2) << " , " << endl;
                  txtfile << "Spacecraft inertia matrix: " << MoI(0,0) << " , " << MoI(0,1) << " , " << MoI(0,2) << " , " << endl;
                  txtfile << "                           " << MoI(1,0) << " , " << MoI(1,1) << " , " << MoI(1,2) << " , " << endl;
                  txtfile << "                           " << MoI(2,0) << " , " << MoI(2,1) << " , " << MoI(2,2) << " , " << "\n" << endl;
                  txtfile << "Drag coefficient: " << SC_Cd << endl;
                  txtfile << "SRP coefficient: " << SC_Cr << endl;
                  txtfile << "Drag area: " << SC_Area_D << endl;
                  txtfile << "SRP area: " << SC_Area_R << "\n" << endl;
                  
                  txtfile << "Spacecraft F_Xplus: " << endl;
                  txtfile << "Area: " << F_Xplus.Area << endl;
                  txtfile << "Material: " << F_Xplus.Material << endl;
                  txtfile << "cP position: " << F_Xplus.cP(0) << " , " << F_Xplus.cP(1) << " , " << F_Xplus.cP(2) << " , " << endl;
                  txtfile << "cA position: " << F_Xplus.cA(0) << " , " << F_Xplus.cA(1) << " , " << F_Xplus.cA(2) << " , " << "\n" << endl;
                  
                  txtfile << "Spacecraft F_Xminus: " << endl;
                  txtfile << "Area: " << F_Xminus.Area << endl;
                  txtfile << "Material: " << F_Xminus.Material << endl;
                  txtfile << "cP position: " << F_Xminus.cP(0) << " , " << F_Xminus.cP(1) << " , " << F_Xminus.cP(2) << " , " << endl;
                  txtfile << "cA position: " << F_Xminus.cA(0) << " , " << F_Xminus.cA(1) << " , " << F_Xminus.cA(2) << " , " << "\n" << endl;
                  
                  txtfile << "Spacecraft F_Yplus: " << endl;
                  txtfile << "Area: " << F_Yplus.Area << endl;
                  txtfile << "Material: " << F_Yplus.Material << endl;
                  txtfile << "cP position: " << F_Yplus.cP(0) << " , " << F_Yplus.cP(1) << " , " << F_Yplus.cP(2) << " , " << endl;
                  txtfile << "cA position: " << F_Yplus.cA(0) << " , " << F_Yplus.cA(1) << " , " << F_Yplus.cA(2) << " , " << "\n" << endl;
                  
                  txtfile << "Spacecraft F_Yminus: " << endl;
                  txtfile << "Area: " << F_Yminus.Area << endl;
                  txtfile << "Material: " << F_Yminus.Material << endl;
                  txtfile << "cP position: " << F_Yminus.cP(0) << " , " << F_Yminus.cP(1) << " , " << F_Yminus.cP(2) << " , " << endl;
                  txtfile << "cA position: " << F_Yminus.cA(0) << " , " << F_Yminus.cA(1) << " , " << F_Yminus.cA(2) << " , " << "\n" << endl;
                  
                  txtfile << "Spacecraft F_Zplus: " << endl;
                  txtfile << "Area: " << F_Zplus.Area << endl;
                  txtfile << "Material: " << F_Zplus.Material << endl;
                  txtfile << "cP position: " << F_Zplus.cP(0) << " , " << F_Zplus.cP(1) << " , " << F_Zplus.cP(2) << " , " << endl;
                  txtfile << "cA position: " << F_Zplus.cA(0) << " , " << F_Zplus.cA(1) << " , " << F_Zplus.cA(2) << " , " << "\n" << endl;
                  
                  txtfile << "Spacecraft F_Zminus: " << endl;
                  txtfile << "Area: " << F_Zminus.Area << endl;
                  txtfile << "Material: " << F_Zminus.Material << endl;
                  txtfile << "cP position: " << F_Zminus.cP(0) << " , " << F_Zminus.cP(1) << " , " << F_Zminus.cP(2) << " , " << endl;
                  txtfile << "cA position: " << F_Zminus.cA(0) << " , " << F_Zminus.cA(1) << " , " << F_Zminus.cA(2) << " , " << endl;
                  
                  txtfile << "\n######################" << endl;
                  txtfile << "ATTITUDE SENSORS AND ACTUATORS" << endl;
                  txtfile << "######################\n" << endl;
                  
                  txtfile << "Sun sensor: " << endl;
                  txtfile << "On/off: " << Sensor_prm_SUN.on_off << endl;
                  txtfile << "Name: " << Sensor_prm_SUN.Name << endl;
                  txtfile << "Mounting matrix: " << Sensor_prm_SUN.SC2SYS(0,0) << " , " << Sensor_prm_SUN.SC2SYS(0,1) << " , " << Sensor_prm_SUN.SC2SYS(0,2) << " , " << endl;
                  txtfile << "                 " << Sensor_prm_SUN.SC2SYS(1,0) << " , " << Sensor_prm_SUN.SC2SYS(1,1) << " , " << Sensor_prm_SUN.SC2SYS(1,2) << " , " << endl;
                  txtfile << "                 " << Sensor_prm_SUN.SC2SYS(2,0) << " , " << Sensor_prm_SUN.SC2SYS(2,1) << " , " << Sensor_prm_SUN.SC2SYS(2,2) << " , " << "\n" << endl;
                  txtfile << "OPS_limits: " << Sensor_prm_SUN.OPS_limits(0) << endl;
                  txtfile << "accuracy: " << Sensor_prm_SUN.Accuracy(0) << " , " << Sensor_prm_SUN.Accuracy(1) << "\n" << endl;
                  
                  txtfile << "Earth sensor: " << endl;
                  txtfile << "On/off: " << Sensor_prm_EARTH.on_off << endl;
                  txtfile << "Name: " << Sensor_prm_EARTH.Name << endl;
                  txtfile << "Mounting matrix: " << Sensor_prm_EARTH.SC2SYS(0,0) << " , " << Sensor_prm_EARTH.SC2SYS(0,1) << " , " << Sensor_prm_EARTH.SC2SYS(0,2) << " , " << endl;
                  txtfile << "                 " << Sensor_prm_EARTH.SC2SYS(1,0) << " , " << Sensor_prm_EARTH.SC2SYS(1,1) << " , " << Sensor_prm_EARTH.SC2SYS(1,2) << " , " << endl;
                  txtfile << "                 " << Sensor_prm_EARTH.SC2SYS(2,0) << " , " << Sensor_prm_EARTH.SC2SYS(2,1) << " , " << Sensor_prm_EARTH.SC2SYS(2,2) << " , " << "\n" << endl;
                  txtfile << "OPS_limits: " << Sensor_prm_EARTH.OPS_limits(0) << endl;
                  txtfile << "accuracy: " << Sensor_prm_EARTH.Accuracy(0) << " , " << Sensor_prm_EARTH.Accuracy(1) << "\n" << endl;
                  
                  txtfile << "Coarse Sun sensor (CSS1): " << endl;
                  txtfile << "On/off: " << Sensor_prm_CSS1.on_off << endl;
                  txtfile << "Name: " << Sensor_prm_CSS1.Name << endl;
                  txtfile << "constparam: " << Sensor_prm_CSS1.ConstPrm(0) << " , " << Sensor_prm_CSS1.ConstPrm(1) << endl;
                  txtfile << "Mounting matrix: " << Sensor_prm_CSS1.SC2SYS(0,0) << " , " << Sensor_prm_CSS1.SC2SYS(0,1) << " , " << Sensor_prm_CSS1.SC2SYS(0,2) << " , " << endl;
                  txtfile << "                 " << Sensor_prm_CSS1.SC2SYS(1,0) << " , " << Sensor_prm_CSS1.SC2SYS(1,1) << " , " << Sensor_prm_CSS1.SC2SYS(1,2) << " , " << endl;
                  txtfile << "                 " << Sensor_prm_CSS1.SC2SYS(2,0) << " , " << Sensor_prm_CSS1.SC2SYS(2,1) << " , " << Sensor_prm_CSS1.SC2SYS(2,2) << " , " << "\n" << endl;
                  txtfile << "accuracy: " << Sensor_prm_CSS1.Accuracy(0) << " , " << Sensor_prm_CSS1.Accuracy(1) << "\n" << endl;
                  
                  txtfile << "Coarse Sun sensor (CSS2): " << endl;
                  txtfile << "On/off: " << Sensor_prm_CSS2.on_off << endl;
                  txtfile << "Name: " << Sensor_prm_CSS2.Name << endl;
                  txtfile << "constparam: " << Sensor_prm_CSS2.ConstPrm(0) << " , " << Sensor_prm_CSS2.ConstPrm(1) << endl;
                  txtfile << "Mounting matrix: " << Sensor_prm_CSS2.SC2SYS(0,0) << " , " << Sensor_prm_CSS2.SC2SYS(0,1) << " , " << Sensor_prm_CSS2.SC2SYS(0,2) << " , " << endl;
                  txtfile << "                 " << Sensor_prm_CSS2.SC2SYS(1,0) << " , " << Sensor_prm_CSS2.SC2SYS(1,1) << " , " << Sensor_prm_CSS2.SC2SYS(1,2) << " , " << endl;
                  txtfile << "                 " << Sensor_prm_CSS2.SC2SYS(2,0) << " , " << Sensor_prm_CSS2.SC2SYS(2,1) << " , " << Sensor_prm_CSS2.SC2SYS(2,2) << " , " << "\n" << endl;
                  txtfile << "accuracy: " << Sensor_prm_CSS2.Accuracy(0) << " , " << Sensor_prm_CSS2.Accuracy(1) << "\n" << endl;
                  
                  txtfile << "Coarse Sun sensor (CSS3): " << endl;
                  txtfile << "On/off: " << Sensor_prm_CSS3.on_off << endl;
                  txtfile << "Name: " << Sensor_prm_CSS3.Name << endl;
                  txtfile << "constparam: " << Sensor_prm_CSS3.ConstPrm(0) << " , " << Sensor_prm_CSS3.ConstPrm(1) << endl;
                  txtfile << "Mounting matrix: " << Sensor_prm_CSS3.SC2SYS(0,0) << " , " << Sensor_prm_CSS3.SC2SYS(0,1) << " , " << Sensor_prm_CSS3.SC2SYS(0,2) << " , " << endl;
                  txtfile << "                 " << Sensor_prm_CSS3.SC2SYS(1,0) << " , " << Sensor_prm_CSS3.SC2SYS(1,1) << " , " << Sensor_prm_CSS3.SC2SYS(1,2) << " , " << endl;
                  txtfile << "                 " << Sensor_prm_CSS3.SC2SYS(2,0) << " , " << Sensor_prm_CSS3.SC2SYS(2,1) << " , " << Sensor_prm_CSS3.SC2SYS(2,2) << " , " << "\n" << endl;
                  txtfile << "accuracy: " << Sensor_prm_CSS3.Accuracy(0) << " , " << Sensor_prm_CSS3.Accuracy(1) << "\n" << endl;
                  
                  txtfile << "Coarse Sun sensor (CSS4): " << endl;
                  txtfile << "On/off: " << Sensor_prm_CSS4.on_off << endl;
                  txtfile << "Name: " << Sensor_prm_CSS4.Name << endl;
                  txtfile << "constparam: " << Sensor_prm_CSS4.ConstPrm(0) << " , " << Sensor_prm_CSS4.ConstPrm(1) << endl;
                  txtfile << "Mounting matrix: " << Sensor_prm_CSS4.SC2SYS(0,0) << " , " << Sensor_prm_CSS4.SC2SYS(0,1) << " , " << Sensor_prm_CSS4.SC2SYS(0,2) << " , " << endl;
                  txtfile << "                 " << Sensor_prm_CSS4.SC2SYS(1,0) << " , " << Sensor_prm_CSS4.SC2SYS(1,1) << " , " << Sensor_prm_CSS4.SC2SYS(1,2) << " , " << endl;
                  txtfile << "                 " << Sensor_prm_CSS4.SC2SYS(2,0) << " , " << Sensor_prm_CSS4.SC2SYS(2,1) << " , " << Sensor_prm_CSS4.SC2SYS(2,2) << " , " << "\n" << endl;
                  txtfile << "accuracy: " << Sensor_prm_CSS4.Accuracy(0) << " , " << Sensor_prm_CSS4.Accuracy(1) << "\n" << endl;
                  
                  txtfile << "Coarse Sun sensor (CSS5): " << endl;
                  txtfile << "On/off: " << Sensor_prm_CSS5.on_off << endl;
                  txtfile << "Name: " << Sensor_prm_CSS5.Name << endl;
                  txtfile << "constparam: " << Sensor_prm_CSS5.ConstPrm(0) << " , " << Sensor_prm_CSS5.ConstPrm(1) << endl;
                  txtfile << "Mounting matrix: " << Sensor_prm_CSS5.SC2SYS(0,0) << " , " << Sensor_prm_CSS5.SC2SYS(0,1) << " , " << Sensor_prm_CSS5.SC2SYS(0,2) << " , " << endl;
                  txtfile << "                 " << Sensor_prm_CSS5.SC2SYS(1,0) << " , " << Sensor_prm_CSS5.SC2SYS(1,1) << " , " << Sensor_prm_CSS5.SC2SYS(1,2) << " , " << endl;
                  txtfile << "                 " << Sensor_prm_CSS5.SC2SYS(2,0) << " , " << Sensor_prm_CSS5.SC2SYS(2,1) << " , " << Sensor_prm_CSS5.SC2SYS(2,2) << " , " << "\n" << endl;
                  txtfile << "accuracy: " << Sensor_prm_CSS5.Accuracy(0) << " , " << Sensor_prm_CSS5.Accuracy(1) << "\n" << endl;
                  
                  txtfile << "Coarse Sun sensor (CSS6): " << endl;
                  txtfile << "On/off: " << Sensor_prm_CSS6.on_off << endl;
                  txtfile << "Name: " << Sensor_prm_CSS6.Name << endl;
                  txtfile << "constparam: " << Sensor_prm_CSS6.ConstPrm(0) << " , " << Sensor_prm_CSS6.ConstPrm(1) << endl;
                  txtfile << "Mounting matrix: " << Sensor_prm_CSS6.SC2SYS(0,0) << " , " << Sensor_prm_CSS6.SC2SYS(0,1) << " , " << Sensor_prm_CSS6.SC2SYS(0,2) << " , " << endl;
                  txtfile << "                 " << Sensor_prm_CSS6.SC2SYS(1,0) << " , " << Sensor_prm_CSS6.SC2SYS(1,1) << " , " << Sensor_prm_CSS6.SC2SYS(1,2) << " , " << endl;
                  txtfile << "                 " << Sensor_prm_CSS6.SC2SYS(2,0) << " , " << Sensor_prm_CSS6.SC2SYS(2,1) << " , " << Sensor_prm_CSS6.SC2SYS(2,2) << " , " << "\n" << endl;
                  txtfile << "accuracy: " << Sensor_prm_CSS6.Accuracy(0) << " , " << Sensor_prm_CSS6.Accuracy(1) << "\n" << endl;
                  
                  txtfile << "Magnetometer (deployed): " << endl;
                  txtfile << "On/off: " << Sensor_prm_MAG.on_off << endl;
                  txtfile << "Name: " << Sensor_prm_MAG.Name << endl;
                  txtfile << "Mounting matrix: " << Sensor_prm_MAG.SC2SYS(0,0) << " , " << Sensor_prm_MAG.SC2SYS(0,1) << " , " << Sensor_prm_MAG.SC2SYS(0,2) << " , " << endl;
                  txtfile << "                 " << Sensor_prm_MAG.SC2SYS(1,0) << " , " << Sensor_prm_MAG.SC2SYS(1,1) << " , " << Sensor_prm_MAG.SC2SYS(1,2) << " , " << endl;
                  txtfile << "                 " << Sensor_prm_MAG.SC2SYS(2,0) << " , " << Sensor_prm_MAG.SC2SYS(2,1) << " , " << Sensor_prm_MAG.SC2SYS(2,2) << " , " << "\n" << endl;
                  txtfile << "accuracy: " << Sensor_prm_MAG.Accuracy(0) << " , " << Sensor_prm_MAG.Accuracy(1) << "\n" << endl;
                  
                  txtfile << "Magnetometer (stowed): " << endl;
                  txtfile << "On/off: " << Sensor_prm_MAGstowed.on_off << endl;
                  txtfile << "Name: " << Sensor_prm_MAGstowed.Name << endl;
                  txtfile << "Mounting matrix: " << Sensor_prm_MAGstowed.SC2SYS(0,0) << " , " << Sensor_prm_MAGstowed.SC2SYS(0,1) << " , " << Sensor_prm_MAGstowed.SC2SYS(0,2) << " , " << endl;
                  txtfile << "                 " << Sensor_prm_MAGstowed.SC2SYS(1,0) << " , " << Sensor_prm_MAGstowed.SC2SYS(1,1) << " , " << Sensor_prm_MAGstowed.SC2SYS(1,2) << " , " << endl;
                  txtfile << "                 " << Sensor_prm_MAGstowed.SC2SYS(2,0) << " , " << Sensor_prm_MAGstowed.SC2SYS(2,1) << " , " << Sensor_prm_MAGstowed.SC2SYS(2,2) << " , " << "\n" << endl;
                  txtfile << "accuracy: " << Sensor_prm_MAGstowed.Accuracy(0) << " , " << Sensor_prm_MAGstowed.Accuracy(1) << "\n" << endl;
                  
                  txtfile << "Rate sensor: " << endl;
                  txtfile << "On/off: " << Sensor_prm_RS.on_off << endl;
                  txtfile << "Name: " << Sensor_prm_RS.Name << endl;
                  txtfile << "Mounting matrix: " << Sensor_prm_RS.SC2SYS(0,0) << " , " << Sensor_prm_RS.SC2SYS(0,1) << " , " << Sensor_prm_RS.SC2SYS(0,2) << " , " << endl;
                  txtfile << "                 " << Sensor_prm_RS.SC2SYS(1,0) << " , " << Sensor_prm_RS.SC2SYS(1,1) << " , " << Sensor_prm_RS.SC2SYS(1,2) << " , " << endl;
                  txtfile << "                 " << Sensor_prm_RS.SC2SYS(2,0) << " , " << Sensor_prm_RS.SC2SYS(2,1) << " , " << Sensor_prm_RS.SC2SYS(2,2) << " , " << "\n" << endl;
                  txtfile << "accuracy: " << Sensor_prm_RS.Accuracy(0) << " , " << Sensor_prm_RS.Accuracy(1) << "\n" << endl;
                  
                  txtfile << "Magnetic torquer: " << endl;
                  txtfile << "On/off: " << Sensor_prm_MAGTRQ.on_off << endl;
                  txtfile << "Name: " << Sensor_prm_MAGTRQ.Name << endl;
                  txtfile << "constparam: " << Sensor_prm_MAGTRQ.ConstPrm(0) << endl;
                  txtfile << "Mounting matrix: " << Sensor_prm_MAGTRQ.SC2SYS(0,0) << " , " << Sensor_prm_MAGTRQ.SC2SYS(0,1) << " , " << Sensor_prm_MAGTRQ.SC2SYS(0,2) << " , " << endl;
                  txtfile << "                 " << Sensor_prm_MAGTRQ.SC2SYS(1,0) << " , " << Sensor_prm_MAGTRQ.SC2SYS(1,1) << " , " << Sensor_prm_MAGTRQ.SC2SYS(1,2) << " , " << endl;
                  txtfile << "                 " << Sensor_prm_MAGTRQ.SC2SYS(2,0) << " , " << Sensor_prm_MAGTRQ.SC2SYS(2,1) << " , " << Sensor_prm_MAGTRQ.SC2SYS(2,2) << " , " << "\n" << endl;
                  txtfile << "accuracy: " << Sensor_prm_MAGTRQ.Accuracy(0) << " , " << Sensor_prm_MAGTRQ.Accuracy(1) << "\n" << endl;
                  
                  txtfile << "Wheel+X: " << endl;
                  txtfile << "On/off: " << Sensor_prm_WHEEL1.on_off << endl;
                  txtfile << "Name: " << Sensor_prm_WHEEL1.Name << endl;
                  txtfile << "constparam: " << Sensor_prm_WHEEL1.ConstPrm(0) << endl;
                  txtfile << "auxparam: " << Sensor_prm_WHEEL1.AuxPrm(0) << endl;
                  txtfile << "Mounting matrix: " << Sensor_prm_WHEEL1.SC2SYS(0,0) << " , " << Sensor_prm_WHEEL1.SC2SYS(0,1) << " , " << Sensor_prm_WHEEL1.SC2SYS(0,2) << " , " << endl;
                  txtfile << "                 " << Sensor_prm_WHEEL1.SC2SYS(1,0) << " , " << Sensor_prm_WHEEL1.SC2SYS(1,1) << " , " << Sensor_prm_WHEEL1.SC2SYS(1,2) << " , " << endl;
                  txtfile << "                 " << Sensor_prm_WHEEL1.SC2SYS(2,0) << " , " << Sensor_prm_WHEEL1.SC2SYS(2,1) << " , " << Sensor_prm_WHEEL1.SC2SYS(2,2) << " , " << "\n" << endl;
                  txtfile << "accuracy: " << Sensor_prm_WHEEL1.Accuracy(0) << " , " << Sensor_prm_WHEEL1.Accuracy(1) << " , " << Sensor_prm_WHEEL1.Accuracy(2) << " , " << Sensor_prm_WHEEL1.Accuracy(3) << "\n" << endl;
                  
                  txtfile << "Wheel+Y: " << endl;
                  txtfile << "On/off: " << Sensor_prm_WHEEL2.on_off << endl;
                  txtfile << "Name: " << Sensor_prm_WHEEL2.Name << endl;
                  txtfile << "constparam: " << Sensor_prm_WHEEL2.ConstPrm(0) << endl;
                  txtfile << "auxparam: " << Sensor_prm_WHEEL2.AuxPrm(0) << endl;
                  txtfile << "Mounting matrix: " << Sensor_prm_WHEEL2.SC2SYS(0,0) << " , " << Sensor_prm_WHEEL2.SC2SYS(0,1) << " , " << Sensor_prm_WHEEL2.SC2SYS(0,2) << " , " << endl;
                  txtfile << "                 " << Sensor_prm_WHEEL2.SC2SYS(1,0) << " , " << Sensor_prm_WHEEL2.SC2SYS(1,1) << " , " << Sensor_prm_WHEEL2.SC2SYS(1,2) << " , " << endl;
                  txtfile << "                 " << Sensor_prm_WHEEL2.SC2SYS(2,0) << " , " << Sensor_prm_WHEEL2.SC2SYS(2,1) << " , " << Sensor_prm_WHEEL2.SC2SYS(2,2) << " , " << "\n" << endl;
                  txtfile << "accuracy: " << Sensor_prm_WHEEL2.Accuracy(0) << " , " << Sensor_prm_WHEEL2.Accuracy(1) << " , " << Sensor_prm_WHEEL2.Accuracy(2) << " , " << Sensor_prm_WHEEL2.Accuracy(3) << "\n" << endl;
                  
                  txtfile << "Wheel+Z: " << endl;
                  txtfile << "On/off: " << Sensor_prm_WHEEL3.on_off << endl;
                  txtfile << "Name: " << Sensor_prm_WHEEL3.Name << endl;
                  txtfile << "constparam: " << Sensor_prm_WHEEL3.ConstPrm(0) << endl;
                  txtfile << "auxparam: " << Sensor_prm_WHEEL3.AuxPrm(0) << endl;
                  txtfile << "Mounting matrix: " << Sensor_prm_WHEEL3.SC2SYS(0,0) << " , " << Sensor_prm_WHEEL3.SC2SYS(0,1) << " , " << Sensor_prm_WHEEL3.SC2SYS(0,2) << " , " << endl;
                  txtfile << "                 " << Sensor_prm_WHEEL3.SC2SYS(1,0) << " , " << Sensor_prm_WHEEL3.SC2SYS(1,1) << " , " << Sensor_prm_WHEEL3.SC2SYS(1,2) << " , " << endl;
                  txtfile << "                 " << Sensor_prm_WHEEL3.SC2SYS(2,0) << " , " << Sensor_prm_WHEEL3.SC2SYS(2,1) << " , " << Sensor_prm_WHEEL3.SC2SYS(2,2) << " , " << "\n" << endl;
                  txtfile << "accuracy: " << Sensor_prm_WHEEL3.Accuracy(0) << " , " << Sensor_prm_WHEEL3.Accuracy(1) << " , " << Sensor_prm_WHEEL3.Accuracy(2) << " , " << Sensor_prm_WHEEL3.Accuracy(3) << "\n" << endl;
                  
                  txtfile << "\n######################" << endl;
                  txtfile << "SOLAR PANELS" << endl;
                  txtfile << "######################\n" << endl;
                  
                  txtfile << "On/off: " << Solarpan1_prm.on_off << endl;
                  txtfile << "Name: " << Solarpan1_prm.Name << endl;
                  txtfile << "constparam: " << Solarpan1_prm.ConstPrm(0) << " , " << Solarpan1_prm.ConstPrm(1) << endl;
                  txtfile << "auxparam: " << Solarpan1_prm.AuxPrm(0) << endl;
                  txtfile << "Mounting matrix: " << Solarpan1_prm.SC2SYS(0,0) << " , " << Solarpan1_prm.SC2SYS(0,1) << " , " << Solarpan1_prm.SC2SYS(0,2) << " , " << endl;
                  txtfile << "                 " << Solarpan1_prm.SC2SYS(1,0) << " , " << Solarpan1_prm.SC2SYS(1,1) << " , " << Solarpan1_prm.SC2SYS(1,2) << " , " << endl;
                  txtfile << "                 " << Solarpan1_prm.SC2SYS(2,0) << " , " << Solarpan1_prm.SC2SYS(2,1) << " , " << Solarpan1_prm.SC2SYS(2,2) << " , " << "\n" << endl;
                  txtfile << "OPS_limits: " << Solarpan1_prm.OPS_limits(0) << endl;
                  
                  txtfile << "On/off: " << Solarpan2_prm.on_off << endl;
                  txtfile << "Name: " << Solarpan2_prm.Name << endl;
                  txtfile << "constparam: " << Solarpan2_prm.ConstPrm(0) << " , " << Solarpan2_prm.ConstPrm(1) << endl;
                  txtfile << "auxparam: " << Solarpan2_prm.AuxPrm(0) << endl;
                  txtfile << "Mounting matrix: " << Solarpan2_prm.SC2SYS(0,0) << " , " << Solarpan2_prm.SC2SYS(0,1) << " , " << Solarpan2_prm.SC2SYS(0,2) << " , " << endl;
                  txtfile << "                 " << Solarpan2_prm.SC2SYS(1,0) << " , " << Solarpan2_prm.SC2SYS(1,1) << " , " << Solarpan2_prm.SC2SYS(1,2) << " , " << endl;
                  txtfile << "                 " << Solarpan2_prm.SC2SYS(2,0) << " , " << Solarpan2_prm.SC2SYS(2,1) << " , " << Solarpan2_prm.SC2SYS(2,2) << " , " << "\n" << endl;
                  txtfile << "OPS_limits: " << Solarpan2_prm.OPS_limits(0) << endl;
                  
                  txtfile << "On/off: " << Solarpan3_prm.on_off << endl;
                  txtfile << "Name: " << Solarpan3_prm.Name << endl;
                  txtfile << "constparam: " << Solarpan3_prm.ConstPrm(0) << " , " << Solarpan3_prm.ConstPrm(1) << endl;
                  txtfile << "auxparam: " << Solarpan3_prm.AuxPrm(0) << endl;
                  txtfile << "Mounting matrix: " << Solarpan3_prm.SC2SYS(0,0) << " , " << Solarpan3_prm.SC2SYS(0,1) << " , " << Solarpan3_prm.SC2SYS(0,2) << " , " << endl;
                  txtfile << "                 " << Solarpan3_prm.SC2SYS(1,0) << " , " << Solarpan3_prm.SC2SYS(1,1) << " , " << Solarpan3_prm.SC2SYS(1,2) << " , " << endl;
                  txtfile << "                 " << Solarpan3_prm.SC2SYS(2,0) << " , " << Solarpan3_prm.SC2SYS(2,1) << " , " << Solarpan3_prm.SC2SYS(2,2) << " , " << "\n" << endl;
                  txtfile << "OPS_limits: " << Solarpan3_prm.OPS_limits(0) << endl;
                  
                  txtfile << "\n######################" << endl;
                  txtfile << "PROPULSION SYSTEMS" << endl;
                  txtfile << "######################\n" << endl;
                  
                  txtfile << "Orbit propulsion system: " << endl;
                  txtfile << "On/off: " << OrbitPropulsion1_prm.on_off << endl;
                  txtfile << "Name: " << OrbitPropulsion1_prm.Name << endl;
                  txtfile << "Thrust: " << OrbitPropulsion1_prm.ConstPrm(0) << endl;
                  txtfile << "Thrust resolution: " << OrbitPropulsion1_prm.ConstPrm(1) << endl;
                  txtfile << "Maximum available thrust: " << OrbitPropulsion1_prm.OPS_limits(0) << endl;
                  txtfile << "dv accuracy: " << OrbitPropulsion1_prm.Accuracy(0) << " , " << OrbitPropulsion1_prm.Accuracy(1) << "\n" << endl;
                  txtfile << "Attitude accuracy: " << OrbitPropulsion1_prm.Accuracy(2) << " , " << OrbitPropulsion1_prm.Accuracy(3) << "\n" << endl;
                  txtfile << "Mounting matrix: " << OrbitPropulsion1_prm.SC2SYS(0,0) << " , " << OrbitPropulsion1_prm.SC2SYS(0,1) << " , " << OrbitPropulsion1_prm.SC2SYS(0,2) << " , " << endl;
                  txtfile << "                 " << OrbitPropulsion1_prm.SC2SYS(1,0) << " , " << OrbitPropulsion1_prm.SC2SYS(1,1) << " , " << OrbitPropulsion1_prm.SC2SYS(1,2) << " , " << endl;
                  txtfile << "                 " << OrbitPropulsion1_prm.SC2SYS(2,0) << " , " << OrbitPropulsion1_prm.SC2SYS(2,1) << " , " << OrbitPropulsion1_prm.SC2SYS(2,2) << " , " << "\n" << endl;
                  
                  txtfile << "Orbit propulsion system: " << endl;
                  txtfile << "On/off: " << OrbitPropulsion2_prm.on_off << endl;
                  txtfile << "Name: " << OrbitPropulsion2_prm.Name << endl;
                  txtfile << "Thrust: " << OrbitPropulsion2_prm.ConstPrm(0) << endl;
                  txtfile << "Thrust resolution: " << OrbitPropulsion2_prm.ConstPrm(1) << endl;
                  txtfile << "Maximum available thrust: " << OrbitPropulsion2_prm.OPS_limits(0) << endl;
                  txtfile << "dv accuracy: " << OrbitPropulsion2_prm.Accuracy(0) << " , " << OrbitPropulsion2_prm.Accuracy(1) << "\n" << endl;
                  txtfile << "Attitude accuracy: " << OrbitPropulsion2_prm.Accuracy(2) << " , " << OrbitPropulsion2_prm.Accuracy(3) << "\n" << endl;
                  txtfile << "Mounting matrix: " << OrbitPropulsion2_prm.SC2SYS(0,0) << " , " << OrbitPropulsion2_prm.SC2SYS(0,1) << " , " << OrbitPropulsion2_prm.SC2SYS(0,2) << " , " << endl;
                  txtfile << "                 " << OrbitPropulsion2_prm.SC2SYS(1,0) << " , " << OrbitPropulsion2_prm.SC2SYS(1,1) << " , " << OrbitPropulsion2_prm.SC2SYS(1,2) << " , " << endl;
                  txtfile << "                 " << OrbitPropulsion2_prm.SC2SYS(2,0) << " , " << OrbitPropulsion2_prm.SC2SYS(2,1) << " , " << OrbitPropulsion2_prm.SC2SYS(2,2) << " , " << "\n" << endl;
                  
                  txtfile << "\n######################" << endl;
                  txtfile << "MANEUVERS" << endl;
                  txtfile << "######################\n" << endl;
                  
                  for(int i = 0; i < maneuvers_num; i++)
                      {
                      txtfile << all_maneuvers[i].name << ": " << endl;
                      txtfile << "maneuver_on: " << all_maneuvers[i].maneuver_on << endl;
                      txtfile << "init_time: " << all_maneuvers[i].init_time << endl;
                      txtfile << "duration: " << all_maneuvers[i].duration << endl;
                      txtfile << "ManVec: " << all_maneuvers[i].ManVec(0) << " , " << all_maneuvers[i].ManVec(1) << " , " << all_maneuvers[i].ManVec(2) << "\n" << endl;
                      }
                  
                  txtfile << fixed << endl;
                  
                  txtfile.close();
                  };
//-------------------------------------------------------------------------------------
// int XML_parser_events(...)
//-------------------------------------------------------------------------------------
/**
 * Load and read XML file for events computation and validate it against an XML schema file
 *
 * @param XML_events_file     Complete path to XML file (e.g. dir1/dir2/dir3/filenae.xml)
 * @param simstep etc         Variables in which read computation parameters are put
 *
 * @return exception against XML schema
 *
 * @note The implementation of this parser is based on CodeSynthesis XSD binding compiler and the xerces library
 * @see http://www.codesynthesis.com/products/xsd/
 * @see https://xerces.apache.org/xerces-c/
 */
//-------------------------------------------------------------------------------------
int XML_parser_events(const string XML_events_file,
                      int& simstep,
                      int& duration,
                      double& FOV_cross,
                      double& FOV_along,
                      int& SC_start,
                      int& SC_end,
                      int& PL_start,
                      int& PL_end,
                      bool& TGs_on,
                      bool& GSs_on,
                      bool& TGs_grid_on,
                      bool& Eclipse_on,
                      Vec4d& TG_grid_limits,
                      double& gridstep,
                      ground::TG* TGs_list,
                      ground::GS* GSs_list,
                      string& Orbit_ephemeris_path,
                      string& Orbit_ephemeris_rootname,
                      string& Data_path,
                      string& planetephemeris,
                      string& eop,
                      string& pck_data,
                      string& leapsecond,
                      string& TG_filename,
                      string& GS_filename,
                      string& Eclipse_filename)
                      {
                      // Instantiate individual parsers.
                      //
                      ::eventsparam_pimpl eventsparam_p;
                      ::fileheader_pimpl fileheader_p;
                      ::xml_schema::string_pimpl string_p;
                      ::license_pimpl license_p;
                      ::xml_schema::uri_pimpl uri_p;
                      ::xml_schema::date_pimpl date_p;
                      ::reference_pimpl reference_p;
                      ::CompParameters_pimpl CompParameters_p;
                      ::durstep_pimpl durstep_p;
                      ::xml_schema::duration_pimpl duration_p;
                      ::Payload_pimpl Payload_p;
                      ::Angle_pimpl Angle_p;
                      ::AngleType_pimpl AngleType_p;
                      ::Spacecraft_pimpl Spacecraft_p;
                      ::xml_schema::positive_integer_pimpl positive_integer_p;
                      ::Compoptions_pimpl Compoptions_p;
                      ::xml_schema::boolean_pimpl boolean_p;
                      ::TGs_pimpl TGs_p;
                      ::TGs_grid_pimpl TGs_grid_p;
                      ::TGs_list_pimpl TGs_list_p;
                      ::TG_pimpl TG_p;
                      ::Altitude_pimpl Altitude_p;
                      ::LengthType_pimpl LengthType_p;
                      ::GSs_pimpl GSs_p;
                      ::GS_pimpl GS_p;
                      ::EventsInputFiles_pimpl EventsInputFiles_p;
                      ::EventsOutputFiles_pimpl EventsOutputFiles_p;
                  
                      // Connect the parsers together.
                      //
                      eventsparam_p.parsers(fileheader_p,
                                            CompParameters_p,
                                            TGs_p,
                                            GSs_p,
                                            EventsInputFiles_p,
                                            EventsOutputFiles_p,
                                            string_p);
                      
                      fileheader_p.parsers (string_p,
                                            string_p,
                                            string_p,
                                            license_p,
                                            string_p,
                                            date_p,
                                            string_p,
                                            string_p,
                                            string_p,
                                            string_p,
                                            reference_p);
                  
                      license_p.parsers (string_p,
                                         uri_p);
                  
                      reference_p.parsers (string_p,
                                           string_p,
                                           string_p,
                                           string_p);
                  
                      CompParameters_p.parsers (durstep_p,
                                                Payload_p,
                                                Spacecraft_p,
                                                Compoptions_p);
                  
                      durstep_p.parsers (duration_p,
                                         duration_p);
                  
                      Payload_p.parsers (Angle_p,
                                         Angle_p);
                  
                      Angle_p.parsers (AngleType_p);
                  
                      Spacecraft_p.parsers (positive_integer_p,
                                            positive_integer_p,
                                            positive_integer_p,
                                            positive_integer_p);
                  
                      Compoptions_p.parsers (boolean_p,
                                             boolean_p,
                                             boolean_p,
                                             boolean_p);
                  
                      TGs_p.parsers (TGs_grid_p,
                                     TGs_list_p);
                  
                      TGs_grid_p.parsers (Angle_p,
                                          Angle_p,
                                          Angle_p,
                                          Angle_p,
                                          Angle_p);
                  
                      TGs_list_p.parsers (TG_p);
                  
                      TG_p.parsers (Angle_p,
                                    Angle_p,
                                    Altitude_p,
                                    string_p);
                  
                      Altitude_p.parsers (LengthType_p);
                  
                      GSs_p.parsers (GS_p);
                  
                      GS_p.parsers (Angle_p,
                                    Angle_p,
                                    Altitude_p,
                                    Angle_p,
                                    string_p);
                  
                      EventsInputFiles_p.parsers (string_p,
                                            string_p,
                                            string_p,
                                            string_p,
                                            string_p,
                                            string_p,
                                            string_p,
                                            string_p);
                  
                      EventsOutputFiles_p.parsers (string_p,
                                             string_p,
                                             string_p,
                                             string_p);
                        
                      try
                        {
                        // Parse the XML document.
                        //
                        ::xml_schema::document doc_p(eventsparam_p, "eventsparam");
                        
                        eventsparam_p.pre ();
                        doc_p.parse(XML_events_file);
                        eventsparam_p.post_eventsparam ();
                        }
                      catch (const ::xml_schema::exception& e)
                        {
                        std::cerr << e << std::endl;
                        return 1;
                        }
                        
                      //////////////////////////////////////////////////////////////////////////
                      ///////////////////////// COMPUTATION PARAMETERS /////////////////////////
                      //////////////////////////////////////////////////////////////////////////
                      // Step of input ephemerides
                      simstep = durstep_p.sim_step;
                      // Time span considered for the events computation
                      duration = durstep_p.sim_duration;
                      // Field of view half angle perpendicolar to the direction of motion
                      FOV_cross = Payload_p.FOV_cross_in;
                      // Field of view half angle along the direction of motion
                      FOV_along = Payload_p.FOV_along_in;
                      // Spacecraft number to start with
                      SC_start = Spacecraft_p.SC_start_in;
                      // Spacecraft number to end with
                      SC_end = Spacecraft_p.SC_end_in;
                      // Orbital plane number to start with
                      PL_start = Spacecraft_p.PL_start_in;
                      // Orbital plane number to end with
                      PL_end = Spacecraft_p.PL_end_in;
                      // Compute contacts of payload with targets
                      TGs_on = Compoptions_p.TGs_on_in;
                      // Compute contacts of spacecraft with ground stations
                      GSs_on = Compoptions_p.GSs_on_in;
                      // Use targets grid for computation of contacts of payload with targets
                      TGs_grid_on = Compoptions_p.TGs_grid_on_in;
                      // Compute eclipse times
                      Eclipse_on = Compoptions_p.Eclipse_on_in;
                      ///////////////////////////////////////////////////////////
                      ///////////////////////// TARGETS /////////////////////////
                      ///////////////////////////////////////////////////////////
                      // Minimum and maximum longitute and latitude of targets grid
                      TG_grid_limits(0) = TGs_grid_p.minlon_in;
                      TG_grid_limits(1) = TGs_grid_p.maxlon_in;
                      TG_grid_limits(2) = TGs_grid_p.minlat_in;
                      TG_grid_limits(3) = TGs_grid_p.maxlat_in;
                      // Step of the targets grid
                      gridstep = TGs_grid_p.gridstep_in;
                      // Targets list
                      unsigned int ind = 0;
                      //size_t listsize = sizeof(TGs_list_p.TGs_list_in)/sizeof(TGs_list_p.TGs_list_in[0]);
                      unsigned int listsize = 1000;
                      
                      while( !TGs_list_p.TGs_list_in[ind].name.empty() && ind < listsize )
                            {
                            TGs_list[ind].name = TGs_list_p.TGs_list_in[ind].name;  
                            TGs_list[ind].lon = TGs_list_p.TGs_list_in[ind].lon;  
                            TGs_list[ind].lat = TGs_list_p.TGs_list_in[ind].lat;  
                            TGs_list[ind].alt = TGs_list_p.TGs_list_in[ind].alt;
                              
                            ind++;
                            }
                      //////////////////////////////////////////////////////////////////////////
                      ///////////////////////// GROUND STATIONS /////////////////////////
                      //////////////////////////////////////////////////////////////////////////
                      // Ground stations list
                      ind = 0;
                      //listsize = sizeof(GSs_p.GSs_list_in)/sizeof(GSs_p.GSs_list_in[0]);
                      
                      while( !GSs_p.GSs_list_in[ind].name.empty() && ind < listsize )
                            {
                            GSs_list[ind].name = GSs_p.GSs_list_in[ind].name;  
                            GSs_list[ind].lon = GSs_p.GSs_list_in[ind].lon;  
                            GSs_list[ind].lat = GSs_p.GSs_list_in[ind].lat;  
                            GSs_list[ind].alt = GSs_p.GSs_list_in[ind].alt;
                            GSs_list[ind].minelev = GSs_p.GSs_list_in[ind].minelev;
                              
                            ind++;
                            }
                      ///////////////////////////////////////////////////////////////
                      ///////////////////////// FILES PATHS /////////////////////////
                      ///////////////////////////////////////////////////////////////
                      // Input files paths
                      Orbit_ephemeris_path = EventsInputFiles_p.Orbit_ephemeris_path_in;
                      Orbit_ephemeris_rootname = EventsInputFiles_p.Orbit_ephemeris_rootname_in;
                      Data_path = EventsInputFiles_p.Data_path_in;
                      planetephemeris = EventsInputFiles_p.planetephemeris_in;
                      pck_data = EventsInputFiles_p.pck_data_in;
                      eop = EventsInputFiles_p.eop_in;
                      leapsecond = EventsInputFiles_p.leapsecond_in;
                      // Output files paths
                      TG_filename = EventsOutputFiles_p.TG_contacts_in;
                      GS_filename = EventsOutputFiles_p.GS_contacts_in;
                      Eclipse_filename = EventsOutputFiles_p.Eclipse_times_in;
                        
                      return(0);  
                      }
//-------------------------------------------------------------------------------------
// void ReadXMLeventstoTXT(...)
//-------------------------------------------------------------------------------------
/**
 * Write events computation parameters read by parser XML_parser_events in a text file for verification purposes
 *
 * @param txt_file        Complete path to txt file (e.g. dir1/dir2/dir3/filenae.txt)
 * @param simstep etc     Computation parameters
 *
 */
//------------------------------------------------------------------------------------- 
void ReadXMLeventstoTXT(const string txt_file,
                        int simstep,
                        int duration,
                        double FOV_cross,
                        double FOV_along,
                        int SC_start,
                        int SC_end,
                        int PL_start,
                        int PL_end,
                        bool TGs_on,
                        bool GSs_on,
                        bool TGs_grid_on,
                        bool Eclipse_on,
                        Vec4d TG_grid_limits,
                        double gridstep,
                        ground::TG* TGs_list,
                        ground::GS* GSs_list,
                        string Orbit_ephemeris_path,
                        string Orbit_ephemeris_rootname,
                        string Data_path,
                        string planetephemeris,
                        string eop,
                        string pck_data,
                        string leapsecond,
                        string TG_filename,
                        string GS_filename,
                        string Eclipse_filename)
                        {
                        ofstream txtfile;
                        txtfile.open(txt_file);
                        
                        txtfile << "\n######################" << endl;
                        txtfile << "COMPUTATION PARAMETERS" << endl;
                        txtfile << "######################\n" << endl;
                        txtfile << "simstep: " << simstep << endl;
                        txtfile << "duration: " << duration << "\n" << endl;
                        txtfile << "FOV_cross: " << FOV_cross << endl;
                        txtfile << "FOV_along: " << FOV_along << "\n" << endl;
                        txtfile << "SC_start: " << SC_start << endl;
                        txtfile << "SC_end: " << SC_end << endl;
                        txtfile << "PL_start: " << PL_start << endl;
                        txtfile << "PL_end: " << PL_end << "\n" << endl;
                        txtfile << "TGs_on: " << TGs_on << endl;
                        txtfile << "GSs_on: " << GSs_on << endl;
                        txtfile << "TGs_grid_on: " << TGs_grid_on << endl;
                        txtfile << "Eclipse_on: " << Eclipse_on << "\n" << endl;
                        txtfile << "\n######################" << endl;
                        txtfile << "TARGETS" << endl;
                        txtfile << "######################\n" << endl;
                        txtfile << "TG_grid_limits(0): " << TG_grid_limits(0) << endl;
                        txtfile << "TG_grid_limits(1): " << TG_grid_limits(1) << endl;
                        txtfile << "TG_grid_limits(2): " << TG_grid_limits(2) << endl;
                        txtfile << "TG_grid_limits(3): " << TG_grid_limits(3) << endl;
                        txtfile << "gridstep: " << gridstep << "\n" << endl;
                        unsigned int ind = 0;
                        unsigned int listsize = 1000;
                        //cout << listsize << endl;
                        while( !TGs_list[ind].name.empty() && ind < listsize )
                              {
                              txtfile << "TG name: " << TGs_list[ind].name << endl;
                              txtfile << "lon: " << TGs_list[ind].lon << endl;
                              txtfile << "lat: " << TGs_list[ind].lat << endl;
                              txtfile << "alt: " << TGs_list[ind].alt << "\n" << endl;
                                
                              ind++;
                              }
                        txtfile << "\n######################" << endl;
                        txtfile << "GROUND STATIONS" << endl;
                        txtfile << "######################\n" << endl;
                        ind = 0;
                        listsize = 1000;
                        
                        while( !GSs_list[ind].name.empty() && ind < listsize )
                              {
                              txtfile << "GS name: " << GSs_list[ind].name << endl;
                              txtfile << "lon: " << GSs_list[ind].lon << endl;
                              txtfile << "lat: " << GSs_list[ind].lat << endl;
                              txtfile << "alt: " << GSs_list[ind].alt << endl;
                              txtfile << "minelev: " << GSs_list[ind].minelev << "\n" << endl;
                                
                              ind++;
                              }
                        txtfile << "######################" << endl;
                        txtfile << "INPUT FILES PATHS" << endl;
                        txtfile << "######################\n" << endl;
                        txtfile << "Orbit Orbit_ephemeris_path: " << Orbit_ephemeris_path << endl;
                        txtfile << "Orbit_ephemeris_rootname: " << Orbit_ephemeris_rootname << endl;
                        txtfile << "Data path: " << Data_path << endl;
                        txtfile << "Planet ephemeris: " << planetephemeris << endl;
                        txtfile << "EOP: " << eop << endl;
                        txtfile << "PCK: " << pck_data << endl;
                        txtfile << "Leap second: " << leapsecond << endl;
                        txtfile << "\n######################" << endl;
                        txtfile << "OUTPUT FILES PATHS" << endl;
                        txtfile << "######################\n" << endl;
                        txtfile << "TG_filename: " << TG_filename << endl;
                        txtfile << "GS_filename: " << GS_filename << endl;
                        txtfile << "Eclipse_filename: " << Eclipse_filename << endl;
                          
                        txtfile << fixed << endl;
                  
                        txtfile.close();
                        }
//-------------------------------------------------------------------------------------
// void RunStatusBar(double t, int simduration, int barwidth)
//-------------------------------------------------------------------------------------
/**
 * Display in terminal the simulation execution status bar
 *
 * @param t             Current simulation time [s]
 * @param simduration   Duration of simulation run [s]
 * @param barwidth      Desired width of status bar at the end of the simulation (100 %)
 *
 */
//------------------------------------------------------------------------------------- 
void RunStatusBar(double t,
                  int simduration,
                  int barwidth)
                  {
                  int barpos = 0;
                  static bool barinit;
                  static double part_dur;
                  static int sim_done;
                  
                  if(!barinit)
                    {
                    part_dur = simduration/10.0;
                    sim_done = 10;
                    barinit = true;
                    }
                  
                  if( (t - part_dur) >= 0)
                    {
                    barpos = barwidth*sim_done/100;
                    for (int i = 0; i < barwidth; ++i) if(i <= barpos) cout << "\u25A0";
                    cout << " " << sim_done << "%" << endl;
                    
                    part_dur = part_dur + simduration/10.0;
                    sim_done = sim_done + 10.0;
                    }  
                    
                    
                    
                    
                    
                    
                  };
                  






//Eigen::MatrixXd read_matfile(const char* filename,
//                             const char* matvarname)
//        {
//        mat_t *mateph;
//        mateph = Mat_Open(filename,MAT_ACC_RDONLY);
//        if(NULL == mateph)
//          {
//          fprintf(stderr,"Error opening MAT file!\n");
//          //return EXIT_FAILURE;
//          }
//          
//        matvar_t *ephvar;
//        
//        ephvar = Mat_VarReadInfo(mateph,matvarname);
//        if(NULL == ephvar) fprintf(stderr,"Variable stateEME_MAIN not found, or error reading MAT file\n");
//        cout << "Dimensions:" << ephvar->dims[0] << " and " << ephvar->dims[1] << endl;
//        
//        int rows = (int)ephvar->dims[0];
//        int cols = (int)ephvar->dims[1];
//        
//        double data[rows*cols];
//        int start[2] = {0, 0};
//        int stride[2] = {1, 1};
//        int edge[2] = {rows, cols};
//      
//        Mat_VarReadData(mateph,ephvar,data,start,stride,edge);
//        
//        //Eigen::Matrix<double, rows, cols> ephem;
//        Eigen::MatrixXd ephem;
//        ephem.resize(rows, cols);
//        
//        
//        for(int k = 0; k < rows; k++)
//            {
//            ephem(k,0) = data[k]; ephem(k,1) = data[k+rows]; ephem(k,2) = data[k+2*rows]; ephem(k,3) = data[k+3*rows]; ephem(k,4) = data[k+4*rows]; ephem(k,5) = data[k+5*rows]; ephem(k,6) = data[k+6*rows]; ephem(k,7) = data[k+7*rows];
//            }
//        
//        cout << "Matrix dimensions:" << ephem.rows() << " and " << ephem.cols() << endl;    
//        //Matrix ephem(5,5);    
//        //
//        //for(int i=0; i< 5; i++)
//        //   for(int j=0; j< cols; j++) ephem(i,j) = temp[i][j];
//        
//         
//        //for(int i=0; i< 5; i++)
//        // {
//        //  for(int j=0; j< cols; j++) cout << temp[i][j] << "   ";
//        //  cout << "\n" << endl;
//        // } 
//        
//        Mat_VarFree(ephvar);
//        Mat_Close(mateph);
//        
//        return ephem;
//        };