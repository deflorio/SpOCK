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

#include <simparam_schema-pimpl.h>

#include <iostream>

#include <VarTypes.h>

int simparam_pimpl::sensact_ind = 0;

double positive_num_value;
double angle_value;
int nMAX_value;
double x_value, y_value, z_value;
double vx_value, vy_value, vz_value;
double face_area;
string face_material;
Vec3d face_cP, face_cA;
double m11_in, m12_in, m13_in, m21_in, m22_in, m23_in, m31_in, m32_in, m33_in;
bool subsystem_on_in;

Mat3x3d SC2SYS_in;
string name_in;

double constpar;
double auxpar;
double opslim;
double accur;

Eigen::VectorXd constparam_in;
Eigen::VectorXd auxparam_in;
Eigen::VectorXd opslimits_in;
Eigen::VectorXd accuracy_in;

//Eigen::Matrix< double, 10, 1 > constparam_in;
//Eigen::Matrix< double, 10, 1 > auxparam_in;
//Eigen::Matrix< double, 10, 1 > opslimits_in;
//Eigen::Matrix< double, 10, 1 > accuracy_in;

string man_name_in = " ";
bool maneuver_on_in = false;
double man_init_time_in= 0.0;
double man_duration_in = 0.0;
Vec3d ManVec_in;
//ManVec_in = Eigen::Vec3d::Zero();

SC::maneuver maneuver_in;


//maneuver_temp.name = man_name_in;
//maneuver_temp.maneuver_on = maneuver_on_in;
//maneuver_temp.init_time = man_init_time_in;
//maneuver_temp.duration = man_duration_in;
//maneuver_temp.ManVec = ManVec_in;

//SC::maneuver maneuver_temp;

//SC::maneuver Maneuvers_pimpl::maneuver_in = maneuver_temp;
//
//vector<SC::maneuver> Maneuvers_pimpl::all_maneuvers = {maneuver_in};




double alt_value;

double TG_lon_in, TG_lat_in, TG_alt_in;
string TG_name_in;

double GS_lon_in, GS_lat_in, GS_alt_in, GS_minelev_in;
string GS_name_in;

//int sensact_ind = 0;

// AreaType_pimpl
//

void AreaType_pimpl::
pre ()
{
}

void AreaType_pimpl::
post_AreaType ()
{
  const ::std::string& v (post_string ());

  // std::cout << "AreaType: " << v << std::endl;
}

// LengthType_pimpl
//

void LengthType_pimpl::
pre ()
{
}

void LengthType_pimpl::
post_LengthType ()
{
  const ::std::string& v (post_string ());

  // std::cout << "LengthType: " << v << std::endl;
}

// InertiaType_pimpl
//

void InertiaType_pimpl::
pre ()
{
}

void InertiaType_pimpl::
post_InertiaType ()
{
  const ::std::string& v (post_string ());

  // std::cout << "InertiaType: " << v << std::endl;
}

// MassType_pimpl
//

void MassType_pimpl::
pre ()
{
}

void MassType_pimpl::
post_MassType ()
{
  const ::std::string& v (post_string ());

  // std::cout << "MassType: " << v << std::endl;
}

// AngleType_pimpl
//

void AngleType_pimpl::
pre ()
{
}

void AngleType_pimpl::
post_AngleType ()
{
  const ::std::string& v (post_string ());

  // std::cout << "AngleType: " << v << std::endl;
}

// PositiveNumber_pimpl
//

void PositiveNumber_pimpl::
pre ()
{
}

void PositiveNumber_pimpl::
post_PositiveNumber ()
{
  double v (post_double ());
  
  positive_num_value = v;
  
  // std::cout << "PositiveNumber: " << v << std::endl;
}

void posV_pimpl::
pre ()
{
}

void posV_pimpl::
x (double x)
{
  // std::cout << "x: " << x << std::endl;
  
  // std::cout << "\nExecute x\n" << std::endl;
  x_value = x;
  // std::cout << x_value << std::endl;
  // std::cout << "\nExecute x\n" << std::endl;
}

void posV_pimpl::
y (double y)
{
  // std::cout << "y: " << y << std::endl;
  y_value = y;
}

void posV_pimpl::
z (double z)
{
  // std::cout << "z: " << z << std::endl;
  z_value = z;
}

void posV_pimpl::
name (const ::std::string& name)
{
  // std::cout << "name: " << name << std::endl;
}

void posV_pimpl::
unit (const ::std::string& unit)
{
}

void posV_pimpl::
post_posV ()
{
}

// velV_pimpl
//

void velV_pimpl::
pre ()
{
}

void velV_pimpl::
vx (double vx)
{
  vx_value = vx;
}

void velV_pimpl::
vy (double vy)
{
  vy_value = vy;
}

void velV_pimpl::
vz (double vz)
{
  vz_value = vz;
}

void velV_pimpl::
name (const ::std::string& name)
{
  // std::cout << "name: " << name << std::endl;
}

void velV_pimpl::
unit (const ::std::string& unit)
{
}

void velV_pimpl::
post_velV ()
{
}

// Vector_pimpl
//

void Vector_pimpl::
pre ()
{
}

void Vector_pimpl::
x (double x)
{
  //std::cout << "x: " << x << std::endl;
  x_value = x;
}

void Vector_pimpl::
y (double y)
{
  //std::cout << "y: " << y << std::endl;
  y_value = y;
}

void Vector_pimpl::
z (double z)
{
  //std::cout << "z: " << z << std::endl;
  z_value = z;
}

void Vector_pimpl::
name (const ::std::string& name)
{
  //std::cout << "name: " << name << std::endl;
}

void Vector_pimpl::
unit (const ::std::string& unit)
{
  //std::cout << "unit: " << unit << std::endl;
}

void Vector_pimpl::
post_Vector ()
{
}

// RotationMatrix_3x3_pimpl
//

void RotationMatrix_3x3_pimpl::
pre ()
{
}

void RotationMatrix_3x3_pimpl::
m11 (double m11)
{
  // std::cout << "m11: " << m11 << std::endl;
  m11_in = m11;
}

void RotationMatrix_3x3_pimpl::
m12 (double m12)
{
  // std::cout << "m12: " << m12 << std::endl;
  m12_in = m12;
}

void RotationMatrix_3x3_pimpl::
m13 (double m13)
{
  // std::cout << "m13: " << m13 << std::endl;
  m13_in = m13;
}

void RotationMatrix_3x3_pimpl::
m21 (double m21)
{
  // std::cout << "m21: " << m21 << std::endl;
  m21_in = m21;
}

void RotationMatrix_3x3_pimpl::
m22 (double m22)
{
  // std::cout << "m22: " << m22 << std::endl;
  m22_in = m22;
}

void RotationMatrix_3x3_pimpl::
m23 (double m23)
{
  // std::cout << "m23: " << m23 << std::endl;
  m23_in = m23;
}

void RotationMatrix_3x3_pimpl::
m31 (double m31)
{
  // std::cout << "m31: " << m31 << std::endl;
  m31_in = m31;
}

void RotationMatrix_3x3_pimpl::
m32 (double m32)
{
  // std::cout << "m32: " << m32 << std::endl;
  m32_in = m32;
}

void RotationMatrix_3x3_pimpl::
m33 (double m33)
{
  // std::cout << "m33: " << m33 << std::endl;
  m33_in = m33;
}

void RotationMatrix_3x3_pimpl::
name (const ::std::string& name)
{
  // std::cout << "name: " << name << std::endl;
}

void RotationMatrix_3x3_pimpl::
post_RotationMatrix_3x3 ()
{
}

// Dimensioned_pimpl
//

void Dimensioned_pimpl::
pre ()
{
}

void Dimensioned_pimpl::
unit (const ::std::string& unit)
{
  //std::cout << "unit: " << unit << std::endl;
}

void Dimensioned_pimpl::
post_Dimensioned ()
{
  post_PositiveNumber ();
}

// Altitude_pimpl
//

void Altitude_pimpl::
pre ()
{
}

void Altitude_pimpl::
unit ()
{
}

void Altitude_pimpl::
post_Altitude ()
{
  double v (post_double ());

  //std::cout << "Altitude: " << v << std::endl;
  alt_value = v;
}

// Angle_pimpl
//

void Angle_pimpl::
pre ()
{
}

void Angle_pimpl::
unit ()
{
}

void Angle_pimpl::
post_Angle ()
{
  double v (post_double ());
  
  angle_value = v;

  //std::cout << "Angle: " << v << std::endl;
}

// MoI_pimpl
//

void MoI_pimpl::
pre ()
{
}

void MoI_pimpl::
unit ()
{
}

void MoI_pimpl::
post_MoI ()
{
  post_PositiveNumber ();
}

// Mass_pimpl
//

void Mass_pimpl::
pre ()
{
}

void Mass_pimpl::
unit ()
{
}

void Mass_pimpl::
post_Mass ()
{
  post_PositiveNumber ();
}

// simparam_pimpl
//

void simparam_pimpl::
pre ()
{
}

void simparam_pimpl::
fileheader ()
{
}

void simparam_pimpl::
SC_Faces ()
{
}

void simparam_pimpl::
SC_properties ()
{
}

void simparam_pimpl::
InputFiles ()
{
}

void simparam_pimpl::
OutputFiles ()
{
}

void simparam_pimpl::
SimParameters ()
{
}

// eventsparam_pimpl
//

void eventsparam_pimpl::
pre ()
{
}

void eventsparam_pimpl::
fileheader ()
{
}

void eventsparam_pimpl::
CompParameters ()
{
}

void eventsparam_pimpl::
TGs ()
{
}

void eventsparam_pimpl::
GSs ()
{
}

void eventsparam_pimpl::
EventsInputFiles ()
{
}

void eventsparam_pimpl::
EventsOutputFiles ()
{
}

void eventsparam_pimpl::
name (const ::std::string& name)
{
  //std::cout << "Name: " << name << std::endl;
}

void eventsparam_pimpl::
post_eventsparam ()
{
}




void simparam_pimpl::
SensorsActuators ()
{
  SensAct_prms_in[sensact_ind].on_off = subsystem_on_in;
  
  SensAct_prms_in[sensact_ind].Name = name_in;
  
  SensAct_prms_in[sensact_ind].ConstPrm = constparam_in;
  
  SensAct_prms_in[sensact_ind].AuxPrm = auxparam_in;
  
  SensAct_prms_in[sensact_ind].OPS_limits = opslimits_in;
  
  SensAct_prms_in[sensact_ind].Accuracy = accuracy_in;
  
  SensAct_prms_in[sensact_ind].SC2SYS = SC2SYS_in;
  
  //cout << sensact_ind << endl;
  //cout << SensAct_prms_in[sensact_ind].Name << endl;
  //cout << SensAct_prms_in[sensact_ind].Accuracy << endl;
  //cout << SensAct_prms_in[sensact_ind].SC2SYS << endl;
  
  sensact_ind++;
  
  SC2SYS_in <<  0.0,  0.0, 0.0,
                0.0,  0.0,  0.0,
                0.0,  0.0,  0.0;
                
  string name_in = " ";
}

void simparam_pimpl::
Maneuvers ()
{
}

void simparam_pimpl::
name (const ::std::string& name)
{
  std::cout << "Name: " << name << std::endl;
}

void simparam_pimpl::
post_simparam ()
{
}

// fileheader_pimpl
//

void fileheader_pimpl::
pre ()
{
}

void fileheader_pimpl::
author (const ::std::string& author)
{
  std::cout << "Author: " << author << std::endl;
}

void fileheader_pimpl::
email (const ::std::string& email)
{
  std::cout << "email: " << email << std::endl;
}

void fileheader_pimpl::
organization (const ::std::string& organization)
{
  std::cout << "Organization: " << organization << std::endl;
}

void fileheader_pimpl::
license ()
{
}

void fileheader_pimpl::
sensitivity (const ::std::string& sensitivity)
{
  std::cout << "Sensitivity: " << sensitivity << std::endl;
}

void fileheader_pimpl::
filecreationdate (const ::xml_schema::date& filecreationdate)
{
  std::cout << "File creation date: "
   << filecreationdate.year () << '-'
   << filecreationdate.month () << '-'
   << filecreationdate.day ();

  if (filecreationdate.zone_present ())
  {
    if (filecreationdate.zone_hours () < 0)
      std::cout << filecreationdate.zone_hours () << ':' << -filecreationdate.zone_minutes ();
    else
      std::cout << '+' << filecreationdate.zone_hours () << ':' << filecreationdate.zone_minutes ();
  }

  std::cout << std::endl;
}

void fileheader_pimpl::
version (const ::std::string& version)
{
  std::cout << "Version: " << version << std::endl;
}

void fileheader_pimpl::
description (const ::std::string& description)
{
  std::cout << "Description: " << description << std::endl;
}

void fileheader_pimpl::
note (const ::std::string& note)
{
  std::cout << "Note: " << note << std::endl;
}

void fileheader_pimpl::
limitation (const ::std::string& limitation)
{
  std::cout << "Limitations: " << limitation << std::endl;
}

void fileheader_pimpl::
reference ()
{
}

void fileheader_pimpl::
post_fileheader ()
{
}

// reference_pimpl
//

void reference_pimpl::
pre ()
{
}

void reference_pimpl::
author (const ::std::string& author)
{
  //std::cout << "author: " << author << std::endl;
}

void reference_pimpl::
date (const ::std::string& date)
{
  //std::cout << "date: " << date << std::endl;
}

void reference_pimpl::
refID (const ::std::string& refID)
{
  //std::cout << "refID: " << refID << std::endl;
}

void reference_pimpl::
title (const ::std::string& title)
{
  //std::cout << "title: " << title << std::endl;
}

void reference_pimpl::
post_reference ()
{
}

// Versor_pimpl
//

void Versor_pimpl::
pre ()
{
}

void Versor_pimpl::
x (double x)
{
  // std::cout << "x: " << x << std::endl;
}

void Versor_pimpl::
y (double y)
{
  // std::cout << "y: " << y << std::endl;
}

void Versor_pimpl::
z (double z)
{
  // std::cout << "z: " << z << std::endl;
}

void Versor_pimpl::
name (const ::std::string& name)
{
  // std::cout << "name: " << name << std::endl;
}

void Versor_pimpl::
post_Versor ()
{
}

// nMAX_pimpl
//

void nMAX_pimpl::
pre ()
{
}

void nMAX_pimpl::
post_nMAX ()
{
  long long v (post_integer ());

  //std::cout << "nMAX: " << v << std::endl;
  
  nMAX_value = v;
}

// SC_Faces_pimpl
//

void SC_Faces_pimpl::
pre ()
{
}

void SC_Faces_pimpl::
Face ()
{
  SC_Face_in[face_ind].Area = face_area;
  SC_Face_in[face_ind].Material = face_material;
  SC_Face_in[face_ind].cP = face_cP;
  SC_Face_in[face_ind].cA = face_cA;
  
  face_ind++;
}

void SC_Faces_pimpl::
post_SC_Faces ()
{
}

// Face_pimpl
//

void Face_pimpl::
pre ()
{
}

void Face_pimpl::
Area ()
{
  face_area = positive_num_value;
}

void Face_pimpl::
Versor ()
{
}

void Face_pimpl::
Material (const ::std::string& Material)
{
  face_material = Material;
  // std::cout << "Material: " << Material << std::endl;
}

void Face_pimpl::
cP_position ()
{
  face_cP(0) = x_value;
  face_cP(1) = y_value;
  face_cP(2) = z_value;
}

void Face_pimpl::
cA_position ()
{
  face_cA(0) = x_value;
  face_cA(1) = y_value;
  face_cA(2) = z_value;
}

void Face_pimpl::
name (const ::std::string& name)
{
  // std::cout << "name: " << name << std::endl;
}

void Face_pimpl::
post_Face ()
{
}

// Length_pimpl
//

void Length_pimpl::
pre ()
{
}

void Length_pimpl::
unit ()
{
}

void Length_pimpl::
post_Length ()
{
  post_PositiveNumber ();
}

// Area_pimpl
//

void Area_pimpl::
pre ()
{
}

void Area_pimpl::
unit ()
{
}

void Area_pimpl::
post_Area ()
{
  post_PositiveNumber ();
}

void SC_properties_pimpl::
pre ()
{
}

void SC_properties_pimpl::
InertiaMatrix ()
{
}

void SC_properties_pimpl::
CoG ()
{
}

void SC_properties_pimpl::
Coefficients ()
{
}

void SC_properties_pimpl::
Areas ()
{
}

void SC_properties_pimpl::
SC_dipole ()
{
Mdip_in(0) = x_value;
Mdip_in(1) = y_value;
Mdip_in(2) = z_value;
}

void SC_properties_pimpl::
post_SC_properties ()
{
}

// InertiaMatrix_pimpl
//

void InertiaMatrix_pimpl::
pre ()
{
}

void InertiaMatrix_pimpl::
Ixx ()
{
  Ixx_in = positive_num_value;
}

void InertiaMatrix_pimpl::
Iyy ()
{
  Iyy_in = positive_num_value;
}

void InertiaMatrix_pimpl::
Izz ()
{
  Izz_in = positive_num_value;
}

void InertiaMatrix_pimpl::
Ixy (double Ixy)
{
  // std::cout << "Ixy: " << Ixy << std::endl;
  Ixy_in = Ixy;
}

void InertiaMatrix_pimpl::
Ixz (double Ixz)
{
  // std::cout << "Ixz: " << Ixz << std::endl;
  Ixz_in = Ixz;
}

void InertiaMatrix_pimpl::
Iyz (double Iyz)
{
  // std::cout << "Iyz: " << Iyz << std::endl;
  Iyz_in = Iyz;
}

void InertiaMatrix_pimpl::
unit ()
{
}

void InertiaMatrix_pimpl::
post_InertiaMatrix ()
{
}

// CoG_pimpl
//

void CoG_pimpl::
pre ()
{
}

void CoG_pimpl::
SC_mass()
{
  //SC_mass_in = *SC_mass_parser_->unit_parser_;
  //std::cout << "\nSC_mass: " << SC_mass_parser_->unit_parser(SC_mass_in) << "\n" << std::endl;
  
  //Mass_pskel* ptr1 = &SC_mass_parser_;
  
  //double** mass_value_ptr;
  //mass_value_ptr = &SC_mass_parser_;
  //double mass_value = **mass_value_ptr;
  
  SC_mass_in = positive_num_value;
  
  //std::cout << "\nExecute SC_mass_parser_\n" << std::endl;
  //std::cout << positive_num_value << std::endl;//SC_mass_parser_->post_PositiveNumber();//post_Mass();//
  //std::cout << "\nExecute SC_mass_parser_\n" << std::endl;
  
  //std::cout << "\nSC_mass: " << SC_mass_in << "\n" << std::endl;
}

void CoG_pimpl::
CoG_pos ()
{

//CoG_pos_parser_->x(x_value);
CoG_pos_vec(0) = x_value;
//CoG_pos_parser_->y(y_value);
CoG_pos_vec(1) = y_value;
//CoG_pos_parser_->z(z_value);
CoG_pos_vec(2) = z_value;

//std::cout << "\nExecute CoG_pos\n" << std::endl;

//std::cout << "\nExecute CoG_pos\n" << std::endl;
//CoG_pos_parser_->x(x_value);//post_Mass();//
//std::cout << x_value << std::endl;
//std::cout << "\nExecute CoG_pos\n" << std::endl;  
 
}

void CoG_pimpl::
post_CoG ()
{
}

// Coefficients_pimpl
//

void Coefficients_pimpl::
pre ()
{
}

void Coefficients_pimpl::
Cd ()
{
  Cd_in = positive_num_value;
}

void Coefficients_pimpl::
Cr ()
{
  Cr_in = positive_num_value;
}

void Coefficients_pimpl::
post_Coefficients ()
{
}

// Areas_pimpl
//

void Areas_pimpl::
pre ()
{
}

void Areas_pimpl::
Area_D ()
{
  Area_D_in = positive_num_value;
}

void Areas_pimpl::
Area_R ()
{
  Area_R_in = positive_num_value;
}

void Areas_pimpl::
post_Areas ()
{
}

// SimParameters_pimpl
//

void SimParameters_pimpl::
pre ()
{
}

void SimParameters_pimpl::
durstep ()
{
}

void SimParameters_pimpl::
ORB_initstate ()
{
}

void SimParameters_pimpl::
ATT_initstate ()
{
}

void SimParameters_pimpl::
simoptions ()
{
}

void SimParameters_pimpl::
post_SimParameters ()
{
}

// durstep_pimpl
//

void durstep_pimpl::
pre ()
{
}

void durstep_pimpl::
simstep (const ::xml_schema::duration& simstep)
{
  if( simstep.years() != 0 || simstep.months() != 0 || simstep.days() != 0 || simstep.hours() != 0 )
    {
    cerr << "Years, months, days and hours are not allowed in the definition of the simulation step" << endl;
    exit(EXIT_FAILURE);  
    }
  
  sim_step = 60*simstep.minutes() + simstep.seconds(); // sim_step is always treated as number of seconds by the simulator software
}

void durstep_pimpl::
simduration (const ::xml_schema::duration& simduration)
{
  if( simduration.years() != 0 || simduration.months() != 0)
    {
    cerr << "Years and months are not allowed in the definition of the simulation duration" << endl;
    exit(EXIT_FAILURE);  
    }
  
  sim_duration = 86400*simduration.days() + 3600*simduration.hours() + 60*simduration.minutes() + simduration.seconds(); // sim_step is always treated as number of seconds by the simulator software
}

void durstep_pimpl::
post_durstep ()
{
}


// ORB_initstate_pimpl
//

void ORB_initstate_pimpl::
pre ()
{
}

void ORB_initstate_pimpl::
Initime (const ::xml_schema::date_time& Initime)
{
  //std::cout << "Initime: "
  UTCdate(0) = Initime.year();
  UTCdate(1) = Initime.month();
  UTCdate(2) = Initime.day();
  UTCdate(3) = Initime.hours();
  UTCdate(4) = Initime.minutes();
  UTCdate(5) = Initime.seconds();

  //if (Initime.zone_present ())
  //{
  //  if (Initime.zone_hours () < 0)
  //    std::cout << Initime.zone_hours () << ':' << -Initime.zone_minutes ();
  //  else
  //    std::cout << '+' << Initime.zone_hours () << ':' << Initime.zone_minutes ();
  //}
  //
  //std::cout << std::endl;
}

void ORB_initstate_pimpl::
Position ()
{
  Pos_vec(0) = x_value;
  Pos_vec(1) = y_value;
  Pos_vec(2) = z_value;
}

void ORB_initstate_pimpl::
Velocity ()
{
  Vel_vec(0) = vx_value;
  Vel_vec(1) = vy_value;
  Vel_vec(2) = vz_value;
}

void ORB_initstate_pimpl::
post_ORB_initstate ()
{
}

// ATT_initstate_pimpl
//

void ATT_initstate_pimpl::
pre ()
{
}

void ATT_initstate_pimpl::
phi ()
{
  phi_in = angle_value;
}

void ATT_initstate_pimpl::
theta ()
{
  theta_in = angle_value;
}

void ATT_initstate_pimpl::
psi ()
{
  psi_in = angle_value;
}

void ATT_initstate_pimpl::
om_x ()
{
  om_x_in = angle_value;
}

void ATT_initstate_pimpl::
om_y ()
{
  om_y_in = angle_value;
}

void ATT_initstate_pimpl::
om_z ()
{
  om_z_in = angle_value;
}

void ATT_initstate_pimpl::
post_ATT_initstate ()
{
}

// simoptions_pimpl
//

void simoptions_pimpl::
pre ()
{
}

void simoptions_pimpl::
initstate_in_RTN (bool initstate_in_RTN)
{
  //std::cout << "initstate_in_RTN: " << initstate_in_RTN << std::endl;
  initstate_in_RTN_in = initstate_in_RTN;
}

void simoptions_pimpl::
realtime (bool realtime)
{
  //std::cout << "realtime: " << realtime << std::endl;
  realtime_in = realtime;
}

void simoptions_pimpl::
realtime_wait ()
{
  realtime_wait_in = positive_num_value;
}

void simoptions_pimpl::
ggrad_on (bool ggrad_on)
{
  //std::cout << "ggrad_on: " << ggrad_on << std::endl;
  ggrad_on_in = ggrad_on;
}

void simoptions_pimpl::
mag_on (bool mag_on)
{
  //std::cout << "mag_on: " << mag_on << std::endl;
  mag_on_in = mag_on;
}

void simoptions_pimpl::
srp_on (bool srp_on)
{
  //std::cout << "srp_on: " << srp_on << std::endl;
  srp_on_in = srp_on;
}

void simoptions_pimpl::
drag_on (bool drag_on)
{
  // std::cout << "drag_on: " << drag_on << std::endl;
  drag_on_in = drag_on;
}

void simoptions_pimpl::
nMAX ()
{
  nMAX_in = nMAX_value;
}

void simoptions_pimpl::
sunmoon_on (bool sunmoon_on)
{
  //std::cout << "sunmoon_on: " << sunmoon_on << std::endl;
  sunmoon_on_in = sunmoon_on;
}

void simoptions_pimpl::
Drag_Model (const ::std::string& Drag_Model)
{
  //std::cout << "Drag_Model: " << Drag_Model << std::endl;
  Drag_Model_in = Drag_Model;
}

void simoptions_pimpl::
SRP_Model (const ::std::string& SRP_Model)
{
  //std::cout << "SRP_Model: " << SRP_Model << std::endl;
  SRP_Model_in = SRP_Model;
}

void simoptions_pimpl::
AttitudeType (const ::std::string& AttitudeType)
{
  //std::cout << "AttitudeType: " << AttitudeType << std::endl;
  AttitudeType_in = AttitudeType;
}

void simoptions_pimpl::
attctrl_on (bool attctrl_on)
{
  //std::cout << "attctrl_on: " << attctrl_on << std::endl;
  attctrl_on_in = attctrl_on;
}

void simoptions_pimpl::
AttCtrlType (const ::std::string& AttCtrlType)
{
  //std::cout << "AttCtrlType: " << AttCtrlType << std::endl;
  AttCtrlType_in = AttCtrlType;
}

void simoptions_pimpl::
orbctrl_on (bool orbctrl_on)
{
  //std::cout << "orbctrl_on: " << orbctrl_on << std::endl;
  orbctrl_on_in = orbctrl_on;
}

void simoptions_pimpl::
OrbCtrlType (const ::std::string& OrbCtrlType)
{
  //std::cout << "OrbCtrlType: " << OrbCtrlType << std::endl;
  OrbCtrlType_in = OrbCtrlType;
}

void simoptions_pimpl::
post_simoptions ()
{
}

// InputFiles_pimpl
//

void InputFiles_pimpl::
pre ()
  {
  }

void InputFiles_pimpl::
Orbit_ephemeris (const ::std::string& Orbit_ephemeris)
  {
  //std::cout << "Orbit_ephemeris: " << Orbit_ephemeris << std::endl;
  Orbit_ephemeris_in = Orbit_ephemeris;
  }
  
void InputFiles_pimpl::
Attitude_ephemeris (const ::std::string& Attitude_ephemeris)
{
  //std::cout << "Attitude_ephemeris: " << Attitude_ephemeris << std::endl;
  Attitude_ephemeris_in = Attitude_ephemeris;
}

void InputFiles_pimpl::
Data_path (const ::std::string& Data_path)
{
  //std::cout << "Data_path: " << Data_path << std::endl;
  Data_path_in = Data_path;
}

void InputFiles_pimpl::
Planet_ephemeris (const ::std::string& Planet_ephemeris)
  {
  //std::cout << "Planet_ephemeris: " << Planet_ephemeris << std::endl;
  planetephemeris_in = Planet_ephemeris;
  }

void InputFiles_pimpl::
EOP_parameters (const ::std::string& EOP_parameters)
  {
  //std::cout << "EOP_parameters: " << EOP_parameters << std::endl;
  eop_in = EOP_parameters;
  }

void InputFiles_pimpl::
PCK_data (const ::std::string& PCK_data)
  {
  //std::cout << "PCK_data: " << PCK_data << std::endl;
  pck_data_in = PCK_data;
  }

void InputFiles_pimpl::
Leap_second (const ::std::string& Leap_second)
  {
  //std::cout << "Leap_second: " << Leap_second << std::endl;
  leapsecond_in = Leap_second;
  }

void InputFiles_pimpl::
Gravity_model (const ::std::string& Gravity_model)
  {
  //std::cout << "Gravity_model: " << Gravity_model << std::endl;
  gravityfield_in = Gravity_model;
  }

void InputFiles_pimpl::
Atmospheric_model (const ::std::string& Atmospheric_model)
  {
  //std::cout << "Atmospheric_model: " << Atmospheric_model << std::endl;
  atmosphere_in = Atmospheric_model;
  }

void InputFiles_pimpl::
Magnetic_model (const ::std::string& Magnetic_model)
  {
  //std::cout << "Magnetic_model: " << Magnetic_model << std::endl;
  magn_model_in = Magnetic_model;
  }

void InputFiles_pimpl::
name (const ::std::string& name)
  {
  // std::cout << "name: " << name << std::endl;
  }

void InputFiles_pimpl::
post_InputFiles ()
  {
  }

// OutputFiles_pimpl
//

void OutputFiles_pimpl::
pre ()
{
}

void OutputFiles_pimpl::
Orbit_ephemeris (const ::std::string& Orbit_ephemeris)
{
  //std::cout << "Orbit_ephemeris: " << Orbit_ephemeris << std::endl;
  orbfile_name_in = Orbit_ephemeris;
}

void OutputFiles_pimpl::
Attitude_ephemeris (const ::std::string& Attitude_ephemeris)
  {
  //std::cout << "Attitude_ephemeris: " << Attitude_ephemeris << std::endl;
  attfile_name_in = Attitude_ephemeris;
  }

void OutputFiles_pimpl::
Sensor_output (const ::std::string& Sensor_output)
  {
  //std::cout << "Sensor_output: " << Sensor_output << std::endl;
  sensors_filename_in = Sensor_output;
  }

void OutputFiles_pimpl::
Torques (const ::std::string& Torques)
  {
  //std::cout << "Torques: " << Torques << std::endl;
  csv_torques_name_in = Torques;
  }
  
void OutputFiles_pimpl::
Accelerations (const ::std::string& Accelerations)
{
  //std::cout << "Accelerations: " << Accelerations << std::endl;
  csv_accelerations_name_in = Accelerations;
}

void OutputFiles_pimpl::
name (const ::std::string& name)
  {
  // std::cout << "name: " << name << std::endl;
  }

void OutputFiles_pimpl::
post_OutputFiles ()
{
}

// SensorsActuators_pimpl
//

void SensorsActuators_pimpl::
pre ()
{
  //sensact_ind = 0;
  constparamind = 0;
  auxparamind = 0;
  opslimitsind = 0;
  accuracyind = 0;
  
  //for(int i = 0; i < 10; i++)
  //  {
  //  constparam_in(i) = 0.0;
  //  auxparam_in(i) = 0.0;
  //  opslimits_in(i) = 0.0;
  //  accuracy_in(i) = 0.0;
  //  }
  
}

void SensorsActuators_pimpl::
subsystem_on (bool subsystem_on)
{
  // std::cout << "subsystem_on: " << subsystem_on << std::endl;
  subsystem_on_in = subsystem_on;
}


void SensorsActuators_pimpl::
constparam ()
{
  constparam_in.conservativeResize(constparamind+1);
  constparam_in(constparamind) = constpar;
  //std::cout << "constpar: " << constpar << std::endl;
  constparamind++;
}

void SensorsActuators_pimpl::
auxparam ()
{
  auxparam_in.conservativeResize(auxparamind+1);
  auxparam_in(auxparamind) = auxpar;
  //std::cout << "auxpar: " << auxpar << std::endl;
  auxparamind++;
}

void SensorsActuators_pimpl::
opslimit ()
{
  opslimits_in.conservativeResize(opslimitsind+1);
  opslimits_in(opslimitsind) = opslim;
  //std::cout << "opslim: " << opslim << std::endl;
  opslimitsind++;
}

void SensorsActuators_pimpl::
accuracy ()
{
  accuracy_in.conservativeResize(accuracyind+1);
  accuracy_in(accuracyind) = accur;
  //std::cout << "accur: " << accur << std::endl;
  accuracyind++;
}



void SensorsActuators_pimpl::
SC2SYS_matrix ()
{
 SC2SYS_in << m11_in, m12_in, m13_in,
              m21_in, m22_in, m23_in,
              m31_in, m32_in, m33_in;
              
  //cout << SC2SYS_in << endl;
  
}


void SensorsActuators_pimpl::
name (const ::std::string& name)
{
  // std::cout << "name: " << name << std::endl;
  name_in = name;
  
}

void SensorsActuators_pimpl::
post_SensorsActuators ()
{
  
  
}

// constparam_pimpl
//

void constparam_pimpl::
pre ()
{
}

void constparam_pimpl::
name (const ::std::string& name)
{
  //std::cout << "name: " << name << std::endl;
}

void constparam_pimpl::
unit (const ::std::string& unit)
{
  //std::cout << "unit: " << unit << std::endl;
}

void constparam_pimpl::
post_constparam ()
{
  double v (post_double ());

  //std::cout << "constparam: " << v << std::endl;
  constpar = v;
}

// auxparam_pimpl
//

void auxparam_pimpl::
pre ()
{
}

void auxparam_pimpl::
name (const ::std::string& name)
{
  //std::cout << "name: " << name << std::endl;
}

void auxparam_pimpl::
unit (const ::std::string& unit)
{
  //std::cout << "unit: " << unit << std::endl;
}

void auxparam_pimpl::
post_auxparam ()
{
  double v (post_double ());

  //std::cout << "auxparam: " << v << std::endl;
  auxpar = v;
}


// opslimit_pimpl
//

void opslimit_pimpl::
pre ()
{
}

void opslimit_pimpl::
name (const ::std::string& name)
{
  //std::cout << "name: " << name << std::endl;
}

void opslimit_pimpl::
unit (const ::std::string& unit)
{
  //std::cout << "unit: " << unit << std::endl;
}

void opslimit_pimpl::
post_opslimit ()
{
  double v (post_double ());
  
  opslim = v;
  //std::cout << "opslimit: " << v << std::endl;
}

// accuracy_pimpl
//

void accuracy_pimpl::
pre ()
{
}

void accuracy_pimpl::
name (const ::std::string& name)
{
  //std::cout << "name: " << name << std::endl;
}

void accuracy_pimpl::
unit (const ::std::string& unit)
{
  //std::cout << "unit: " << unit << std::endl;
}

void accuracy_pimpl::
post_accuracy ()
{
  double v (post_double ());
  
  accur = v;
  //std::cout << "accuracy: " << v << std::endl;
}

// Maneuvers_pimpl
//

void Maneuvers_pimpl::
pre ()
    {
    man_name_in = " ";;
    maneuver_on_in = false;
    man_init_time_in = 0.0;
    man_duration_in = 0.0;
    ManVec_in = Vec3d::Zero();
    }

void Maneuvers_pimpl::
Man ()
    {
    maneuver_in.name = man_name_in;
    maneuver_in.maneuver_on = maneuver_on_in;
    maneuver_in.init_time = man_init_time_in;
    maneuver_in.duration = man_duration_in;
    maneuver_in.ManVec = ManVec_in;
      
    all_maneuvers.push_back(maneuver_in);
    
    //cout << all_maneuvers.size() << endl;
    
    man_name_in = " ";;
    maneuver_on_in = false;
    man_init_time_in = 0.0;
    man_duration_in = 0.0;
    ManVec_in = Vec3d::Zero(); 
    
    
    }

void Maneuvers_pimpl::
post_Maneuvers ()
{
}

// Man_pimpl
//

void Man_pimpl::
pre ()
{
}

void Man_pimpl::
maneuver_on (bool maneuver_on)
{
  //std::cout << "maneuver_on: " << maneuver_on << std::endl;
  maneuver_on_in = maneuver_on;
}

void Man_pimpl::
init_time (double init_time)
{
  //std::cout << "init_time: " << init_time << std::endl;
  man_init_time_in = init_time;
}

void Man_pimpl::
duration (double duration)
{
  //std::cout << "duration: " << duration << std::endl;
  man_duration_in = duration;
}

void Man_pimpl::
ManVec ()
{
  ManVec_in(0) = x_value;
  ManVec_in(1) = y_value;
  ManVec_in(2) = z_value;
}

void Man_pimpl::
name ()
{
}

void Man_pimpl::
post_Man ()
{
}

// CompParameters_pimpl
//

void CompParameters_pimpl::
pre ()
{
}

void CompParameters_pimpl::
durstep ()
{
}

void CompParameters_pimpl::
Payload ()
{
}

void CompParameters_pimpl::
Spacecraft ()
{
}

void CompParameters_pimpl::
Compoptions ()
{
}

void CompParameters_pimpl::
post_CompParameters ()
{
}

// Payload_pimpl
//

void Payload_pimpl::
pre ()
{
}

void Payload_pimpl::
FOV_cross ()
{
  FOV_cross_in = angle_value;
}

void Payload_pimpl::
FOV_along ()
{
  FOV_along_in = angle_value;
}

void Payload_pimpl::
post_Payload ()
{
}

// Spacecraft_pimpl
//

void Spacecraft_pimpl::
pre ()
{
}

void Spacecraft_pimpl::
SC_start (unsigned long long SC_start)
{
  //std::cout << "SC_start: " << SC_start << std::endl;
  SC_start_in = SC_start;
}

void Spacecraft_pimpl::
SC_end (unsigned long long SC_end)
{
  //std::cout << "SC_end: " << SC_end << std::endl;
  SC_end_in = SC_end;
}

void Spacecraft_pimpl::
PL_start (unsigned long long PL_start)
{
  //std::cout << "PL_start: " << PL_start << std::endl;
  PL_start_in = PL_start;
}

void Spacecraft_pimpl::
PL_end (unsigned long long PL_end)
{
  //std::cout << "PL_end: " << PL_end << std::endl;
  PL_end_in = PL_end;
}

void Spacecraft_pimpl::
post_Spacecraft ()
{
}

// Compoptions_pimpl
//

void Compoptions_pimpl::
pre ()
{
}

void Compoptions_pimpl::
TGs_on (bool TGs_on)
{
  //std::cout << "TGs_on: " << TGs_on << std::endl;
  TGs_on_in = TGs_on;
}

void Compoptions_pimpl::
GSs_on (bool GSs_on)
{
  //std::cout << "GSs_on: " << GSs_on << std::endl;
  GSs_on_in = GSs_on;
}

void Compoptions_pimpl::
TGs_grid_on (bool TGs_grid_on)
{
  //std::cout << "TGs_grid_on: " << TGs_grid_on << std::endl;
  TGs_grid_on_in = TGs_grid_on;
}

void Compoptions_pimpl::
Eclipse_on (bool Eclipse_on)
{
  //std::cout << "Eclipse_on: " << Eclipse_on << std::endl;
  Eclipse_on_in = Eclipse_on;
}

void Compoptions_pimpl::
post_Compoptions ()
{
}

// TGs_pimpl
//

void TGs_pimpl::
pre ()
{
}

void TGs_pimpl::
TGs_grid ()
{
}

void TGs_pimpl::
TGs_list ()
{
}

void TGs_pimpl::
post_TGs ()
{
}

// TGs_grid_pimpl
//

void TGs_grid_pimpl::
pre ()
{
}

void TGs_grid_pimpl::
minlon ()
{
  minlon_in = angle_value;
}

void TGs_grid_pimpl::
maxlon ()
{
  maxlon_in = angle_value;
}

void TGs_grid_pimpl::
minlat ()
{
  minlat_in = angle_value;
}

void TGs_grid_pimpl::
maxlat ()
{
  maxlat_in = angle_value;
}

void TGs_grid_pimpl::
gridstep ()
{
  gridstep_in = angle_value;
}

void TGs_grid_pimpl::
post_TGs_grid ()
{
}

// TGs_list_pimpl
//

void TGs_list_pimpl::
pre ()
{
  TG_ind = 0;
}

void TGs_list_pimpl::
TG ()
{
  TGs_list_in[TG_ind].name = TG_name_in;
  TGs_list_in[TG_ind].lon = TG_lon_in;
  TGs_list_in[TG_ind].lat = TG_lat_in;
  TGs_list_in[TG_ind].alt = TG_alt_in;
  
  TG_ind++;
}

void TGs_list_pimpl::
post_TGs_list ()
{
}

// TG_pimpl
//

void TG_pimpl::
pre ()
{
}

void TG_pimpl::
lon ()
{
  TG_lon_in = angle_value;
}

void TG_pimpl::
lat ()
{
  TG_lat_in = angle_value;
}

void TG_pimpl::
alt ()
{
  TG_alt_in = alt_value;
}

void TG_pimpl::
name (const ::std::string& name)
{
  //std::cout << "name: " << name << std::endl;
  TG_name_in = name;
}

void TG_pimpl::
post_TG ()
{
}

// GSs_pimpl
//

void GSs_pimpl::
pre ()
{
  GS_ind = 0;
}

void GSs_pimpl::
GS ()
{
  GSs_list_in[GS_ind].name = GS_name_in;
  GSs_list_in[GS_ind].lon = GS_lon_in;
  GSs_list_in[GS_ind].lat = GS_lat_in;
  GSs_list_in[GS_ind].alt = GS_alt_in;
  GSs_list_in[GS_ind].minelev = GS_minelev_in;
  
  GS_ind++;
}

void GSs_pimpl::
post_GSs ()
{
}

// GS_pimpl
//

void GS_pimpl::
pre ()
{
}

void GS_pimpl::
lon ()
{
  GS_lon_in = angle_value;
}

void GS_pimpl::
lat ()
{
  GS_lat_in = angle_value;
}

void GS_pimpl::
alt ()
{
  GS_alt_in = alt_value;
}

void GS_pimpl::
minelev ()
{
  GS_minelev_in = angle_value;
}

void GS_pimpl::
name (const ::std::string& name)
{
  //std::cout << "name: " << name << std::endl;
  GS_name_in = name;
}

void GS_pimpl::
post_GS ()
{
}

// EventsInputFiles_pimpl
//

void EventsInputFiles_pimpl::
pre ()
{
}

void EventsInputFiles_pimpl::
Orbit_ephemeris_path (const ::std::string& Orbit_ephemeris_path)
{
  //std::cout << "Orbit_ephemeris_path: " << Orbit_ephemeris_path << std::endl;
  Orbit_ephemeris_path_in = Orbit_ephemeris_path;
}

void EventsInputFiles_pimpl::
Orbit_ephemeris_rootname (const ::std::string& Orbit_ephemeris_rootname)
{
  //std::cout << "Orbit_ephemeris_rootname: " << Orbit_ephemeris_rootname << std::endl;
  Orbit_ephemeris_rootname_in = Orbit_ephemeris_rootname;
}

void EventsInputFiles_pimpl::
Data_path (const ::std::string& Data_path)
{
  //std::cout << "Data_path: " << Data_path << std::endl;
  Data_path_in = Data_path;
}

void EventsInputFiles_pimpl::
Planet_ephemeris (const ::std::string& Planet_ephemeris)
{
  //std::cout << "Planet_ephemeris: " << Planet_ephemeris << std::endl;
  planetephemeris_in = Planet_ephemeris;
}

void EventsInputFiles_pimpl::
EOP_parameters (const ::std::string& EOP_parameters)
{
  //std::cout << "EOP_parameters: " << EOP_parameters << std::endl;
  eop_in = EOP_parameters;
}

void EventsInputFiles_pimpl::
PCK_data (const ::std::string& PCK_data)
{
  //std::cout << "PCK_data: " << PCK_data << std::endl;
  pck_data_in = PCK_data;
}

void EventsInputFiles_pimpl::
Leap_second (const ::std::string& Leap_second)
{
  //std::cout << "Leap_second: " << Leap_second << std::endl;
  leapsecond_in = Leap_second;
}

void EventsInputFiles_pimpl::
name (const ::std::string& name)
{
  //std::cout << "name: " << name << std::endl;
}

void EventsInputFiles_pimpl::
post_EventsInputFiles ()
{
}

// EventsOutputFiles_pimpl
//

void EventsOutputFiles_pimpl::
pre ()
{
}

void EventsOutputFiles_pimpl::
TG_contacts (const ::std::string& TG_contacts)
{
  //std::cout << "TG_contacts: " << TG_contacts << std::endl;
  TG_contacts_in = TG_contacts;
}

void EventsOutputFiles_pimpl::
GS_contacts (const ::std::string& GS_contacts)
{
  //std::cout << "GS_contacts: " << GS_contacts << std::endl;
  GS_contacts_in = GS_contacts;
}

void EventsOutputFiles_pimpl::
Eclipse_times (const ::std::string& Eclipse_times)
{
  //std::cout << "Eclipse_times: " << Eclipse_times << std::endl;
  Eclipse_times_in = Eclipse_times;
}

void EventsOutputFiles_pimpl::
name (const ::std::string& name)
{
  //std::cout << "name: " << name << std::endl;
}

void EventsOutputFiles_pimpl::
post_EventsOutputFiles ()
{
}

// license_pimpl
//

void license_pimpl::
pre ()
{
}

void license_pimpl::
licenseName (const ::std::string& licenseName)
{
  std::cout << "License name: " << licenseName << std::endl;
}

void license_pimpl::
licenseURL (const ::std::string& licenseURL)
{
  std::cout << "License URL: " << licenseURL << std::endl;
}

void license_pimpl::
post_license ()
{
}

// name_pimpl
//

void name_pimpl::
pre ()
{
}

void name_pimpl::
post_name ()
{
  const ::std::string& v (post_string ());

  //std::cout << "name: " << v << std::endl;
  
  man_name_in = v;
}





