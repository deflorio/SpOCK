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

#ifndef SIMPARAM_SCHEMA_PIMPL_HXX
#define SIMPARAM_SCHEMA_PIMPL_HXX

#include <simparam_schema-pskel.h>

#include <Eigen/Core>
#include <VarTypes.h>

using namespace std;

class AreaType_pimpl: public virtual AreaType_pskel,
  public ::xml_schema::string_pimpl
{
  public:
  virtual void
  pre ();

  virtual void
  post_AreaType ();
};

class LengthType_pimpl: public virtual LengthType_pskel,
  public ::xml_schema::string_pimpl
{
  public:
  virtual void
  pre ();

  virtual void
  post_LengthType ();
};

class InertiaType_pimpl: public virtual InertiaType_pskel,
  public ::xml_schema::string_pimpl
{
  public:
  virtual void
  pre ();

  virtual void
  post_InertiaType ();
};

class MassType_pimpl: public virtual MassType_pskel,
  public ::xml_schema::string_pimpl
{
  public:
  virtual void
  pre ();

  virtual void
  post_MassType ();
};

class AngleType_pimpl: public virtual AngleType_pskel,
  public ::xml_schema::string_pimpl
{
  public:
  virtual void
  pre ();

  virtual void
  post_AngleType ();
};

class PositiveNumber_pimpl: public virtual PositiveNumber_pskel,
  public ::xml_schema::double_pimpl
{
  public:
  virtual void
  pre ();

  virtual void
  post_PositiveNumber ();
};

class posV_pimpl: public virtual posV_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  x (double);

  virtual void
  y (double);

  virtual void
  z (double);

  virtual void
  name (const ::std::string&);

  virtual void
  unit (const ::std::string&);

  virtual void
  post_posV ();
};

class velV_pimpl: public virtual velV_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  vx (double);

  virtual void
  vy (double);

  virtual void
  vz (double);

  virtual void
  name (const ::std::string&);

  virtual void
  unit (const ::std::string&);

  virtual void
  post_velV ();
};

class Vector_pimpl: public virtual Vector_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  x (double);

  virtual void
  y (double);

  virtual void
  z (double);

  virtual void
  name (const ::std::string&);
  
  virtual void
  unit (const ::std::string&);

  virtual void
  post_Vector ();
};

class RotationMatrix_3x3_pimpl: public virtual RotationMatrix_3x3_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  m11 (double);

  virtual void
  m12 (double);

  virtual void
  m13 (double);

  virtual void
  m21 (double);

  virtual void
  m22 (double);

  virtual void
  m23 (double);

  virtual void
  m31 (double);

  virtual void
  m32 (double);

  virtual void
  m33 (double);

  virtual void
  name (const ::std::string&);

  virtual void
  post_RotationMatrix_3x3 ();
};

class Dimensioned_pimpl: public virtual Dimensioned_pskel,
  public ::PositiveNumber_pimpl
{
  public:
  virtual void
  pre ();

  virtual void
  unit (const ::std::string&);

  virtual void
  post_Dimensioned ();
};

class Altitude_pimpl: public virtual Altitude_pskel,
  public ::xml_schema::double_pimpl
{
  public:
  virtual void
  pre ();

  virtual void
  unit ();

  virtual void
  post_Altitude ();
};

class Angle_pimpl: public virtual Angle_pskel,
  public ::xml_schema::double_pimpl
{
  public:
  virtual void
  pre ();

  virtual void
  unit ();

  virtual void
  post_Angle ();
};

class MoI_pimpl: public virtual MoI_pskel,
  public ::PositiveNumber_pimpl
{
  public:
  virtual void
  pre ();

  virtual void
  unit ();

  virtual void
  post_MoI ();
};

class Mass_pimpl: public virtual Mass_pskel,
  public ::PositiveNumber_pimpl
{
  public:
  virtual void
  pre ();

  virtual void
  unit ();

  virtual void
  post_Mass ();
};

class simparam_pimpl: public virtual simparam_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  fileheader ();

  virtual void
  SC_Faces ();

  virtual void
  SC_properties ();

  virtual void
  InputFiles ();

  virtual void
  OutputFiles ();

  virtual void
  SimParameters ();

  virtual void
  SensorsActuators ();
  
  virtual void
  Maneuvers ();

  virtual void
  name (const ::std::string&);

  virtual void
  post_simparam ();
  
  public:
  
  //VectorNd<2> constprm, opslimits;
  //Vec4d accuracy;
  SC::SYS_params SensAct_prms_in[30];
  static int sensact_ind;
  
};

class eventsparam_pimpl: public virtual eventsparam_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  fileheader ();

  virtual void
  CompParameters ();

  virtual void
  TGs ();

  virtual void
  GSs ();

  virtual void
  EventsInputFiles ();

  virtual void
  EventsOutputFiles ();

  virtual void
  name (const ::std::string&);

  virtual void
  post_eventsparam ();
};

class fileheader_pimpl: public virtual fileheader_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  author (const ::std::string&);

  virtual void
  email (const ::std::string&);

  virtual void
  organization (const ::std::string&);

  virtual void
  license ();

  virtual void
  sensitivity (const ::std::string&);

  virtual void
  filecreationdate (const ::xml_schema::date&);

  virtual void
  version (const ::std::string&);

  virtual void
  description (const ::std::string&);

  virtual void
  note (const ::std::string&);

  virtual void
  limitation (const ::std::string&);

  virtual void
  reference ();

  virtual void
  post_fileheader ();
};

class reference_pimpl: public virtual reference_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  author (const ::std::string&);

  virtual void
  date (const ::std::string&);

  virtual void
  refID (const ::std::string&);

  virtual void
  title (const ::std::string&);

  virtual void
  post_reference ();
};

class Versor_pimpl: public virtual Versor_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  x (double);

  virtual void
  y (double);

  virtual void
  z (double);

  virtual void
  name (const ::std::string&);

  virtual void
  post_Versor ();
};

class nMAX_pimpl: public virtual nMAX_pskel,
  public ::xml_schema::integer_pimpl
{
  public:
  virtual void
  pre ();

  virtual void
  post_nMAX ();
};

class SC_Faces_pimpl: public virtual SC_Faces_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  Face ();

  virtual void
  post_SC_Faces ();
  
  public:
  
  SC::Face SC_Face_in[6];
  int face_ind;
};

class Face_pimpl: public virtual Face_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  Area ();

  virtual void
  Versor ();

  virtual void
  Material (const ::std::string&);

  virtual void
  cP_position ();

  virtual void
  cA_position ();

  virtual void
  name (const ::std::string&);

  virtual void
  post_Face ();
  
};

class Length_pimpl: public virtual Length_pskel,
  public ::PositiveNumber_pimpl
{
  public:
  virtual void
  pre ();

  virtual void
  unit ();

  virtual void
  post_Length ();
};

class Area_pimpl: public virtual Area_pskel,
  public ::PositiveNumber_pimpl
{
  public:
  virtual void
  pre ();

  virtual void
  unit ();

  virtual void
  post_Area ();
};

class SC_properties_pimpl: public virtual SC_properties_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  InertiaMatrix ();

  virtual void
  CoG ();

  virtual void
  Coefficients ();
  
  virtual void
  Areas ();

  virtual void
  SC_dipole ();

  virtual void
  post_SC_properties ();
  
  public:
  
  Vec3d Mdip_in;
};

class InertiaMatrix_pimpl: public virtual InertiaMatrix_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  Ixx ();

  virtual void
  Iyy ();

  virtual void
  Izz ();

  virtual void
  Ixy (double);

  virtual void
  Ixz (double);

  virtual void
  Iyz (double);
  
  virtual void
  unit ();

  virtual void
  post_InertiaMatrix ();
  
  public:
  
  double Ixx_in, Iyy_in, Izz_in, Ixy_in, Ixz_in, Iyz_in;
};

class CoG_pimpl: public virtual CoG_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  SC_mass ();

  virtual void
  CoG_pos ();

  virtual void
  post_CoG ();
  
  public:
    
  double SC_mass_in;
  Vec3d CoG_pos_vec;
  
};

class Coefficients_pimpl: public virtual Coefficients_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  Cd ();

  virtual void
  Cr ();

  virtual void
  post_Coefficients ();
  
  public:
    
  double Cd_in, Cr_in;
};

class Areas_pimpl: public virtual Areas_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  Area_D ();

  virtual void
  Area_R ();

  virtual void
  post_Areas ();
  
  double Area_D_in, Area_R_in;
};

class SimParameters_pimpl: public virtual SimParameters_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  durstep ();
  
  virtual void
  ORB_initstate ();

  virtual void
  ATT_initstate ();

  virtual void
  simoptions ();

  virtual void
  post_SimParameters ();
};

class durstep_pimpl: public virtual durstep_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  simstep (const ::xml_schema::duration&);

  virtual void
  simduration (const ::xml_schema::duration&);

  virtual void
  post_durstep ();
  
  public:
    
  int sim_step;
  int sim_duration;
    
};

class ORB_initstate_pimpl: public virtual ORB_initstate_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  Initime (const ::xml_schema::date_time&);

  virtual void
  Position ();

  virtual void
  Velocity ();

  virtual void
  post_ORB_initstate ();
  
  public:
  
  Vector6d UTCdate;
  Vec3d Pos_vec, Vel_vec;
};

class ATT_initstate_pimpl: public virtual ATT_initstate_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  phi ();

  virtual void
  theta ();

  virtual void
  psi ();

  virtual void
  om_x ();

  virtual void
  om_y ();

  virtual void
  om_z ();

  virtual void
  post_ATT_initstate ();
  
  public:
  
  double phi_in, theta_in, psi_in, om_x_in, om_y_in, om_z_in;
};

class simoptions_pimpl: public virtual simoptions_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  initstate_in_RTN (bool);

  virtual void
  realtime (bool);
  
  virtual void
  realtime_wait ();

  virtual void
  ggrad_on (bool);

  virtual void
  mag_on (bool);

  virtual void
  srp_on (bool);

  virtual void
  drag_on (bool);
  
  virtual void
  nMAX ();

  virtual void
  sunmoon_on (bool);
  
  virtual void
  Drag_Model (const ::std::string&);

  virtual void
  SRP_Model (const ::std::string&);

  virtual void
  AttitudeType (const ::std::string&);
  
  virtual void
  attctrl_on (bool);

  virtual void
  AttCtrlType (const ::std::string&);

  virtual void
  orbctrl_on (bool);

  virtual void
  OrbCtrlType (const ::std::string&);

  virtual void
  post_simoptions ();
  
  public:
  
  bool initstate_in_RTN_in, realtime_in, ggrad_on_in, mag_on_in, drag_on_in, srp_on_in, sunmoon_on_in, attctrl_on_in, orbctrl_on_in;
  double realtime_wait_in;
  int nMAX_in;
  string AttitudeType_in, Drag_Model_in, SRP_Model_in, AttCtrlType_in, OrbCtrlType_in;
};

class InputFiles_pimpl: public virtual InputFiles_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  Orbit_ephemeris (const ::std::string&);
  
  virtual void
  Attitude_ephemeris (const ::std::string&);
  
  virtual void
  Data_path (const ::std::string&);

  virtual void
  Planet_ephemeris (const ::std::string&);

  virtual void
  EOP_parameters (const ::std::string&);

  virtual void
  PCK_data (const ::std::string&);

  virtual void
  Leap_second (const ::std::string&);

  virtual void
  Gravity_model (const ::std::string&);

  virtual void
  Atmospheric_model (const ::std::string&);

  virtual void
  Magnetic_model (const ::std::string&);

  virtual void
  name (const ::std::string&);

  virtual void
  post_InputFiles ();
  
  
  public:
  
  string Orbit_ephemeris_in, Attitude_ephemeris_in, Data_path_in, planetephemeris_in, eop_in, pck_data_in, leapsecond_in, magn_model_in, gravityfield_in, atmosphere_in;
};

class OutputFiles_pimpl: public virtual OutputFiles_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  Orbit_ephemeris (const ::std::string&);
  
  virtual void
  Attitude_ephemeris (const ::std::string&);

  virtual void
  Sensor_output (const ::std::string&);

  virtual void
  Torques (const ::std::string&);
  
  virtual void
  Accelerations (const ::std::string&);

  virtual void
  name (const ::std::string&);

  virtual void
  post_OutputFiles ();
  
  public:
    
  string orbfile_name_in, attfile_name_in, csv_attfile_name_in, sensors_filename_in, csv_torques_name_in, csv_accelerations_name_in;
};

class SensorsActuators_pimpl: public virtual SensorsActuators_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  subsystem_on (bool);

  virtual void
  constparam ();

  virtual void
  auxparam ();

  virtual void
  opslimit ();

  virtual void
  accuracy ();

  virtual void
  SC2SYS_matrix ();

  virtual void
  name (const ::std::string&);

  virtual void
  post_SensorsActuators ();
  
  int auxparamind, constparamind, opslimitsind, accuracyind;
};

class constparam_pimpl: public virtual constparam_pskel,
  public ::xml_schema::double_pimpl
{
  public:
  virtual void
  pre ();

  virtual void
  name (const ::std::string&);

  virtual void
  unit (const ::std::string&);

  virtual void
  post_constparam ();
};

class auxparam_pimpl: public virtual auxparam_pskel,
  public ::xml_schema::double_pimpl
{
  public:
  virtual void
  pre ();

  virtual void
  name (const ::std::string&);

  virtual void
  unit (const ::std::string&);

  virtual void
  post_auxparam ();
};

class opslimit_pimpl: public virtual opslimit_pskel,
  public ::xml_schema::double_pimpl
{
  public:
  virtual void
  pre ();

  virtual void
  name (const ::std::string&);

  virtual void
  unit (const ::std::string&);

  virtual void
  post_opslimit ();
};

class accuracy_pimpl: public virtual accuracy_pskel,
  public ::xml_schema::double_pimpl
{
  public:
  virtual void
  pre ();

  virtual void
  name (const ::std::string&);

  virtual void
  unit (const ::std::string&);

  virtual void
  post_accuracy ();
};

class Maneuvers_pimpl: public virtual Maneuvers_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  Man ();

  virtual void
  post_Maneuvers ();
  
  //SC::maneuver maneuver_in;
  vector<SC::maneuver> all_maneuvers; // struct maneuver defined in VarTypes.h
  //vector<SC::maneuver>::iterator man_ind;
  
};

class Man_pimpl: public virtual Man_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  maneuver_on (bool);

  virtual void
  init_time (double);

  virtual void
  duration (double);

  virtual void
  ManVec ();

  virtual void
  name ();

  virtual void
  post_Man ();
};

class CompParameters_pimpl: public virtual CompParameters_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  durstep ();

  virtual void
  Payload ();

  virtual void
  Spacecraft ();

  virtual void
  Compoptions ();

  virtual void
  post_CompParameters ();
};

class Payload_pimpl: public virtual Payload_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  FOV_cross ();

  virtual void
  FOV_along ();

  virtual void
  post_Payload ();
  
  double FOV_cross_in, FOV_along_in;
};

class Spacecraft_pimpl: public virtual Spacecraft_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  SC_start (unsigned long long);

  virtual void
  SC_end (unsigned long long);

  virtual void
  PL_start (unsigned long long);

  virtual void
  PL_end (unsigned long long);

  virtual void
  post_Spacecraft ();
  
  int SC_start_in, SC_end_in, PL_start_in, PL_end_in;
};

class Compoptions_pimpl: public virtual Compoptions_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  TGs_on (bool);

  virtual void
  GSs_on (bool);

  virtual void
  TGs_grid_on (bool);
  
  virtual void
  Eclipse_on (bool);

  virtual void
  post_Compoptions ();
  
  bool TGs_on_in, GSs_on_in, TGs_grid_on_in, Eclipse_on_in;
};

class TGs_pimpl: public virtual TGs_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  TGs_grid ();

  virtual void
  TGs_list ();

  virtual void
  post_TGs ();
};

class TGs_grid_pimpl: public virtual TGs_grid_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  minlon ();

  virtual void
  maxlon ();

  virtual void
  minlat ();

  virtual void
  maxlat ();

  virtual void
  gridstep ();

  virtual void
  post_TGs_grid ();
  
  double minlon_in, maxlon_in, minlat_in, maxlat_in, gridstep_in;
};

class TGs_list_pimpl: public virtual TGs_list_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  TG ();

  virtual void
  post_TGs_list ();
  
  ground::TG TGs_list_in[1000];
  int TG_ind;
};

class TG_pimpl: public virtual TG_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  lon ();

  virtual void
  lat ();

  virtual void
  alt ();

  virtual void
  name (const ::std::string&);

  virtual void
  post_TG ();
};

class GSs_pimpl: public virtual GSs_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  GS ();

  virtual void
  post_GSs ();
  
  ground::GS GSs_list_in[1000];
  int GS_ind;
};

class GS_pimpl: public virtual GS_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  lon ();

  virtual void
  lat ();

  virtual void
  alt ();

  virtual void
  minelev ();

  virtual void
  name (const ::std::string&);

  virtual void
  post_GS ();
};

class EventsInputFiles_pimpl: public virtual EventsInputFiles_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  Orbit_ephemeris_path (const ::std::string&);

  virtual void
  Orbit_ephemeris_rootname (const ::std::string&);
  
  virtual void
  Data_path (const ::std::string&);

  virtual void
  Planet_ephemeris (const ::std::string&);

  virtual void
  EOP_parameters (const ::std::string&);

  virtual void
  PCK_data (const ::std::string&);

  virtual void
  Leap_second (const ::std::string&);

  virtual void
  name (const ::std::string&);

  virtual void
  post_EventsInputFiles ();
  
  public:
  
  string Orbit_ephemeris_path_in, Orbit_ephemeris_rootname_in, Data_path_in, planetephemeris_in, eop_in, pck_data_in, leapsecond_in;
};

class EventsOutputFiles_pimpl: public virtual EventsOutputFiles_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  TG_contacts (const ::std::string&);

  virtual void
  GS_contacts (const ::std::string&);
  
  virtual void
  Eclipse_times (const ::std::string&);

  virtual void
  name (const ::std::string&);

  virtual void
  post_EventsOutputFiles ();
  
  public:
  
  string TG_contacts_in, GS_contacts_in, Eclipse_times_in;
};

class license_pimpl: public virtual license_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  licenseName (const ::std::string&);

  virtual void
  licenseURL (const ::std::string&);

  virtual void
  post_license ();
};

class name_pimpl: public virtual name_pskel,
  public ::xml_schema::string_pimpl
{
  public:
  virtual void
  pre ();

  virtual void
  post_name ();
};

#endif // SIMPARAM_SCHEMA_PIMPL_HXX
