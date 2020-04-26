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

#ifndef CONTACTS_H_
#define CONTACTS_H_

#include <string>
#include <fstream>
#include <Eigen/Core>
#include <Eigen/Geometry> 
#include <VarTypes.h>

#include <Solarsys.h>
      
using namespace std;
using namespace math;
using namespace Eigen;

using namespace solarsystem;

//------------------------------------------------------------------------------
// Group of functions for spacecraft contacts with targets and ground stations
//------------------------------------------------------------------------------

// Compute sensor contacts of one or more spacecraft with one or more targets
void TGsContacts(int SC_num, int PL_num, Vec4d& orbel, const MatrixXd& orbstateECEF, VectorXd& TGs_grid_lons, VectorXd& TGs_grid_lats, ground::TG* TGs_list, string output_path, double FOV, double FOVdir, double comp_duration, int simstep, bool TGs_on, bool TGs_grid_on, ofstream& AllContacts_file);
// Compute contacts of one or more spacecraft with one or more ground stations
void GSsContacts(int SC_num, int PL_num, Vec4d& orbel, const MatrixXd& orbstateECEF, ground::GS* GSs_list, string output_path, double comp_duration, int simstep, ofstream& AllContacts_file);
// Compute umbra and penumbra entry and exit times of one or more spacecraft
void Umbras(int SC_num, int PL_num, SOLSYS Solar, const MatrixXd& orbpos, string output_path, double comp_duration, int simstep, ofstream& AllUmbras_file);

//// Compute sensor contacts of one or more spacecraft with one or more targets
//Matrix<double,Dynamic,6> TGsContacts(int SC_num, int PL_num, Vec4d& orbel, const MatrixXd& orbstateECEF, VectorXd& TGs_grid_lons, VectorXd& TGs_grid_lats, ground::TG* TGs_list, string output_path, double FOV, double FOVdir, double duration, int simstep, bool TGs_grid_on, ofstream& AllContacts_file);
//// Compute contacts of one or more spacecraft with one or more ground stations
//Matrix<double,Dynamic,6> GSsContacts(int SC_num, int PL_num, Vec4d& orbel, const MatrixXd& orbstateECEF, ground::GS* GSs_list, string output_path, double comp_duration, int simstep, ofstream& AllContacts_file);

#endif // PROPAGATION_H_