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

#include <SpaceEnvironment.h>
#include <Constants.h>

using namespace std;
using namespace math;
using namespace constants;

namespace spaceenvironment
    {
    //------------------------------------------------------------------------------
    // SPACEENV implementation
    //------------------------------------------------------------------------------
    int SPACEENV::n_max = 0;
    // Constructor
    //------------------------------------------------------------------------------
    SPACEENV::SPACEENV()
       {
        modelname = " ";
        modelfilepath = " ";
        refsys = " ";
        Area = 0.0;
        CF = 0.0;
        rho1 = 0.0;
        rho2 = 0.0;
        rho3 = 0.0;
        n_max = 0;
       };
    //------------------------------------------------------------------------------
    // Destructor
    //------------------------------------------------------------------------------
    SPACEENV::~SPACEENV() {};
    //------------------------------------------------------------------------------
    // Method setfilespath(const string model_filename)
    //------------------------------------------------------------------------------
    /**
     * Set model's file path. In this version of the function, the name of the model's
     * file is put in class member 'modelname'
     *
     * @param model_filepath   Path of model's file
     *
     */
    //------------------------------------------------------------------------------  
    void SPACEENV::setfilespath(const string model_filepath)
                                {
                                modelfilepath = model_filepath;
                                // Find model name
                                size_t last_slash, last_point;
                                last_slash = modelfilepath.rfind("/");
                                last_point = modelfilepath.rfind(".");
                                int len = last_point - last_slash - 1;
                                
                                modelname = modelfilepath.substr(last_slash + 1, len);
                                };
    //------------------------------------------------------------------------------
    // Method void setfilespath(const string model_filename, const string model_filepath)
    //------------------------------------------------------------------------------
    /**
     * Set model file's folder path and name
     *
     * @param model_name       Name of the model
     * @param model_filepath   Path of the folder containing the model's file
     *
     */
    //------------------------------------------------------------------------------
    void SPACEENV::setfilespath(const string model_name, const string model_filepath)
                                {
                                modelname = model_name;
                                modelfilepath = model_filepath;
                                };
    //------------------------------------------------------------------------------
    // Method void SetReferenceFrame(const string ref_sys)
    //------------------------------------------------------------------------------
    /**
     * Set the reference frame in which the field's vector is computed
     *
     * @param ref_sys       Name of reference frame
     *
     */
    //------------------------------------------------------------------------------
    void SPACEENV::SetReferenceFrame(const string ref_sys) { refsys = ref_sys;};
    
}; // End of namespace spaceenvironment



