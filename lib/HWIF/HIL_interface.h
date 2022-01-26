//------------------------------------------------------------------------------
//
// HIL_interface.h
// 
// Purpose: 
//
//   Utilities to read and write from host computer 
//
// Last modified:
//
//   2016/11/24  PS  Created
//   2016/12/13  SDF  Added function send_sensorsTCs
//
// (c) Sergio De Florio
//
//------------------------------------------------------------------------------

#ifndef HILIF_H_
#define HILIF_H_

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/iostreams/stream.hpp>
#include <sstream>
#include <cassert>
#include <iostream>
#include <boost/asio.hpp>
#include <string>
//#include <string.h>

#include <VarTypes.h>

using namespace math;
using namespace Eigen;
using namespace std;
/*
sim -> interface: sensor values (TCP:Json+\0)
interface -> adcs: trigger HIL loop
interface -> adcs: wait for complete
interface -> adcs: get telemetry
adcs -> interface: actuator telemetry
interface -> sim: actuator values (TCP:Json)
*/

//using boost::asio::ip::tcp;
namespace pt = boost::property_tree;

//void print(boost::property_tree::ptree const& pt);

//int HIL_step(pt::ptree &send_json, pt::ptree &rcv_json);

Vec3d actions_vec(pt::ptree const& tree, pt::ptree::key_type const& key);

void send_receiveTCs(Vec3d& magnetotorquer, Vec3d& reaction_wheel, sensorTCs* sensTCs);

//void send_receiveTCs(const Ref<const VectorXi>& sensorTC, const string& TC_name);
//void send_receiveTCs(const Ref< const Matrix<unsigned int, Dynamic, 1> >& sensorTC, const string& TC_name);
//void send_receiveTCs(unsigned int sensorTC, const string& TC_name);

#endif // HILIF_H_
