//------------------------------------------------------------------------------
//
// HIL_interface.cpp
// 
// Purpose: 
//
//   Utilities to read and write from host computer 
//
// Last modified:
//
//   2016/11/24  PS  Created
//
//
// (c) Sergio De Florio
//
//------------------------------------------------------------------------------

#define BOOST_BIND_GLOBAL_PLACEHOLDERS

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/iostreams/stream.hpp>
#include <sstream>
#include <cassert>
#include <iostream>
#include <boost/asio.hpp>
#include <string>
//#include <string.h>

#include <HIL_interface.hpp>
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

//void print(boost::property_tree::ptree const& pt)
//    {
//    using boost::property_tree::ptree;
//    ptree::const_iterator end = pt.end();
//    for (ptree::const_iterator it = pt.begin(); it != end; ++it)
//        {
//        std::cout << it->first << ": " << it->second.get_value<std::string>() << std::endl;
//        print(it->second);
//        }
//    }
//
//
//int HIL_step(pt::ptree &send_json, pt::ptree &rcv_json)
//    {
//    tcp::iostream socket("localhost", "9999");
//    if (!socket) {
//      std::cout << "Unable to connect: " << socket.error().message() << std::endl;
//      return -1;
//    }
//    std::stringstream json_buf;
//    write_json(json_buf, send_json);
//    socket << json_buf.str() << '\0';
//    pt::read_json(socket, rcv_json);
//
//    return 0;
//    }


Vec3d actions_vec(pt::ptree const& tree, pt::ptree::key_type const& key)
    {
    Vec3d r;
    int ind = 0;
    
    for (auto& item : tree.get_child(key))
        {
        r(ind) = item.second.get_value<double>();
        ind++;
        }
        
    return(r);
    };
    
    
void send_receiveTCs(Vec3d& magnetotorquer, Vec3d& reaction_wheel, sensorTCs* sensTCs)
    {
    pt::ptree send_json;
    pt::ptree rcv_json;
  
    //send_json.put<int>("foo", 23);
    
    send_json.put<unsigned int>("UnixTime", sensTCs->UnixTime);
    
    for(int i = 0; i < 10; i++)
        {
        string TC_name = "CssRaw" + to_string(i+1);
        send_json.put<unsigned int>(TC_name, sensTCs->CssRaw(i));
        }
        
    send_json.put<int>("SunRawX", sensTCs->Cam1Raw(0));
    send_json.put<int>("SunRawY", sensTCs->Cam1Raw(1));
    send_json.put<unsigned int>("SunBusy", sensTCs->Cam1Busy);
    send_json.put<unsigned int>("SunResult", sensTCs->Cam1Result);
    
    send_json.put<int>("NadirRawX", sensTCs->Cam2Raw(0));
    send_json.put<int>("NadirRawY", sensTCs->Cam2Raw(1));
    send_json.put<unsigned int>("NadirBusy", sensTCs->Cam2Busy);
    send_json.put<unsigned int>("NadirResult", sensTCs->Cam2Result);
    
    send_json.put<int>("MagRawX", sensTCs->MagRaw(0));
    send_json.put<int>("MagRawY", sensTCs->MagRaw(1));
    send_json.put<int>("MagRawZ", sensTCs->MagRaw(2));
    
    send_json.put<int>("RateRawX", sensTCs->RateRaw(0));
    send_json.put<int>("RateRawY", sensTCs->RateRaw(1));
    send_json.put<int>("RateRawZ", sensTCs->RateRaw(2));
    
    send_json.put<int>("WheelRawX", sensTCs->WheelRaw(0));
    send_json.put<int>("WheelRawY", sensTCs->WheelRaw(1));
    send_json.put<int>("WheelRawZ", sensTCs->WheelRaw(2));
    
    send_json.put<int>("Star1CameraX", sensTCs->Star1Camera(0));
    send_json.put<int>("Star1CameraY", sensTCs->Star1Camera(1));
    send_json.put<int>("Star1CameraZ", sensTCs->Star1Camera(2));
    
    send_json.put<int>("Star1InertialX", sensTCs->Star1Inertial(0));
    send_json.put<int>("Star1InertialY", sensTCs->Star1Inertial(1));
    send_json.put<int>("Star1InertialZ", sensTCs->Star1Inertial(2));
        
    send_json.put<int>("Star2CameraX", sensTCs->Star2Camera(0));
    send_json.put<int>("Star2CameraY", sensTCs->Star2Camera(1));
    send_json.put<int>("Star2CameraZ", sensTCs->Star2Camera(2));
    
    send_json.put<int>("Star2InertialX", sensTCs->Star2Inertial(0));
    send_json.put<int>("Star2InertialY", sensTCs->Star2Inertial(1));
    send_json.put<int>("Star2InertialZ", sensTCs->Star2Inertial(2));
    
    send_json.put<int>("Star3CameraX", sensTCs->Star3Camera(0));
    send_json.put<int>("Star3CameraY", sensTCs->Star3Camera(1));
    send_json.put<int>("Star3CameraZ", sensTCs->Star3Camera(2));
    
    send_json.put<int>("Star3InertialX", sensTCs->Star3Inertial(0));
    send_json.put<int>("Star3InertialY", sensTCs->Star3Inertial(1));
    send_json.put<int>("Star3InertialZ", sensTCs->Star3Inertial(2));

    
    //vector<float> rate = {1, 2, 3};
    //Vec3d rate(1, 2, 3);
    //
    //pt::ptree rate_node;
    ////for (auto &rate_i : rate)
    //for(int i = 0; i < rate.size(); i++)
    //  {
    //  pt::ptree rate_i_node;
    //  rate_i_node.put("", rate(i));
    //  rate_node.push_back(make_pair("", rate_i_node));
    //  }
    //send_json.add_child("rate", rate_node);
    
    cout << endl << "sending: " << endl;
    print(send_json);
  
    HIL_step(send_json, rcv_json);
  
    //Vec3d magnetotorquer, reaction_wheel;
    
    magnetotorquer = actions_vec(rcv_json, "magnetotorquer");
    cout << "magnetotorquer:" << endl;
    cout << magnetotorquer(0) << endl;
    cout << magnetotorquer(1) << endl;
    cout << magnetotorquer(2) << endl;
    
    reaction_wheel = actions_vec(rcv_json, "reaction_wheel");
    cout << "reaction_wheel:" << endl;
    cout << reaction_wheel(0) << endl;
    cout << reaction_wheel(1) << endl;
    cout << reaction_wheel(2) << endl;
    };
    
    
    
  void send_receiveTCs(const Ref<const VectorXi>& sensorTC, const string& TC_name)
    {
    pt::ptree send_json;
    pt::ptree rcv_json;
  
    //send_json.put<int>("foo", 23);
  
    //vector<float> rate = {1, 2, 3};
    //Vec3d sensorTC(1, 2, 3);
    
    pt::ptree rate_node;
    //for (auto &rate_i : sensorTC)
    for(int i = 0; i < sensorTC.size(); i++)
      {
      pt::ptree rate_i_node;
      rate_i_node.put("", sensorTC(i));
      rate_node.push_back(make_pair("", rate_i_node));
      }
    send_json.add_child(TC_name, rate_node);
  
    cout << endl << "sending: " << endl;
    print(send_json);
  
    HIL_step(send_json, rcv_json);
    };
    
    
    void send_receiveTCs(const Ref< const Matrix<unsigned int, Dynamic, 1> >& sensorTC, const string& TC_name)
    {
    pt::ptree send_json;
    pt::ptree rcv_json;
  
    //send_json.put<int>("foo", 23);
  
    //vector<float> rate = {1, 2, 3};
    //Vec3d sensorTC(1, 2, 3);
    
    pt::ptree rate_node;
    //for (auto &rate_i : sensorTC)
    for(int i = 0; i < sensorTC.size(); i++)
      {
      pt::ptree rate_i_node;
      rate_i_node.put("", sensorTC(i));
      rate_node.push_back(make_pair("", rate_i_node));
      }
    send_json.add_child(TC_name, rate_node);
  
    cout << endl << "sending: " << endl;
    print(send_json);
  
    HIL_step(send_json, rcv_json);
    };
    
    
    
    
    void send_receiveTCs(unsigned int sensorTC, const string& TC_name)
    {
    pt::ptree send_json;
    pt::ptree rcv_json;
  
    //send_json.put<int>("foo", 23);
  
    //vector<float> rate = {1, 2, 3};
    //Vec3d sensorTC(1, 2, 3);
    
    pt::ptree rate_node;
    //for (auto &rate_i : sensorTC)
    pt::ptree rate_i_node;
    rate_i_node.put("", sensorTC);
    rate_node.push_back(make_pair("", rate_i_node));
    
    send_json.add_child(TC_name, rate_node);
  
    cout << endl << "sending: " << endl;
    print(send_json);
  
    HIL_step(send_json, rcv_json);
    };  
    
    
    
