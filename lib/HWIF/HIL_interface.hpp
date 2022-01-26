#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/iostreams/stream.hpp>
#include <sstream>
#include <cassert>
#include <iostream>
#include <boost/asio.hpp>
#include <string>

#include <string.h>

/*
sim -> interface: sensor values (TCP:Json+\0)
interface -> adcs: trigger HIL loop
interface -> adcs: wait for complete
interface -> adcs: get telemetry
adcs -> interface: actuator telemetry
interface -> sim: actuator values (TCP:Json)
*/


using boost::asio::ip::tcp;
namespace pt = boost::property_tree;


void print(boost::property_tree::ptree const& pt)
{
    using boost::property_tree::ptree;
    ptree::const_iterator end = pt.end();
    for (ptree::const_iterator it = pt.begin(); it != end; ++it) {
        std::cout << it->first << ": " << it->second.get_value<std::string>() << std::endl;
        print(it->second);
    }
}


int HIL_step(pt::ptree &send_json, pt::ptree &rcv_json)
{
    tcp::iostream socket("localhost", "9999");
    if (!socket) {
      std::cout << "Unable to connect: " << socket.error().message() << std::endl;
      return -1;
    }
    std::stringstream json_buf;
    write_json(json_buf, send_json);
    socket << "[\"adcs_step\", " << json_buf.str() << "]" << '\0';
    pt::read_json(socket, rcv_json);

    return 0;
}
