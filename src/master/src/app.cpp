/*********************************************************************
 *
 *  Software License Agreement
 *
 *  Copyright (c) 2020,
 *  TU Dortmund - Institute of Control Theory and Systems Engineering.
 *  All rights reserved.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  Authors: Christoph Rösmann
 *********************************************************************/

#include <corbo-master/master.h>

#include <corbo-controllers/corbo_controllers.h>
#include <corbo-core/corbo_core.h>
#include <corbo-numerics/corbo_numerics.h>
#include <corbo-observers/corbo_observers.h>
#include <corbo-plants/corbo_plants.h>
#include <corbo-systems/corbo_systems.h>
#include <corbo-tasks/corbo_tasks.h>

#include <algorithm>
#include <string>

bool has_option(char** argv_begin, char** argv_end, const std::string& option) { return std::find(argv_begin, argv_end, option) != argv_end; }

std::string get_option(char** argv_begin, char** argv_end, const std::string& option)
{
    char** it = std::find(argv_begin, argv_end, option);
    if (it != argv_end && ++it != argv_end) return std::string(*it);
    return std::string();
}

void print_help()
{
    PRINT_INFO("=== corbo-master ===");
    PRINT_INFO("The following options are available:");
    PRINT_INFO("-s\tSet gRPC server hostname, e.g. -s localhost:50051");
    PRINT_INFO("-m\tLoad protobuf message file of type corboParameters, e.g. -m path/to/file.cparams");
}

int main(int argc, char** argv)
{
    // parse arguments for options
    if (has_option(argv, argv + argc, "-h"))
    {
        print_help();
        return 0;
    }

    std::string hostname           = get_option(argv, argv + argc, "-s");
    if (hostname.empty()) hostname = "localhost:50051";

    std::string proto_msg_path = get_option(argv, argv + argc, "-m");

    corbo::Master master;

    if (!proto_msg_path.empty())
    {
        if (master.loadFromFile(proto_msg_path))
        {
            PRINT_INFO("Loaded default parameters from file " << proto_msg_path);
        }
        else
        {
            PRINT_ERROR("Cannot load default parameters from file " << proto_msg_path);
        }
    }

    PRINT_DEBUG_WARN("corbo IS COMPILED WITH DEBUG FLAGS: PLEASE EXPECT LARGE COMPUTATION TIMES!");

    master.start(hostname);
    return 0;
}
