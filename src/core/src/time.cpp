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

#include <corbo-core/time.h>

namespace corbo {

// implementation based on ros::Rate
// (http://docs.ros.org/diamondback/api/rostime/html/rate_8cpp_source.html#l00040)
bool Rate::sleep()
{
    Time expected_end = _start + _cycle_time;
    Time actual_end   = Time::now();

    // backwards jump in time
    if (actual_end < _start) expected_end = actual_end + _cycle_time;

    Duration sleep_duration = expected_end - actual_end;
    _last_cycle_time        = actual_end - _start;

    _start = expected_end;

    // only sleep if positive sleep duration
    if (sleep_duration <= Duration(0))
    {
        if (actual_end > expected_end + _cycle_time) _start = actual_end;
        return !(sleep_duration < Duration(0));
    }
    sleep_duration.sleep();
    return true;
}

Time Duration::toTime(double basis_time) { return Time(basis_time) + *this; }

}  // namespace corbo
