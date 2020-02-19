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

#include <corbo-systems/time_value_buffer.h>
#include <limits>

namespace corbo {

void TimeValueBuffer::getValues(double ts, double dt, std::vector<std::pair<double, Eigen::VectorXd>>& useq_out)
{
    useq_out.clear();
    if (_ucache.empty())
    {
        if (_uinit.rows() == 0) PRINT_ERROR_NAMED("_uinit not initialized. Call setInitialValue().");
        _ucache.emplace_back(std::numeric_limits<double>::lowest(), _uinit);
    }

    // now check which controls we need to output
    // find index corresponding to the current time in the control sequence which considers deadtime
    int start_idx = 0;
    for (; start_idx < (int)_ucache.size(); ++start_idx)
    {
        if (ts < _ucache[start_idx].first) break;
    }
    start_idx -= 1;  // note, the time stamps in the cache are the beginning of each interval!

    // now export the control sequence of total length dt
    double cur_t = ts;
    int idx      = start_idx;
    for (; idx < (int)_ucache.size() - 1; ++idx)
    {
        double dti = _ucache[idx + 1].first - cur_t;

        if (dti + cur_t < ts + dt)
        {
            useq_out.emplace_back(dti, _ucache[idx].second);
        }
        else
        {
            useq_out.emplace_back(ts + dt - cur_t, _ucache[idx].second);
            break;
        }

        cur_t = _ucache[idx + 1].first;
    }

    if (idx == (int)_ucache.size() - 1)
    {
        if (!useq_out.empty()) cur_t = _ucache.back().first;

        useq_out.emplace_back(ts + dt - cur_t, _ucache.back().second);
    }

    // cleanup values that are not longer required
    if (start_idx - 1 > 0)
    {
        _ucache.erase(_ucache.begin(), _ucache.begin() + start_idx - 1);  // TODO(roesmann): better std::deque?
    }

#ifndef NDEBUG
    // check if dt is correct
    double dt_test = 0;
    for (int i = 0; i < useq_out.size(); ++i) dt_test += useq_out[i].first;
    PRINT_ERROR_COND(std::abs(dt_test - dt) > 1e-10, "Deadtime: Computed dt (" << dt_test << ") does not match original dt (" << dt << ")");
#endif
}

void TimeValueBuffer::appendValues(double t, const Eigen::Ref<const Eigen::VectorXd>& u)
{
    if (_ucache.empty())
    {
        _ucache.emplace_back(std::numeric_limits<double>::lowest(), _uinit);
    }

    if (t < _ucache.back().first) PRINT_WARNING_NAMED("t can not be less than the last value. We can not change the past.");

    if (_ucache.back().first == t)  // make sure that t is unique
    {
        _ucache.back().second = u;
    }
    else
    {
        _ucache.emplace_back(t, u);
    }
}

}  // namespace corbo
