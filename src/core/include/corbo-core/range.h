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

#ifndef SRC_CORE_INCLUDE_CORBO_CORE_RANGE_H_
#define SRC_CORE_INCLUDE_CORBO_CORE_RANGE_H_

#include <corbo-core/console.h>

#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace corbo {

class Range
{
 public:
    explicit Range(double single_val, bool force_include_end = false)
        : _start(single_val), _step(0.0), _end(single_val), _include_end(force_include_end)
    {
        update();
    }
    explicit Range(double start, double step, double end, bool force_include_end = false)
        : _start(start), _step(step), _end(end), _include_end(force_include_end)
    {
        update();
    }
    explicit Range(double start, double step, int num, bool force_include_end = false)
        : _start(start), _step(step), _end(start + (double)num * step), _include_end(force_include_end)
    {
        _remainder = 0;
        _n         = num;
    }
    explicit Range(const std::array<double, 3>& interval, bool force_include_end = false)
        : _start(interval[0]), _step(interval[1]), _end(interval[2]), _include_end(force_include_end)
    {
        update();
    }

    inline double getStart() const { return _start; }
    inline double getStep() const { return _step; }
    inline double getEnd() const { return _end; }

    inline bool includeEnd() const { return _include_end && getRemainder() > 1e-8; }

    inline double getEndPlusEps() const { return _end + std::numeric_limits<double>::epsilon(); }
    inline double getRemainder() const { return _remainder; }

    inline int getNumInRange() const { return _n; }  // including start
    
    inline double getProgressFactor(double t) const { return _step == 0.0 ? 0.0 : (t - _start) / ( _end - _start); }
    
    void getGrid(std::vector<double>& values)
    {
        int n = getNumInRange();
        for (int i = 0; i < n; ++i) values.push_back(_start + (double)i * _step);
        if (includeEnd()) values.push_back(_end);
    }

    void getGrid(std::vector<double>& values, double offset)
    {
        int n = getNumInRange();
        for (int i = 0; i < n; ++i) values.push_back(_start + (double)i * _step + offset);
        if (includeEnd()) values.push_back(_end + offset);
    }

    void printGrid()
    {
        int n = getNumInRange();
        PRINT_INFO("=== Grid start ===");
        for (int i = 0; i < n; ++i) PRINT_INFO(_start + (double)i * _step);
        if (includeEnd()) PRINT_INFO(_end);
        PRINT_INFO("=== Grid end ===");
    }

 protected:
    void update()
    {
        if (std::abs(_step) < 1e-15)
        {
            _n         = 1;
            _remainder = 0;
            return;
        }
        _n         = (1 + (int)((_end - _start) / _step));
        _remainder = std::fmod(_end - _start, _step);
    }

 private:
    double _start;
    double _step;
    double _end;
    bool _include_end = false;

    double _remainder;
    int _n;
};

// class Slice
//{
// public:
//    explicit Slice(int start, int stop, int step) : _start(start), _stop(stop), _step(step) {}
//    explicit Slice(int start, int num) : _start(start), _stop(start + num), _step(1) {}

// private:
//    int _start;
//    int _stop;
//    int _step;
//};

}  // namespace corbo

#endif  // SRC_CORE_INCLUDE_CORBO_CORE_RANGE_H_
