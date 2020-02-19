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

#ifndef SRC_CORE_INCLUDE_CORBO_CORE_TYPES_H_
#define SRC_CORE_INCLUDE_CORBO_CORE_TYPES_H_

#include <Eigen/Core>
#include <limits>
#include <memory>

namespace corbo {

namespace property {

//! Dynamic property
constexpr const int DYNAMIC = -32767;

//! Inherit property
constexpr const int INHERITED = -32766;
}  // namespace property

//! Maximum double value
constexpr const double CORBO_MAX_DBL = std::numeric_limits<double>::max();
//! Minimum (negative) double value
constexpr const double CORBO_MIN_DBL = std::numeric_limits<double>::lowest();
//! Maximum integer value
constexpr const int CORBO_MAX_INT = std::numeric_limits<int>::max();
//! Minimum (negative) integer value
constexpr const int CORBO_MIN_INT = std::numeric_limits<int>::lowest();

//! Representation for infinity (double version)
constexpr const double CORBO_INF_DBL = 2e30;  // 0.9 * CORBO_MAX_DBL;
//! Representation for infinity (integer version)
constexpr const int CORBO_INF_INT = CORBO_MAX_INT - 1000;

}  // namespace corbo

#endif  // SRC_CORE_INCLUDE_CORBO_CORE_TYPES_H_
