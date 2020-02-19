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

#ifndef SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_UTILITIES_H_
#define SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_UTILITIES_H_

#include <cmath>

namespace corbo {

namespace util {

/**
 * @brief Continuous approximation of the max-operator
 *
 * Approximates \c max(x,y).
 * @param[in] x    first value
 * @param[in] y    second value
 * @param[in] eps  coefficient
 * @returns approximation of max(x,y)
 */
inline double max_approx_continuous(double x, double y, double eps) { return 0.5 * (std::sqrt(std::pow(x - y, 2) + eps * eps) + x + y); }

}  // namespace util

}  // namespace corbo

#endif  // SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_UTILITIES_H_
