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

#ifndef SRC_CORE_INCLUDE_CORBO_CORE_VALUE_COMPARISON_H_
#define SRC_CORE_INCLUDE_CORBO_CORE_VALUE_COMPARISON_H_

#include <algorithm>
#include <cmath>
#include <complex>

namespace corbo {

/**
 * @brief Check if two doubles are similar (but with an absolute threshold!)
 * @ingroup core
 * @warning better use relative comparisons to be much more robust against different magnitutes
 * @param a         first value
 * @param b         second value
 * @param epsilon   minimum treshold to state that a and b are similar (|a-b| < epsilon)
 * @return \c true if <c> fabs(a-b)< epsilon </c>, otherwise \c false
 */
inline bool approx_equal_abs(double a, double b, double epsilon = 1e-6) { return std::abs(a - b) < epsilon; }

// The double/float comparison definitions can be found in "The art of computer programming by Knuth"

/**
 * @brief Check if two doubles are appoximately equal
 * @ingroup core
 * @warning use approx_zero() if you want to compare against zero
 * @param a         first value
 * @param b         second value
 * @param epsilon   minimum treshold to state that a and b are similar (|a-b| <= epsilon*max(a,b))
 * @return \c true if <c> fabs(a-b)<= epsilon*max(a,b) </c>, otherwise \c false
 */
inline bool approx_equal(double a, double b, double epsilon = 1e-6) { return std::abs(a - b) <= std::max(std::abs(a), std::abs(b)) * epsilon; }

/**
 * @brief Check if two complex numbers are appoximately equal by comparing their real and imaginary parts
 * @ingroup core
 * @param a         first complex value
 * @param b         second complex value
 * @param epsilon   minimum treshold to state that a and b are similar (|a-b| <= epsilon*max(a,b))
 * @return \c true if <c> fabs(a-b)<= epsilon*max(a,b) </c>, otherwise \c false
 */
inline bool approx_equal(const std::complex<double>& a, const std::complex<double>& b, double epsilon = 1e-6)
{
    return approx_equal(a.real(), b.real(), epsilon) && approx_equal(a.imag(), b.imag(), epsilon);
}

/**
 * @brief Check if two doubles are essentially equal
 * @ingroup core
 * @warning use approx_zero() if you want to compare against zero
 * @param a         first value
 * @param b         second value
 * @param epsilon   minimum treshold to state that a and b are similar (|a-b| <= epsilon*min(a,b))
 * @return \c true if <c> fabs(a-b)<= epsilon*min(a,b) </c>, otherwise \c false
 */
inline bool essentially_equal(double a, double b, double epsilon = 1e-6) { return std::abs(a - b) <= std::min(std::abs(a), std::abs(b)) * epsilon; }

/**
 * @brief Check if two compelx numbers are essentially equal by comparing their real and imaginary parts
 * @ingroup core
 * @param a         first complex value
 * @param b         second complex value
 * @param epsilon   minimum treshold to state that a and b are similar (|a-b| <= epsilon*min(a,b))
 * @return \c true if <c> fabs(a-b)<= epsilon*min(a,b) </c>, otherwise \c false
 */
inline bool essentially_equal(const std::complex<double>& a, const std::complex<double>& b, double epsilon = 1e-6)
{
    return essentially_equal(a.real(), b.real(), epsilon) && essentially_equal(a.imag(), b.imag(), epsilon);
}

/**
 * @brief Check if a double is appoximately zero
 * @ingroup core
 * @param val       first value
 * @param epsilon   minimum treshold to state that val is similar (|a| <= epsilon)
 * @return \c true if <c> fabs(val)<= epsilon) </c>, otherwise \c false
 */
inline bool approx_zero(double val, double epsilon = 1e-10) { return std::abs(val) <= epsilon; }

/**
 * @brief Check if a value is definetly greater than another value
 * @ingroup core
 * @param a         first value
 * @param b         second value
 * @param epsilon   minimum treshold to state that a and b are similar (|a-b| > epsilon*max(a,b))
 * @return \c true if <c> fabs(a-b)> epsilon*min(a,b) </c>, otherwise \c false
 */
inline bool definitely_greater_than(double a, double b, double epsilon = 1e-6) { return (a - b) > std::max(std::abs(a), std::abs(b)) * epsilon; }

/**
 * @brief Check if a value is definetly less than another value
 * @ingroup core
 * @param a         first value
 * @param b         second value
 * @param epsilon   minimum treshold to state that a and b are similar (|b-a| > epsilon*min(a,b))
 * @return \c true if <c> fabs(b-a)> epsilon*max(a,b) </c>, otherwise \c false
 */
inline bool definitely_less_than(double a, double b, double epsilon = 1e-6) { return (b - a) > std::min(std::abs(a), std::abs(b)) * epsilon; }

// almostEqual2sComplement(double A, double B, int maxUlps)
// {
//   assert(maxUlps > 0 && maxUlps < 4 * 1024 * 1024);
//   int64_t aLong = *reinterpret_cast<int64_t*>( &A );
//   if (aLong < 0) aLong = 0x8000000000000000 - aLong; int64_t bLong = *reinterpret_cast<int64_t*>( &B );
//   if (bLong < 0) bLong = 0x8000000000000000 - bLong; int64_t longDiff = (aLong - bLong) & 0x7FFFFFFFFFFFFFFF;
//   if (longDiff <= maxUlps) return true;
//   return false;
// }

}  // namespace corbo

#endif  // SRC_CORE_INCLUDE_CORBO_CORE_VALUE_COMPARISON_H_
