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

#ifndef SRC_CORE_INCLUDE_CORBO_CORE_UTILITIES_H_
#define SRC_CORE_INCLUDE_CORBO_CORE_UTILITIES_H_

#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>

namespace corbo {

namespace util {

#if __cplusplus > 201402L
using std::clamp;
#else
/**
 * @brief Bound a value to a lower and upper limit if exceeded
 * @param[in]   v   value
 * @param[in]   lo  lower bound
 * @param[in]   hi  upper bound
 * @tparam T    value and bound type which must be less-comparable (<)
 * @returns v if lo <= v <= hi, lo if v < lo and hi if v > hi.
 */
template <class T>
constexpr const T& clamp(const T& v, const T& lo, const T& hi)
{
    return std::max(lo, std::min(v, hi));
}
#endif

/**
 * @brief Check if a value is inside the interval [lo, hi]
 * @param[in]   v   value
 * @param[in]   lo  lower bound
 * @param[in]   hi  upper bound
 * @tparam T    value and bound type which must be comparable (<, >)
 * @returns true if lo <= v <= hi, false otherwise.
 */
template <class T>
constexpr const bool is_in_bounds(const T& v, const T& lo, const T& hi)
{
    return !(v < lo) && !(v > hi);
}

/**
 * @brief Check if all components of a container are inside the interval [lo, hi]
 * @param[in]  first    iterator pointing to the begin of a container
 * @param[in]  last     iterator pointing to the end of a container
 * @param[in]  o        lower bound
 * @param[in]  hi       upper bound
 *
 * @tparam T          value and bound type which must be comparable (<, >)
 * @tparam Iterator   Container iterator type
 * @returns true if lo <= v <= hi, false otherwise.
 */
template <class T, class Iterator>
inline const bool is_in_bounds_all(Iterator first, Iterator last, const T& lo, const T& hi)
{
    for (Iterator it = first; it != last; ++it)
        if (*it < lo || *it > hi) return false;
    return true;
}

/**
 * @brief Signum function
 *
 * if \c val < 0: -1,
 * if \c val > 0: 1,
 * if \c val == 0: 0.
 *
 * @param[in] val  value of type T
 * @tparam    T    value type with T(0) as zero value and the type must support "<".
 * @returns integer {-1, 0, 1} according to the rules specified above.
 */
template <class T>
inline int sign(T val)
{
    return (T(0) < val) - (val < T(0));
}

/**
 * @brief Return value a with sign of value b
 *
 * If b > 0 then the result is |a|, else it is -|a|.
 * @param[in] a   value which is subject to sign of b
 * @param[in] b   value from which the sign is extracted
 * @tparam    ValT  Type of a
 * @tparam    SignT Type of b (must support SignT(0))
 * @returns |a| with sign of b and type ValT.
 */
template <class ValT, class SignT>
inline ValT sign(ValT a, SignT b)
{
    return (b > SignT(0) ? std::abs(a) : -std::abs(a));
}

//! Compute the l2 norm of a 2d vector [a,b] as sqrt(a*a+b*b)
inline double l2_norm_2d(double a, double b) { return std::sqrt(a * a + b * b); }

/**
 * @brief Constructs an object of type T and wraps it in a std::unique_ptr.
 *
 * Until C++14 is not widespread, we define our own make_unique.
 * @param args Arbitrary argument list that is compatible with the constructor list required by type \c T
 * @return std::unique_ptr of type \c T constructed with \c args
 */
#if __cplusplus >= 201300
using std::make_unique;
#else
template <typename T, typename... Args>
inline std::unique_ptr<T> make_unique(Args&&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}
#endif

/// @cond 0 // exclude from doxygen
template <typename... A>
struct variadic_temp_equal : std::false_type
{
};

/**
 * @brief This template construct checks whether two variadic template parameter packs have the same type (pairwise).
 *
 * Application for parameter packs A... and B... :
 * @code
 * 	static_assert(variadic_temp_equal<std::tuple<A...>, std::tuple<B...>>::value, "Type mismatch for both template parameter packs.");
 * @endcode
 */
template <typename A1, typename... Aother, typename B1, typename... Bother>
struct variadic_temp_equal<std::tuple<A1, Aother...>, std::tuple<B1, Bother...>>
{
    static const bool value = std::is_same<A1, B1>::value && variadic_temp_equal<std::tuple<Aother...>, std::tuple<Bother...>>::value;
};

template <typename... B>
struct variadic_temp_equal<std::tuple<>, std::tuple<B...>> : std::true_type
{
};
/// @endcond

}  // namespace util

}  // namespace corbo

#endif  // SRC_CORE_INCLUDE_CORBO_CORE_UTILITIES_H_
