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

#ifndef SRC_CORE_INCLUDE_CORBO_CORE_MACROS_H_
#define SRC_CORE_INCLUDE_CORBO_CORE_MACROS_H_

namespace corbo {

// Convert text to string
#define corbo_STRINGIZE_I(text) #text
#define corbo_STRINGIZE(text) corbo_STRINGIZE_I(text)

// Concatenate text
#define corbo_CAT_I(a, b) a##b
#define corbo_CAT(a, b) corbo_CAT_I(a, b)

// Extensions for GTEST
#define EXPECT_EQ_MATRIX(A, B, tol)                                                                                          \
    {                                                                                                                        \
        EXPECT_EQ((A).rows(), (B).rows());                                                                                   \
        EXPECT_EQ((A).cols(), (B).cols());                                                                                   \
        EXPECT_TRUE((A).isApprox(B, tol)) << corbo_STRINGIZE(A) ":\n" << A << "\n" << corbo_STRINGIZE(B) ":\n" << B << ". "; \
    }

#define EXPECT_NEQ_MATRIX(A, B, tol)                                                                                          \
    {                                                                                                                         \
        EXPECT_EQ((A).rows(), (B).rows());                                                                                    \
        EXPECT_EQ((A).cols(), (B).cols());                                                                                    \
        EXPECT_FALSE((A).isApprox(B, tol)) << corbo_STRINGIZE(A) ":\n" << A << "\n" << corbo_STRINGIZE(B) ":\n" << B << ". "; \
    }

#define EXPECT_ZERO_MATRIX(A, tol)                         \
    EXPECT_TRUE(A.isZero(tol)) << corbo_STRINGIZE(A) ":\n" \
                               << A << "\n"                \
                               << ". "

// Static if for C++17
#if __cplusplus > 201402L
#define static_if \
    if            \
    constexpr
#else
#define static_if if
#endif

#if defined(__GNUC__)
#define corbo_DEPRECATED __attribute__((deprecated))
#define corbo_FORCE_INLINE __attribute__((always_inline))
#elif defined(_MSC_VER)
#define corbo_DEPRECATED
#define corbo_FORCE_INLINE __forceinline
#else
#define corbo_DEPRECATED
#define corbo_FORCE_INLINE inline
#endif

}  // namespace corbo

#endif  // SRC_CORE_INCLUDE_CORBO_CORE_MACROS_H_
