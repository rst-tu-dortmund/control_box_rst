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

#ifndef SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_LYAPUNOV_CONTINUOUS_H_
#define SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_LYAPUNOV_CONTINUOUS_H_

#include <corbo-core/console.h>
#include <corbo-core/types.h>

namespace corbo {

/**
 * @brief Methods for dealing with continuous-time Lyapunov equations
 *
 * @ingroup numerics
 *
 * The continuous-time Lyapunov equation is given by
 * \f[
 *       A X + X A^T + Q = 0
 * \f]
 * If Q is symmetric, solution X is symmetric as well.
 * A and Q must be square.
 *
 * @remarks The Lyapunov equation is a special case of the more general Sylvester equation.
 *                  We still provide this particular class in order to provide a more dedicated implementation
 *                  in terms of efficiency.
 *
 * @see LyapunovDiscrete SylvesterContinuous SylvesterDiscrete AlgebraicRiccatiContinuous AlgebraicRiccatiDiscrete
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 *
 * @todo Add support for the generalized Lyapunov equations \f$ A X E^T + E X A^T + Q = 0 \f$.
 * @todo Allow the user to precompute the Schur decomposition for subsequent calls of solve() with varying \f$ Q \f$.
 */
class LyapunovContinuous
{
 public:
    /**
     * @brief Solve continuous-time Lyapunov equation
     *
     * Solve \f$ A X + X A^T + Q = 0 \f$ w.r.t. \f$ X \f$.
     * If Q is symmetric, solution X is symmetric as well.
     * A and Q must be square.
     * The solution is unique if for Eigenvalues of A it is:
     * \f$ \lambda_i \neq -\lambda_j, \forall i,j. \f$
     *
     * The solution is obtained via Schur decomposition [1-4].
     * We solve the transformed equation: \f$ T Y + Y T^T + F = 0 \f$
     * with
     * \f[
     *   T = U^T A U, \\
     *   F = U^T Q U, \\
     *   Y = U^T X U.
     * \f]
     * In this transformed representation the system of equations can be solved
     * by backward substitution.
     *
     * [1] R.H. Bartels, G. W. Stewart. "Algorithm 432: Solution of the matrix equation AX + XB = C".
     *     Comm. ACM. 15 (9): 820–826, 1972.
     * [2] V. Simoncini. "Computational Methods for Linear Matrix Equations".
     *     SIAM Review 58 (3): 377-441, 2016.
     * [3] https://people.kth.se/~eliasj/NLA/matrixeqs.pdf
     * [4] https://github.com/ajt60gaibb/freeLYAP
     *
     * @warning Input matrix sizes are checked only in debug mode and occuring issues immediately break execution
     *
     * @todo Reimplement using Real Schur form in order to improve efficiency (if suitable)
     *
     * @param[in]  A    Square matrix
     * @param[in]  Q    Square matrix with same size as A
     * @param[out] X    Solution with same size as A and Q (size(X) must not be preallocated a-priori)
     * @returns true if a finite solution is found, false otherwise.
     */
    static bool solve(const Eigen::Ref<const Eigen::MatrixXd>& A, const Eigen::Ref<const Eigen::MatrixXd>& Q, Eigen::MatrixXd& X);

    /**
     * @brief Determine if the Lyapunov equation exhibits a unique solution
     *
     * The solution is unique if for Eigenvalues of A it is:
     * \f$ \lambda_i \neq -\lambda_j, \forall i,j. \f$
     *
     * @param[in]  A    Matrix (might be non-square)
     * @return true if the solution is unique, false otherwise (also if matrix A is not square).
     */
    static bool hasUniqueSolution(const Eigen::Ref<const Eigen::MatrixXd>& A);
};

}  // namespace corbo

#endif  // SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_LYAPUNOV_CONTINUOUS_H_
