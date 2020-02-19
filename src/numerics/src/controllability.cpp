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

#include <corbo-numerics/controllability.h>

#include <corbo-numerics/matrix_utilities.h>

#include <corbo-core/console.h>

// #include <Eigen/LU>
#include <Eigen/QR>

namespace corbo {

bool Controllability::checkLinearTimeInvariantSystem(const Eigen::Ref<const Eigen::MatrixXd>& A, const Eigen::Ref<const Eigen::MatrixXd>& B,
                                                     int* rank)
{
    assert(is_square(A));
    assert(A.rows() == B.rows());
    int n = A.rows();
    int m = B.cols();
    if (n < 1 || m < 1) return false;

    Eigen::MatrixXd ctrl_mat(n, n * m);

    ctrl_mat.leftCols(m) = B;
    if (n > 1) ctrl_mat.middleCols(m, m) = A * B;

    if (n > 2)
    {
        Eigen::MatrixXd A_aux = A;
        int current_col       = 2 * m;
        for (int i = 2; i < n; ++i)
        {
            A_aux *= A;
            ctrl_mat.middleCols(current_col, m) = A_aux * B;
            current_col += m;
        }
    }
    // get rank using FullPivLU or ColPivHousehodler. Note, rank-revealing LU is the most reliable way with floating values
    // see https://en.wikipedia.org/wiki/Rank_%28linear_algebra%29#Computing_the_rank_of_a_matrix

    // Eigen::FullPivLU<Eigen::MatrixXd> decomp(ctrl_mat);
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> decomp(ctrl_mat);
    int ctrl_mat_rank = decomp.rank();
    if (rank) *rank   = ctrl_mat_rank;
    return (ctrl_mat_rank == n);
}

}  // namespace corbo
