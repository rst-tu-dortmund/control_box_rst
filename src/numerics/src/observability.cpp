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

#include <corbo-numerics/observability.h>

#include <corbo-core/console.h>

#include <corbo-numerics/matrix_utilities.h>

// #include <Eigen/LU>
#include <Eigen/QR>

namespace corbo {

bool Observability::checkLinearTimeInvariantSystem(const Eigen::Ref<const Eigen::MatrixXd>& A, const Eigen::Ref<const Eigen::MatrixXd>& C, int* rank)
{
    assert(is_square(A));
    assert(A.rows() == C.cols());
    int n = A.rows();
    int p = C.rows();
    if (n < 1 || p < 1) return false;

    Eigen::MatrixXd observ_mat(n * p, n);

    observ_mat.topRows(p) = C;
    if (n > 1) observ_mat.middleRows(p, p) = C * A;

    if (n > 2)
    {
        Eigen::MatrixXd A_aux = A;
        int current_row       = 2 * p;
        for (int i = 2; i < n; ++i)
        {
            A_aux *= A;
            observ_mat.middleRows(current_row, p) = C * A_aux;
            current_row += p;
        }
    }

    // get rank using FullPivLU or ColPivHousehodler. Note, rank-revealing LU is the most reliable way with floating values
    // see https://en.wikipedia.org/wiki/Rank_%28linear_algebra%29#Computing_the_rank_of_a_matrix

    // Eigen::FullPivLU<Eigen::MatrixXd> decomp(ctrl_mat);
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> decomp(observ_mat);
    int observ_mat_rank = decomp.rank();
    if (rank) *rank     = observ_mat_rank;
    return (observ_mat_rank == n);
}

}  // namespace corbo
