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

#include <corbo-core/value_comparison.h>
#include <corbo-numerics/matrix_utilities.h>
#include <corbo-numerics/schur.h>

#include <utility>
#include <vector>

namespace corbo {

template <class Predicate>
bool reorder_schur_blocks(Eigen::Ref<Eigen::MatrixXd> T, Eigen::Ref<Eigen::MatrixXd> Q, Predicate predicate, int* subspace_dim, bool standardize)
{
    assert(is_square(T));
    assert(have_equal_size(T, Q));

    if (T.rows() < 2) return true;

    std::vector<std::pair<int, int>> blocks;  // row_idx, block size
    blocks.reserve(T.rows());
    std::vector<int> selected;
    selected.reserve(T.rows());

    int subspace_dim_internal = 0;

    // search for blocks and determine if they are selected (predicate == true)
    int row = 0;
    while (row < T.rows())
    {
        // check, if we have a 1 x 1 or 2 x 2 block
        int p;
        if (row + 1 >= T.rows())
            p = 1;
        else if (approx_zero(T(row + 1, row), 1e-10))
            p = 1;
        else
            p = 2;

        if (predicate(T.block(row, row, p, p)))
        {
            selected.push_back((int)blocks.size());  // add current block index
            subspace_dim_internal += p;
        }

        blocks.emplace_back(row, p);

        row += p;
    }

    if (selected.empty())
    {
        if (subspace_dim) *subspace_dim = 0;
        return true;
    }

    int top = 0;
    for (int l = 0; l < (int)selected.size(); ++l)
    {
        int current_idx = selected[l];  // current_idx is is increasing only with l
        for (int i = current_idx; i >= top + 1; --i)
        {
            // swap A_{i-1} with A_{i}
            if (!swap_schur_blocks(T, blocks[i - 1].first, blocks[i - 1].second, blocks[i].second, Q, standardize)) return false;
            // update index
            int p_prev           = blocks[i - 1].second;
            blocks[i - 1].second = blocks[i].second;
            blocks[i].first      = blocks[i - 1].first + blocks[i].second;
            blocks[i].second     = p_prev;
        }
        ++top;
    }

    if (subspace_dim) *subspace_dim = subspace_dim_internal;
    return true;
}

}  // namespace corbo
