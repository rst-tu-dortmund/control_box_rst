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

#ifndef SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_MISC_H_
#define SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_MISC_H_

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace corbo {

void convert_triplet(const Eigen::Ref<const Eigen::VectorXi>& i_row, const Eigen::Ref<const Eigen::VectorXi>& j_col,
                     const Eigen::Ref<const Eigen::VectorXd>& values, std::vector<Eigen::Triplet<double>>& triplets);

}  // namespace corbo

#endif  // SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_MISC_H_
