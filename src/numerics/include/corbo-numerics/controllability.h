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

#ifndef SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_CONTROLLABILITY_H_
#define SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_CONTROLLABILITY_H_

#include <corbo-core/types.h>

namespace corbo {

/**
 * @brief Methods for checking controllability of dynamic systems
 *
 * @ingroup numerics
 **
 * https://en.wikipedia.org/wiki/Controllability
 * http://www.me.umn.edu/courses/me8281/notes/Chapter%203%20Controllability%20Observability.pdf
 *
 * @see Observability
 *
 * @todo add stabilizability checks (https://de.mathworks.com/help/control/ref/stabsep.html)
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class Controllability
{
 public:
    /**
     * @brief Check controllability of linear time invariant system
     *
     * @param[in]  A       [n x n] matrix
     * @param[in]  B       [n x m] matrix
     * @param[out] rank Retreive rank of the observability matrix [optional]
     * @returns true if system is controllable, false otherwise
     */
    static bool checkLinearTimeInvariantSystem(const Eigen::Ref<const Eigen::MatrixXd>& A, const Eigen::Ref<const Eigen::MatrixXd>& B,
                                               int* rank = nullptr);

 protected:
};

}  // namespace corbo

#endif  // SRC_NUMERICS_INCLUDE_CORBO_NUMERICS_CONTROLLABILITY_H_
