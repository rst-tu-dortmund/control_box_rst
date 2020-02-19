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

#ifndef SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_CORBO_OPTIMAL_CONTROL_H_
#define SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_CORBO_OPTIMAL_CONTROL_H_

#include <corbo-optimal-control/functions/final_state_constraints.h>
#include <corbo-optimal-control/functions/final_state_cost.h>
#include <corbo-optimal-control/functions/hybrid_cost.h>
#include <corbo-optimal-control/functions/minimum_time.h>
#include <corbo-optimal-control/functions/quadratic_control_cost.h>
#include <corbo-optimal-control/functions/quadratic_cost.h>
#include <corbo-optimal-control/functions/quadratic_state_cost.h>
#include <corbo-optimal-control/structured_ocp/discretization_grids/finite_differences_grid.h>
#include <corbo-optimal-control/structured_ocp/discretization_grids/finite_differences_grid_move_blocking.h>
#include <corbo-optimal-control/structured_ocp/discretization_grids/finite_differences_variable_grid.h>
#include <corbo-optimal-control/structured_ocp/discretization_grids/multiple_shooting_grid.h>
#include <corbo-optimal-control/structured_ocp/discretization_grids/multiple_shooting_variable_grid.h>
#include <corbo-optimal-control/structured_ocp/discretization_grids/non_uniform_finite_differences_variable_grid.h>
#include <corbo-optimal-control/structured_ocp/discretization_grids/non_uniform_multiple_shooting_variable_grid.h>
#include <corbo-optimal-control/structured_ocp/structured_optimal_control_problem.h>

#endif  // SRC_OPTIMAL_CONTROL_INCLUDE_CORBO_OPTIMAL_CONTROL_CORBO_OPTIMAL_CONTROL_H_
