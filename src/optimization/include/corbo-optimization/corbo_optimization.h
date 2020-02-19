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

#ifndef SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_CORBO_OPTIMIZATION_H_
#define SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_CORBO_OPTIMIZATION_H_

#include <corbo-optimization/hyper_graph/hyper_graph_optimization_problem_edge_based.h>
#include <corbo-optimization/hyper_graph/hyper_graph_optimization_problem_ignore_structure.h>
#include <corbo-optimization/hyper_graph/hyper_graph_optimization_problem_vertex_based.h>
#include <corbo-optimization/simple_optimization_problem.h>

#include <corbo-optimization/solver/levenberg_marquardt_dense.h>
#include <corbo-optimization/solver/levenberg_marquardt_sparse.h>
#include <corbo-optimization/solver/nlp_solver_ipopt.h>
#include <corbo-optimization/solver/qp_solver_osqp.h>

#endif  // SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_CORBO_OPTIMIZATION_H_
