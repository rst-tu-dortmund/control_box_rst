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

#ifndef SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_HYPER_GRAPH_OPTIMIZATION_PROBLEM_IGNORE_STRUCTURE_H_
#define SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_HYPER_GRAPH_OPTIMIZATION_PROBLEM_IGNORE_STRUCTURE_H_

#include <corbo-optimization/hyper_graph/hyper_graph_optimization_problem_base.h>

#include <memory>

namespace corbo {

class HyperGraphOptimizationProblemIgnoreStructure : public BaseHyperGraphOptimizationProblem
{
 public:
    HyperGraphOptimizationProblemIgnoreStructure() { _warn_if_not_specialized = false; }
    HyperGraphOptimizationProblemIgnoreStructure(OptimizationEdgeSet::Ptr edges, VertexSetInterface::Ptr vertices)
        : BaseHyperGraphOptimizationProblem(edges, vertices)
    {
        _warn_if_not_specialized = false;  // we want to use default implementations even though there are slow
    }

    BaseHyperGraphOptimizationProblem::Ptr getInstance() const override { return std::make_shared<HyperGraphOptimizationProblemIgnoreStructure>(); }
};

FACTORY_REGISTER_HYPER_GRAPH_OPTIMIZATION_PROBLEM(HyperGraphOptimizationProblemIgnoreStructure)

}  // namespace corbo

#endif  // SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_HYPER_GRAPH_OPTIMIZATION_PROBLEM_IGNORE_STRUCTURE_H_
