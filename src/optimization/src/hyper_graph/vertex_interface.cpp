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

#include <corbo-optimization/hyper_graph/vertex_interface.h>

#include <corbo-optimization/hyper_graph/edge.h>

namespace corbo {

void VertexInterface::registerObjectiveEdge(BaseEdge* edge) { _edges_objective.insert(edge); }
void VertexInterface::registerLsqObjectiveEdge(BaseEdge* edge) { _edges_lsq_objective.insert(edge); }
void VertexInterface::registerEqualityEdge(BaseEdge* edge) { _edges_equalities.insert(edge); }
void VertexInterface::registerInequalityEdge(BaseEdge* edge) { _edges_inequalities.insert(edge); }
void VertexInterface::registerMixedEdge(BaseMixedEdge* edge) { _edges_mixed.insert(edge); }

/*
int VertexInterface::getNumObjectiveEdgesWithCustomJacobian() const
{
    int num = 0;
    for (EdgeInterface* edge : _edges_objective) num += edge->providesJacobian();
    return num;
}

int VertexInterface::getNumEqualityEdgesWithCustomJacobian() const
{
    int num = 0;
    for (EdgeInterface* edge : _edges_equalities) num += edge->providesJacobian();
    return num;
}

int VertexInterface::getNumInequalityEdgesWithCustomJacobian() const
{
    int num = 0;
    for (EdgeInterface* edge : _edges_inequalities) num += edge->providesJacobian();
    return num;
}
*/

}  // namespace corbo
