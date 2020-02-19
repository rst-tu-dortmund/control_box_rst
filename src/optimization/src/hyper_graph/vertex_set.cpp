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

#include <corbo-optimization/hyper_graph/vertex_set.h>

namespace corbo {

int VertexSetInterface::getParameterDimension()
{
    int param_dim = 0;

    for (VertexInterface* vertex : getActiveVertices()) param_dim += vertex->getDimensionUnfixed();

    return param_dim;
}

double VertexSetInterface::getParameterValue(int idx)
{
    if (isModified()) computeVertexIndices();

    assert(idx < getParameterDimension());

    // TODO(roesmann) replace linear search strategy by a more efficient one
    for (VertexInterface* vertex : getActiveVertices())
    {
        int vtx_idx = vertex->getVertexIdx();
        if (vertex->getDimensionUnfixed() == vertex->getDimension())
        {
            if (vtx_idx + vertex->getDimension() > idx)
            {
                return vertex->getData()[idx - vtx_idx];
            }
        }
        else
        {
            int param_idx = 0;
            for (int i = 0; i < vertex->getDimension(); ++i)
            {
                if (!vertex->isFixedComponent(i))
                {
                    if (vtx_idx + param_idx == idx) return vertex->getData()[i];
                    ++param_idx;
                }
            }
        }
    }
    return CORBO_MAX_DBL;
}

void VertexSetInterface::setParameterValue(int idx, double x)
{
    if (isModified()) computeVertexIndices();

    assert(idx < getParameterDimension());

    for (VertexInterface* vertex : getActiveVertices())
    {
        int vtx_idx = vertex->getVertexIdx();
        if (vertex->getDimensionUnfixed() == vertex->getDimension())
        {
            if (vtx_idx + vertex->getDimension() > idx)
            {
                vertex->setData(idx - vtx_idx, x);
                return;
            }
        }
        else
        {
            int param_idx = 0;
            for (int i = 0; i < vertex->getDimension(); ++i)
            {
                if (!vertex->isFixedComponent(i))
                {
                    if (vtx_idx + param_idx == idx)
                    {
                        vertex->setData(i, x);
                        return;
                    }
                    ++param_idx;
                }
            }
        }
    }
}

void VertexSetInterface::getParameterVector(Eigen::Ref<Eigen::VectorXd> x)
{
    if (isModified()) computeVertexIndices();

    assert(x.size() == getParameterDimension());

    for (VertexInterface* vertex : getActiveVertices())
    {
        int x_idx = vertex->getVertexIdx();
        if (vertex->hasFixedComponents())
        {
            for (int i = 0; i < vertex->getDimension(); ++i)
            {
                if (vertex->isFixedComponent(i)) continue;

                x[x_idx] = vertex->getData()[i];
                ++x_idx;
            }
        }
        else
        {
            x.segment(x_idx, vertex->getDimension()) = vertex->getDataMap();
        }
    }
}

void VertexSetInterface::setParameterVector(const Eigen::Ref<const Eigen::VectorXd>& x)
{
    if (isModified()) computeVertexIndices();

    assert(x.size() == getParameterDimension());

    for (VertexInterface* vertex : getActiveVertices())
    {
        int x_idx = vertex->getVertexIdx();
        if (vertex->hasFixedComponents())
        {
            for (int i = 0; i < vertex->getDimension(); ++i)
            {
                if (vertex->isFixedComponent(i)) continue;

                vertex->setData(i, x[x_idx]);
                ++x_idx;
            }
        }
        else
        {
            vertex->getDataRawMap() = x.segment(x_idx, vertex->getDimension());
        }
    }
}

void VertexSetInterface::getBounds(Eigen::Ref<Eigen::VectorXd> lb, Eigen::Ref<Eigen::VectorXd> ub)
{
    if (isModified()) computeVertexIndices();

    assert(lb.size() == getParameterDimension());
    assert(ub.size() == getParameterDimension());

    for (const VertexInterface* vertex : getActiveVertices())
    {
        if (vertex->getDimensionUnfixed() == vertex->getDimension())
        {
            lb.segment(vertex->getVertexIdx(), vertex->getDimension()) = vertex->getLowerBoundsMap();
            ub.segment(vertex->getVertexIdx(), vertex->getDimension()) = vertex->getUpperBoundsMap();
        }
        else
        {
            int param_idx = 0;
            for (int i = 0; i < vertex->getDimension(); ++i)
            {
                if (!vertex->isFixedComponent(i))
                {
                    lb(vertex->getVertexIdx() + param_idx) = vertex->getLowerBounds()[i];
                    ub(vertex->getVertexIdx() + param_idx) = vertex->getUpperBounds()[i];
                    ++param_idx;
                }
            }
        }
    }
}

void VertexSetInterface::setBounds(const Eigen::Ref<const Eigen::VectorXd>& lb, const Eigen::Ref<const Eigen::VectorXd>& ub)
{
    if (isModified()) computeVertexIndices();

    assert(lb.size() == getParameterDimension());
    assert(ub.size() == getParameterDimension());

    for (VertexInterface* vertex : getActiveVertices())
    {
        if (vertex->getDimensionUnfixed() == vertex->getDimension())
        {
            vertex->setLowerBounds(lb.segment(vertex->getVertexIdx(), vertex->getDimension()));
            vertex->setUpperBounds(ub.segment(vertex->getVertexIdx(), vertex->getDimension()));
        }
        else
        {
            int param_idx = 0;
            for (int i = 0; i < vertex->getDimension(); ++i)
            {
                if (!vertex->isFixedComponent(i))
                {
                    vertex->setLowerBound(i, lb(vertex->getVertexIdx() + param_idx));
                    vertex->setUpperBound(i, ub(vertex->getVertexIdx() + param_idx));
                    ++param_idx;
                }
            }
        }
    }
}

double VertexSetInterface::getLowerBound(int idx)
{
    if (isModified()) computeVertexIndices();

    assert(idx < getParameterDimension());

    // TODO(roesmann) replace linear search strategy by a more efficient one
    for (VertexInterface* vertex : getActiveVertices())
    {
        int vtx_idx = vertex->getVertexIdx();
        if (vertex->getDimensionUnfixed() == vertex->getDimension())
        {
            if (vtx_idx + vertex->getDimension() > idx)
            {
                return vertex->getLowerBounds()[idx - vtx_idx];
            }
        }
        else
        {
            int param_idx = 0;
            for (int i = 0; i < vertex->getDimension(); ++i)
            {
                if (!vertex->isFixedComponent(i))
                {
                    if (vtx_idx + param_idx == idx) return vertex->getLowerBounds()[i];
                    ++param_idx;
                }
            }
        }
    }
    return CORBO_MAX_DBL;
}

double VertexSetInterface::getUpperBound(int idx)
{
    if (isModified()) computeVertexIndices();

    assert(idx < getParameterDimension());

    // TODO(roesmann) replace linear search strategy by a more efficient one
    for (VertexInterface* vertex : getActiveVertices())
    {
        int vtx_idx = vertex->getVertexIdx();
        if (vertex->getDimensionUnfixed() == vertex->getDimension())
        {
            if (vtx_idx + vertex->getDimension() > idx)
            {
                return vertex->getUpperBounds()[idx - vtx_idx];
            }
        }
        else
        {
            int param_idx = 0;
            for (int i = 0; i < vertex->getDimension(); ++i)
            {
                if (!vertex->isFixedComponent(i))
                {
                    if (vtx_idx + param_idx == idx) return vertex->getUpperBounds()[i];
                    ++param_idx;
                }
            }
        }
    }
    return CORBO_MAX_DBL;
}

void VertexSetInterface::setLowerBound(int idx, double lb)
{
    if (isModified()) computeVertexIndices();

    assert(idx < getParameterDimension());

    for (VertexInterface* vertex : getActiveVertices())
    {
        int vtx_idx = vertex->getVertexIdx();
        if (vertex->getDimensionUnfixed() == vertex->getDimension())
        {
            if (vtx_idx + vertex->getDimension() > idx)
            {
                vertex->setLowerBound(idx - vtx_idx, lb);
                return;
            }
        }
        else
        {
            int param_idx = 0;
            for (int i = 0; i < vertex->getDimension(); ++i)
            {
                if (!vertex->isFixedComponent(i))
                {
                    if (vtx_idx + param_idx == idx)
                    {
                        vertex->setLowerBound(i, lb);
                        return;
                    }
                    ++param_idx;
                }
            }
        }
    }
}

void VertexSetInterface::setUpperBound(int idx, double ub)
{
    if (isModified()) computeVertexIndices();

    assert(idx < getParameterDimension());

    for (VertexInterface* vertex : getActiveVertices())
    {
        int vtx_idx = vertex->getVertexIdx();
        if (vertex->getDimensionUnfixed() == vertex->getDimension())
        {
            if (vtx_idx + vertex->getDimension() > idx)
            {
                vertex->setUpperBound(idx - vtx_idx, ub);
                return;
            }
        }
        else
        {
            int param_idx = 0;
            for (int i = 0; i < vertex->getDimension(); ++i)
            {
                if (!vertex->isFixedComponent(i))
                {
                    if (vtx_idx + param_idx == idx)
                    {
                        vertex->setUpperBound(i, ub);
                        return;
                    }
                    ++param_idx;
                }
            }
        }
    }
}

void VertexSetInterface::applyIncrementNonFixed(const Eigen::Ref<const Eigen::VectorXd>& increment)
{
    if (isModified()) computeVertexIndices();

    assert((int)increment.size() == getParameterDimension());

    for (VertexInterface* vertex : getActiveVertices())
    {
        if (vertex->getDimensionUnfixed() != 0) vertex->plusUnfixed(increment.segment(vertex->getVertexIdx(), vertex->getDimensionUnfixed()).data());
    }
}

void VertexSetInterface::applyIncrementNonFixed(int idx, double increment)
{
    if (isModified()) computeVertexIndices();

    assert(idx < getParameterDimension());

    for (VertexInterface* vertex : getActiveVertices())
    {
        int vtx_idx = vertex->getVertexIdx();
        if (vertex->getDimensionUnfixed() == vertex->getDimension())
        {
            if (vtx_idx + vertex->getDimension() > idx)
            {
                vertex->plus(idx - vtx_idx, increment);
                return;
            }
        }
        else
        {
            int param_idx = 0;
            for (int i = 0; i < vertex->getDimension(); ++i)
            {
                if (!vertex->isFixedComponent(i))
                {
                    if (vtx_idx + param_idx == idx)
                    {
                        vertex->plus(i, increment);
                        return;
                    }
                    ++param_idx;
                }
            }
        }
    }
}

void VertexSetInterface::computeVertexIndices()
{
    if (isModified()) computeActiveVertices();

    std::vector<VertexInterface*>& vertices = getActiveVertices();

    if (vertices.empty()) return;

    setVertexIdx(*vertices.front(), 0);
    for (int i = 1; i < (int)vertices.size(); ++i)
    {
        setVertexIdx(*vertices[i], vertices[i - 1]->getVertexIdx() + vertices[i - 1]->getDimensionUnfixed());
    }
}

void VertexSetInterface::clearConnectedEdges()
{
    if (isModified()) computeActiveVertices();

    std::vector<VertexInterface*>& vertices = getActiveVertices();

    if (vertices.empty()) return;

    for (VertexInterface* vtx : getActiveVertices()) vtx->clearConnectedEdges();
}

void VertexSetInterface::backupParametersActiveVertices()
{
    for (VertexInterface* vertex : getActiveVertices()) vertex->push();
}

void VertexSetInterface::restoreBackupParametersActiveVertices(bool keep_backup)
{
    if (keep_backup)
    {
        for (VertexInterface* vertex : getActiveVertices()) vertex->top();
    }
    else
    {
        for (VertexInterface* vertex : getActiveVertices()) vertex->pop();
    }
}

void VertexSetInterface::discardBackupParametersActiveVertices(bool all)
{
    if (all)
    {
        for (VertexInterface* vertex : getActiveVertices())
        {
            vertex->clearBackups();
        }
    }
    else
    {
        for (VertexInterface* vertex : getActiveVertices())
        {
            vertex->discardTop();
        }
    }
}

// ############# VertexSet ##############

std::vector<VertexInterface*>& VertexSet::getActiveVertices()
{
    if (isModified()) computeActiveVertices();
    return _active_vertices;
}

void VertexSet::getVertices(std::vector<VertexInterface*>& vertices)
{
    vertices.clear();
    for (const VertexInterface::Ptr& vertex : _vertices) vertices.push_back(vertex.get());
}

VertexSet::VertexSet(std::initializer_list<VertexInterface::Ptr> vertices) { _vertices.insert(_vertices.end(), vertices.begin(), vertices.end()); }

void VertexSet::computeActiveVertices()
{
    _active_vertices.clear();
    _active_vertices.reserve(_vertices.size());

    for (VertexInterface::Ptr& vertex : _vertices)
    {
        if (!vertex->isFixed()) _active_vertices.push_back(vertex.get());
    }
}

void VertexSet::clear()
{
    _vertices.clear();
    setModified(true);
}

}  // namespace corbo
