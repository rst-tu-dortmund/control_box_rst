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

#ifndef SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_VERTEX_SET_H_
#define SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_VERTEX_SET_H_

#include <corbo-optimization/hyper_graph/vertex_interface.h>

#include <initializer_list>
#include <memory>
#include <vector>

namespace corbo {

/**
 * @brief Abstract class representing a set of vertices
 *
 * @ingroup hyper-graph
 *
 * @see HyperGraph VertexInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class VertexSetInterface
{
 public:
    using Ptr = std::shared_ptr<VertexSetInterface>;

    VertexSetInterface() {}

    //! Virtual destructor
    virtual ~VertexSetInterface() {}

    int getParameterDimension();

    virtual std::vector<VertexInterface*>& getActiveVertices()        = 0;
    virtual void getVertices(std::vector<VertexInterface*>& vertices) = 0;

    virtual void computeActiveVertices() = 0;
    //! Precompute vertex indices in the hyper-graph (e.g. for the Jacobian or Hessian structure)
    void computeVertexIndices();
    void clearConnectedEdges();

    //! Active vertices related methods
    void applyIncrementNonFixed(const Eigen::Ref<const Eigen::VectorXd>& increment);
    void applyIncrementNonFixed(int idx, double increment);
    double getParameterValue(int idx);
    void setParameterValue(int idx, double x);
    void getParameterVector(Eigen::Ref<Eigen::VectorXd> x);
    void setParameterVector(const Eigen::Ref<const Eigen::VectorXd>& x);
    void getBounds(Eigen::Ref<Eigen::VectorXd> lb, Eigen::Ref<Eigen::VectorXd> ub);
    void setBounds(const Eigen::Ref<const Eigen::VectorXd>& lb, const Eigen::Ref<const Eigen::VectorXd>& ub);
    double getLowerBound(int idx);
    double getUpperBound(int idx);
    void setLowerBound(int idx, double lb);
    void setUpperBound(int idx, double ub);
    void backupParametersActiveVertices();
    void restoreBackupParametersActiveVertices(bool keep_backup);
    void discardBackupParametersActiveVertices(bool all = false);

    void setModified(bool modified) { _modified = modified; }
    bool isModified() const { return _modified; }

    virtual void clear() = 0;

 protected:
    void setVertexIdx(VertexInterface& vertex, int idx) { vertex._vertex_idx = idx; }

    bool _modified = true;
};

class VertexSet : public VertexSetInterface
{
 public:
    using Ptr = std::shared_ptr<VertexSet>;

    VertexSet() = default;
    VertexSet(std::initializer_list<VertexInterface::Ptr> vertices);

    // implements interface method
    std::vector<VertexInterface*>& getActiveVertices() override;
    // implements interface method
    void getVertices(std::vector<VertexInterface*>& vertices) override;
    // implements interface method
    void computeActiveVertices() override;

    // implements interface method
    void clear() override;

    void addVertex(VertexInterface::Ptr vertex)
    {
        setModified(true);
        _vertices.push_back(vertex);
    }

    void addVertices(std::initializer_list<VertexInterface::Ptr> vertices)
    {
        setModified(true);
        _vertices.insert(_vertices.end(), vertices.begin(), vertices.end());
    }

    const std::vector<VertexInterface::Ptr>& getVertices() { return _vertices; }
    std::vector<VertexInterface::Ptr>& getVerticesRef()
    {
        setModified(true);  // Worst case: the user could change the set
        return _vertices;
    }

 private:
    std::vector<VertexInterface::Ptr> _vertices;
    std::vector<VertexInterface*> _active_vertices;
};

}  // namespace corbo

#endif  // SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_VERTEX_SET_H_
