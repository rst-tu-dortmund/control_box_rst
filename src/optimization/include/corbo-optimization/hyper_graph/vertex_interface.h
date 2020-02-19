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

#ifndef SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_VERTEX_INTERFACE_H_
#define SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_VERTEX_INTERFACE_H_

#include <corbo-core/types.h>

#include <memory>
#include <set>

namespace corbo {

class EdgeInterface;  // Forward definition, include in child headers
class BaseEdge;       // Forward definition, include in child headers
class BaseMixedEdge;  // Forward definition, include in child headers

/**
 * @brief Generic interface class for vertices.
 *
 * @ingroup optimization hyper-graph
 *
 * This abstract class defines the interface for vertices.
 * The optimization problem might be formulated as a hyper-graph
 * in which optimization parameters are represented as vertices.
 * Cost functions and constraints are represented as edges (refer to EdgeInterface).
 *
 * @see HyperGraph EdgeInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class VertexInterface
{
    friend class VertexSetInterface;

 public:
    using Ptr  = std::shared_ptr<VertexInterface>;
    using UPtr = std::unique_ptr<VertexInterface>;

    //! Virtual destructor
    virtual ~VertexInterface() {}

    //! Return number of elements/values/components stored in this vertex.
    virtual int getDimension() const = 0;
    //! Return number of unfixed elements (unfixed elements are skipped as parameters in the Hessian and Jacobian.
    virtual int getDimensionUnfixed() const = 0;
    //! Check if all components are fixed
    virtual bool isFixed() const { return getDimensionUnfixed() == 0; }

    //! Register an adjacent objective edge
    void registerObjectiveEdge(BaseEdge* edge);
    //! Register an adjacent least-squares objective edge
    void registerLsqObjectiveEdge(BaseEdge* edge);
    //! Register an adjacent equality constraint edge
    void registerEqualityEdge(BaseEdge* edge);
    //! Register an adjacent inequality constraint edge
    void registerInequalityEdge(BaseEdge* edge);
    //! Register an adjacent mixed edge
    void registerMixedEdge(BaseMixedEdge* edge);
    //! Clear all connected edges
    void clearConnectedEdges()
    {
        _edges_objective.clear();
        _edges_lsq_objective.clear();
        _edges_equalities.clear();
        _edges_inequalities.clear();
        _edges_mixed.clear();
    }

    //! Add value to a specific component of the vertex: x[idx] += inc
    virtual void plus(int idx, double inc) = 0;
    //! Define the increment for the vertex: x += inc with dim(inc)=getDimension().
    virtual void plus(const double* inc) = 0;
    //! Define the increment for the unfixed components of the vertex: x += inc with dim(inc)=getDimensionUnfixed().
    virtual void plusUnfixed(const double* inc) = 0;

    //! Get read-only raw access to the values of the vertex
    virtual const double* getData() const = 0;
    //! Get write access to the values of the vertex
    virtual double* getDataRaw() = 0;
    //! Get a read-only Eigen::Map to the values of the vertex
    Eigen::Map<const Eigen::VectorXd> getDataMap() const { return Eigen::Map<const Eigen::VectorXd>(getData(), getDimension()); }
    //! Get a Eigen::Map to the values of the vertex
    Eigen::Map<Eigen::VectorXd> getDataRawMap() { return Eigen::Map<Eigen::VectorXd>(getDataRaw(), getDimension()); }

    //! Write data to to a specific component
    virtual void setData(int idx, double data) = 0;

    //! Check if the vertex has fixed components
    virtual bool hasFixedComponents() const = 0;

    //! Check if individual components are fixed or unfixed
    virtual bool isFixedComponent(int idx) const = 0;

    //! Define lower bounds on the vertex values [getDimension() x 1]
    virtual void setLowerBounds(const Eigen::Ref<const Eigen::VectorXd>& lb) = 0;
    //! Define upper bounds on the vertex values [getDimension() x 1]
    virtual void setUpperBounds(const Eigen::Ref<const Eigen::VectorXd>& ub) = 0;
    //! Define lower bound on a single component of the vertex (0 <= idx < getDimension())
    virtual void setLowerBound(int idx, double lb) = 0;
    //! Define upper bound on a single component of the vertex (0 <= idx < getDimension())
    virtual void setUpperBound(int idx, double ub) = 0;
    //! Check if finite bounds (lower or upper) are provided
    virtual bool hasFiniteBounds() const = 0;
    //! Check if finite lower bounds are provided
    virtual bool hasFiniteLowerBounds() const = 0;
    //! Check if finite upper bounds are provided
    virtual bool hasFiniteUpperBounds() const = 0;
    //! Check if finite lower bound for a single component is provided
    virtual bool hasFiniteLowerBound(int idx) const = 0;
    //! Check if finite upper bound for a single component is provided
    virtual bool hasFiniteUpperBound(int idx) const = 0;
    //! Get number of finite lower bounds
    virtual int getNumberFiniteLowerBounds(bool unfixed_only) const = 0;
    //! Get number of finite upper bounds
    virtual int getNumberFiniteUpperBounds(bool unfixed_only) const = 0;
    //! Get number of finite upper bounds (either upper or lower must be finite)
    virtual int getNumberFiniteBounds(bool unfixed_only) const = 0;

    //! Read-only raw access to lower bounds [getDimension() x 1]
    virtual const double* getLowerBounds() const = 0;
    //! Read-only Eigen::Map for lower bounds [getDimension() x 1]
    Eigen::Map<const Eigen::VectorXd> getLowerBoundsMap() const { return Eigen::Map<const Eigen::VectorXd>(getLowerBounds(), getDimension()); }
    //! Read-only raw access to upper bounds [getDimension() x 1]
    virtual const double* getUpperBounds() const = 0;
    //! Read-only Eigen::Map for upper bounds [getDimension() x 1]
    Eigen::Map<const Eigen::VectorXd> getUpperBoundsMap() const { return Eigen::Map<const Eigen::VectorXd>(getUpperBounds(), getDimension()); }

    //! Raw access for connected objective edges
    const std::set<BaseEdge*>& getConnectedObjectiveEdgesRef() const { return _edges_objective; }
    //! Raw access for connected least-squares objective edges
    const std::set<BaseEdge*>& getConnectedLsqObjectiveEdgesRef() const { return _edges_lsq_objective; }
    //! Raw access for connected equality constraint edges
    const std::set<BaseEdge*>& getConnectedEqualityEdgesRef() const { return _edges_equalities; }
    //! Raw access for connected inequality constraint edges
    const std::set<BaseEdge*>& getConnectedInequalityEdgesRef() const { return _edges_inequalities; }
    //! Raw access for connected mixed edges
    const std::set<BaseMixedEdge*>& getConnectedMixedEdgesRef() const { return _edges_mixed; }

    // Backup values
    virtual void push()       = 0;  //!< Store all values into a internal backup stack.
    virtual void pop()        = 0;  //!< Restore the previously stored values of the backup stack and remove them from the stack.
    virtual void top()        = 0;  //!< Restore the previously stored values of the backup stack WITHOUT removing them from the stack.
    virtual void discardTop() = 0;  //!< Delete the previously made backup from the stack without restoring it.
    virtual void clear()      = 0;  //!< Clear complete backup container.
    virtual void clearBackups() { clear(); }
    // virtual void copyBackup(int k, double* values) = 0;   //!< Make sure that values complies with dimension()
    virtual int getNumBackups() const = 0;  //!< Return the current size/number of backups of the backup stack.

    //! Retrieve current edge index (warning, this value is determined within the related VertexSetInterface implementation)
    int getVertexIdx() const { return _vertex_idx; }
    //! Retrieve number of connected objective edges with custom Jacobian
    // int getNumObjectiveEdgesWithCustomJacobian() const;
    //! Retrieve number of connected equality constraint edges with custom Jacobian
    // int getNumEqualityEdgesWithCustomJacobian() const;
    //! Retrieve number of connected inequality constraint edges with custom Jacobian
    // int getNumInequalityEdgesWithCustomJacobian() const;

 private:
    std::set<BaseEdge*> _edges_objective;      //!< connected objective edges
    std::set<BaseEdge*> _edges_lsq_objective;  //!< connected least-squares objective edges
    std::set<BaseEdge*> _edges_equalities;     //!< connected equality constraint edges
    std::set<BaseEdge*> _edges_inequalities;   //!< connected inequality constraint edges
    std::set<BaseMixedEdge*> _edges_mixed;     //!< connected mixed edges

    int _vertex_idx = 0;  //!< vertex index in jacobian or hessian (determined by friend class HyperGraph).
    // int _num_edges_with_custom_jacob = 0;  //!< we can skip jacobian computation for this vertex if all edges provide their own jacobian
};

}  // namespace corbo

#endif  // SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_VERTEX_INTERFACE_H_
