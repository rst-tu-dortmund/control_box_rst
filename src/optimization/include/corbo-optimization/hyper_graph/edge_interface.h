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

#ifndef SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_EDGE_INTERFACE_H_
#define SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_EDGE_INTERFACE_H_

#include <corbo-core/types.h>
#include <corbo-optimization/hyper_graph/edge_cache.h>

#include <memory>

namespace corbo {

class VertexInterface;     // Forward declaration, include in subclass!
class VertexSetInterface;  // Forward declaration, include in subclass!
class EdgeSetInterface;    // Forward declaration for friend access

/**
 * @brief Generic interface class for edges.
 *
 * @ingroup optimization hyper-graph
 *
 * This abstract class defines the interface for edges for the HyperGraph.
 * Cost functions and constraints (that depend on vertices) are formulated as multi-edges.
 * Multi-edges can connect an arbitary number of vertices. The resulting graph
 * is called hyper-graph. If the number of vertices attached to a single
 * edge is low, the resulting optimization problem is likely to be sparse.
 *
 * @see HyperGraph, VertexInterface, BaseEdge
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class EdgeInterface
{
    friend class EdgeSetInterface;

 public:
    using Ptr  = std::shared_ptr<EdgeInterface>;
    using UPtr = std::unique_ptr<EdgeInterface>;

    //! Virtual destructor
    virtual ~EdgeInterface() {}

    //! Get dimension of the edge (dimension of the cost-function/constraint value vector)
    virtual int getDimension() const = 0;

    /**
     * @brief Compute function values
     *
     * Here, the actual cost/constraint function values are computed:
     * - objective in non-least-squares form: e(x) (hereby, the actual cost is f(x) = e(x)^T e(x))
     * - objective in least-squares form: f(x)
     * - equality constraints: ceq(x)  (in case constraints are satisfied: ceq(x) = 0)
     * - inequality constraints: c(x)  (in case constraints are satisfied: c(x) < 0)
     * @param[in] values     values should be stored here according to getDimension().
     */
    virtual void computeValues(Eigen::Ref<Eigen::VectorXd> values) = 0;

    virtual double computeSumOfValues()
    {
        Eigen::VectorXd values(getDimension());
        computeValues(values);
        return values.sum();
    }

    virtual double computeSquaredNormOfValues()
    {
        Eigen::VectorXd values(getDimension());
        computeValues(values);
        return values.squaredNorm();
    }

    //! Return number of attached vertices
    virtual int getNumVertices() const = 0;
    //! Return the combined dimension of all attached vertices (excluding fixed vertex components)
    virtual int verticesDimension() const = 0;

    //! Get access to vertex with index \c idx (0 <= idx < numVertices)
    virtual VertexInterface* getVertexRaw(int idx)          = 0;
    virtual const VertexInterface* getVertex(int idx) const = 0;

    int getNumFiniteVerticesLowerBounds() const;
    int getNumFiniteVerticesUpperBounds() const;
};

class BaseEdge : public EdgeInterface
{
    friend class EdgeSetInterface;

 public:
    using Ptr  = std::shared_ptr<BaseEdge>;
    using UPtr = std::unique_ptr<BaseEdge>;

    //! Get dimension of the edge (dimension of the cost-function/constraint value vector)
    int getDimension() const override = 0;

    /**
     * @brief Defines if the edge is formulated as Least-Squares form
     *
     * Least-squares cost terms are defined as \f$ f(x) = e(x)^T e(x) \f$ and the
     * function values and Jacobian are computed for \f$ e(x) \f$
     * rather than for \f$ f(x) \f$. Specialiezed least-squares solvers
     * require the optimization problem to be defined in this particular form.
     * Other solvers can automatically compute the square of least-squares edges
     * if required. However, the other way round is more restrictive:
     * general solvers might not cope with non-least-squares forms.
     *
     * Note, in the LS-form case computeValues() computes e(x) and
     * otherwise f(x).
     *
     * @returns true if the edge is given in LS-form
     */
    virtual bool isLeastSquaresForm() const { return false; }

    //! Return true if the edge is linear (and hence its Hessian is always zero)
    virtual bool isLinear() const { return false; }
    //! Return true if a custom Jacobian is provided (e.g. computeJacobian() is overwritten)
    virtual bool providesJacobian() const { return false; }
    //! Return true if a custom Hessian is provided (e.g. computeHessian() is overwritten)
    virtual bool providesHessian() const { return false; }

    void computeValues(Eigen::Ref<Eigen::VectorXd> values) override = 0;

    /**
     * @brief Compute edge block jacobian for a given vertex
     *
     * This interface class provides a numerical computation of the Jacobian using finite differences.
     * However, a user-defined implementation might be specified by overwriting this method.
     * In that case do not forget to also overwrite providesJacobian().
     *
     * @param[in]   vtx_idx          Vertex number for which the block jacobian should be computed
     * @param[out]  block_jacobian   Resulting block jacobian [dimension() x VertexInterface::dimensionUnfixed()] (must be preallocated)
     * @param[in]   multipliers      Scale each value component by an optional multiplier if argument is not null
     * null.
     */
    virtual void computeJacobian(int vtx_idx, Eigen::Ref<Eigen::MatrixXd> block_jacobian, const double* multipliers = nullptr);

    /**
     * @brief Compute edge block Hessian for a given vertex pair
     *
     * This interface class provides a numerical computation of the Hessian using finite differences.
     * However, a user-defined implementation might be specified by overwriting this method.
     * In that case do not forget to also overwrite providesHessian().
     *
     * Given two vertices i and j, four matrices can be computed:
     * Hii, Hij, Hji, Hjj. However, usually it is Hij = Hji.
     *
     * @param[in]   vtx_idx_i           First vertex number
     * @param[in]   vtx_idx_j           Second vertex number
     * @param[in]   block_jacobian_i    Block jacobian for vertex i which might be taken into account:
     *                                  [dimension() x VertexInterface::dimensionUnfixed()] (must be preallocated)
     * @param[out]  block_hessian_ij    Resulting block hessian for vertex pair i and j:
     *                                  [VertexInterface::dimensionUnfixed() x VertexInterface::dimensionUnfixed()]
     * @param[in]   multipliers         Scale each value component by an optional multiplier if argument is not null
     */
    virtual void computeHessianInc(int vtx_idx_i, int vtx_idx_j, const Eigen::Ref<const Eigen::MatrixXd>& block_jacobian_i,
                                   Eigen::Ref<Eigen::MatrixXd> block_hessian_ij, const double* multipliers = nullptr, double weight = 1.0);

    virtual void computeHessian(int vtx_idx_i, int vtx_idx_j, const Eigen::Ref<const Eigen::MatrixXd>& block_jacobian_i,
                                Eigen::Ref<Eigen::MatrixXd> block_hessian_ij, const double* multipliers = nullptr, double weight = 1.0);

    virtual void computeHessianInc(int vtx_idx_i, int vtx_idx_j, Eigen::Ref<Eigen::MatrixXd> block_hessian_ij, const double* multipliers = nullptr,
                                   double weight = 1.0);

    /** @name Caching of edge values and block  jacobians */
    //@{

    void reserveCacheMemory(int num_value_vectors, int num_jacobians)
    {
        reserveValuesCacheMemory(num_value_vectors);
        reserveJacobiansCacheMemory(num_jacobians);
    }
    void reserveValuesCacheMemory(int num_value_vectors) { _cache.reserveMemoryValues(num_value_vectors); }
    void reserveJacobiansCacheMemory(int num_jacobians) { _cache.reserveMemoryJacobians(num_jacobians); }

    //! Call computeValues() and store result to previously allocated internal cache (call allocateInteralValuesCache() first once)
    void computeValuesCached() { computeValues(_cache.pushValues(getDimension())); }

    //! compute the specialied squared-norm method for computing the values (note only the first element in the values cache is used)
    void computeSquaredNormOfValuesCached() { _cache.pushValues(1)[0] = computeSquaredNormOfValues(); }

    //! Retreive values computed previously via computeValuesCached()
    EdgeCache& getCache() { return _cache; }
    const EdgeCache& getCache() const { return _cache; }

    //@}

    //! Retrieve current edge index (warning, this value is determined within the related HyperGraph)
    int getEdgeIdx() const { return _edge_idx; }

 protected:
    int _edge_idx = 0;
    EdgeCache _cache;
};

class BaseMixedEdge : public EdgeInterface
{
    friend class EdgeSetInterface;

 public:
    using Ptr  = std::shared_ptr<BaseMixedEdge>;
    using UPtr = std::unique_ptr<BaseMixedEdge>;

    int getDimension() const override { return getObjectiveDimension() + getEqualityDimension() + getInequalityDimension(); }

    virtual int getObjectiveDimension() const  = 0;
    virtual int getEqualityDimension() const   = 0;
    virtual int getInequalityDimension() const = 0;

    virtual bool isObjectiveLeastSquaresForm() const = 0;

    //! Return true if the edge is linear (and hence its Hessian is always zero)
    virtual bool isObjectiveLinear() const { return false; }
    virtual bool isEqualityLinear() const { return false; }
    virtual bool isInequalityLinear() const { return false; }

    virtual void precompute()                                                     = 0;
    virtual void computeObjectiveValues(Eigen::Ref<Eigen::VectorXd> obj_values)   = 0;
    virtual void computeEqualityValues(Eigen::Ref<Eigen::VectorXd> eq_values)     = 0;
    virtual void computeInequalityValues(Eigen::Ref<Eigen::VectorXd> ineq_values) = 0;

    void computeValues(Eigen::Ref<Eigen::VectorXd> values) final
    {
        precompute();
        computeObjectiveValues(values);
        computeEqualityValues(values.segment(getObjectiveDimension(), getEqualityDimension()));
        computeInequalityValues(values.tail(getInequalityDimension()));
    }

    virtual double computeSumOfObjectiveValues()
    {
        Eigen::VectorXd values(getObjectiveDimension());
        computeObjectiveValues(values);
        return values.sum();
    }

    virtual double computeSquaredNormOfObjectiveValues()
    {
        Eigen::VectorXd values(getObjectiveDimension());
        computeObjectiveValues(values);
        return values.squaredNorm();
    }

    virtual void computeObjectiveJacobian(int vtx_idx, Eigen::Ref<Eigen::MatrixXd> block_jacobian, const double* multipliers = nullptr);
    virtual void computeEqualityJacobian(int vtx_idx, Eigen::Ref<Eigen::MatrixXd> block_jacobian, const double* multipliers = nullptr);
    virtual void computeInequalityJacobian(int vtx_idx, Eigen::Ref<Eigen::MatrixXd> block_jacobian, const double* multipliers = nullptr);
    virtual void computeConstraintJacobians(int vtx_idx, Eigen::Ref<Eigen::MatrixXd> eq_jacobian, Eigen::Ref<Eigen::MatrixXd> ineq_jacobian,
                                            const double* eq_multipliers = nullptr, const double* ineq_multipliers = nullptr);
    virtual void computeJacobians(int vtx_idx, Eigen::Ref<Eigen::MatrixXd> obj_jacobian, Eigen::Ref<Eigen::MatrixXd> eq_jacobian,
                                  Eigen::Ref<Eigen::MatrixXd> ineq_jacobian, const double* obj_multipliers = nullptr,
                                  const double* eq_multipliers = nullptr, const double* ineq_multipliers = nullptr);

    virtual void computeObjectiveHessian(int vtx_idx_i, int vtx_idx_j, const Eigen::Ref<const Eigen::MatrixXd>& block_jacobian_i,
                                         Eigen::Ref<Eigen::MatrixXd> block_hessian_ij, const double* multipliers = nullptr, double weight = 1.0);
    virtual void computeEqualityHessian(int vtx_idx_i, int vtx_idx_j, const Eigen::Ref<const Eigen::MatrixXd>& block_jacobian_i,
                                        Eigen::Ref<Eigen::MatrixXd> block_hessian_ij, const double* multipliers = nullptr, double weight = 1.0);
    virtual void computeInequalityHessian(int vtx_idx_i, int vtx_idx_j, const Eigen::Ref<const Eigen::MatrixXd>& block_jacobian_i,
                                          Eigen::Ref<Eigen::MatrixXd> block_hessian_ij, const double* multipliers = nullptr, double weight = 1.0);
    virtual void computeHessians(int vtx_idx_i, int vtx_idx_j, const Eigen::Ref<const Eigen::MatrixXd>& obj_jacobian_i,
                                 const Eigen::Ref<const Eigen::MatrixXd>& eq_jacobian_i, const Eigen::Ref<const Eigen::MatrixXd>& ineq_jacobian_i,
                                 Eigen::Ref<Eigen::MatrixXd> obj_hessian_ij, Eigen::Ref<Eigen::MatrixXd> eq_hessian_ij,
                                 Eigen::Ref<Eigen::MatrixXd> ineq_hessian_ij, const double* multipliers_eq = nullptr,
                                 const double* multipliers_ineq = nullptr, double weight_obj = 1.0);
    virtual void computeConstraintHessians(int vtx_idx_i, int vtx_idx_j, const Eigen::Ref<const Eigen::MatrixXd>& eq_jacobian_i,
                                           const Eigen::Ref<const Eigen::MatrixXd>& ineq_jacobian_i, Eigen::Ref<Eigen::MatrixXd> eq_hessian_ij,
                                           Eigen::Ref<Eigen::MatrixXd> ineq_hessian_ij, const double* multipliers_eq = nullptr,
                                           const double* multipliers_ineq = nullptr);

    virtual void computeObjectiveHessianInc(int vtx_idx_i, int vtx_idx_j, const Eigen::Ref<const Eigen::MatrixXd>& block_jacobian_i,
                                            Eigen::Ref<Eigen::MatrixXd> block_hessian_ij, const double* multipliers = nullptr, double weight = 1.0);
    virtual void computeEqualityHessianInc(int vtx_idx_i, int vtx_idx_j, const Eigen::Ref<const Eigen::MatrixXd>& block_jacobian_i,
                                           Eigen::Ref<Eigen::MatrixXd> block_hessian_ij, const double* multipliers = nullptr, double weight = 1.0);
    virtual void computeInequalityHessianInc(int vtx_idx_i, int vtx_idx_j, const Eigen::Ref<const Eigen::MatrixXd>& block_jacobian_i,
                                             Eigen::Ref<Eigen::MatrixXd> block_hessian_ij, const double* multipliers = nullptr, double weight = 1.0);
    virtual void computeHessiansInc(int vtx_idx_i, int vtx_idx_j, const Eigen::Ref<const Eigen::MatrixXd>& obj_jacobian_i,
                                    const Eigen::Ref<const Eigen::MatrixXd>& eq_jacobian_i, const Eigen::Ref<const Eigen::MatrixXd>& ineq_jacobian_i,
                                    Eigen::Ref<Eigen::MatrixXd> obj_hessian_ij, Eigen::Ref<Eigen::MatrixXd> eq_hessian_ij,
                                    Eigen::Ref<Eigen::MatrixXd> ineq_hessian_ij, const double* multipliers_eq = nullptr,
                                    const double* multipliers_ineq = nullptr, double weight_obj = 1.0);
    virtual void computeConstraintHessiansInc(int vtx_idx_i, int vtx_idx_j, const Eigen::Ref<const Eigen::MatrixXd>& eq_jacobian_i,
                                              const Eigen::Ref<const Eigen::MatrixXd>& ineq_jacobian_i, Eigen::Ref<Eigen::MatrixXd> eq_hessian_ij,
                                              Eigen::Ref<Eigen::MatrixXd> ineq_hessian_ij, const double* multipliers_eq = nullptr,
                                              const double* multipliers_ineq = nullptr);

    /** @name Caching of edge values and block  jacobians */
    //@{

    void reserveCacheMemory(int num_value_vectors, int num_jacobians)
    {
        reserveValuesCacheMemory(num_value_vectors, num_value_vectors, num_value_vectors);
        reserveJacobiansCacheMemory(num_jacobians, num_jacobians, num_jacobians);
    }
    void reserveValuesCacheMemory(int num_obj_values, int num_eq_values, int num_ineq_values)
    {
        _objective_cache.reserveMemoryValues(num_obj_values);
        _equality_cache.reserveMemoryValues(num_eq_values);
        _inequality_cache.reserveMemoryValues(num_ineq_values);
    }
    void reserveJacobiansCacheMemory(int num_obj_jacobians, int num_eq_jacobians, int num_ineq_jacobians)
    {
        _objective_cache.reserveMemoryJacobians(num_obj_jacobians);
        _equality_cache.reserveMemoryJacobians(num_eq_jacobians);
        _inequality_cache.reserveMemoryJacobians(num_ineq_jacobians);
    }

    //! Call computeObjectiveValues() and store result to the internal cache
    void computeObjectiveValuesCached() { computeObjectiveValues(_objective_cache.pushValues(getObjectiveDimension())); }
    //! compute the specialied squared-norm method for computing the values (note only the first element in the values cache is used)
    void computeSquaredNormOfObjectiveValuesCached() { _objective_cache.pushValues(1)[0] = computeSquaredNormOfObjectiveValues(); }
    //! Call computeEqualityValues() and store result to the internal cache
    void computeEqualityValuesCached() { computeEqualityValues(_equality_cache.pushValues(getEqualityDimension())); }
    //! Call computeInequalityValues() and store result to the internal cache
    void computeInequalityValuesCached() { computeInequalityValues(_inequality_cache.pushValues(getInequalityDimension())); }

    EdgeCache& getObjectiveCache() { return _objective_cache; }
    const EdgeCache& getObjectiveCache() const { return _objective_cache; }

    EdgeCache& getEqualityCache() { return _equality_cache; }
    const EdgeCache& getEqualityCache() const { return _equality_cache; }

    EdgeCache& getInequalityCache() { return _inequality_cache; }
    const EdgeCache& getInequalityCache() const { return _inequality_cache; }

    //@}

    //! Retrieve current edge index (warning, this value is determined within the related HyperGraph)
    int getEdgeObjectiveIdx() const { return _edge_idx_obj; }
    int getEdgeEqualityIdx() const { return _edge_idx_eq; }
    int getEdgeInequalityIdx() const { return _edge_idx_ineq; }

 protected:
    int _edge_idx_obj  = 0;
    int _edge_idx_eq   = 0;
    int _edge_idx_ineq = 0;

    EdgeCache _objective_cache;
    EdgeCache _equality_cache;
    EdgeCache _inequality_cache;
};

}  // namespace corbo

#endif  // SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_EDGE_INTERFACE_H_
