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

#ifndef SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_HYPER_GRAPH_OPTIMIZATION_PROBLEM_BASE_H_
#define SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_HYPER_GRAPH_OPTIMIZATION_PROBLEM_BASE_H_

#include <corbo-optimization/optimization_problem_interface.h>

#include <corbo-optimization/hyper_graph/hyper_graph.h>

#include <corbo-core/factory.h>

#include <functional>
#include <memory>

#if MESSAGE_SUPPORT
#include <corbo-communication/messages/optimization/hyper_graph_optimization_problem.pb.h>
#endif

namespace corbo {

/**
 * @brief Hyper-graph optimization problem formulation
 *
 * @ingroup optimization hyper-graph
 *
 * For the mathematical description refer to the OptimizationProblemInterface.
 *
 * Edges in a hyper-graph are cost terms. Vertices optimization parameters
 * (e.g. states or control inputs). The complete cost-function is defined
 * in terms of the sum of all edge cost-terms.
 * The major advantage of a hyper-graph formulation is that
 * the graph is modular and directly mimics the sparsity structure of the optimization
 * problem. Hence, by iterating vertices and edges only the non-zero portions of the
 * values/Jacobian/Hessian are computed.
 *
 * Note, the graph can be operated in Least-Squares mode (all least-squares
 * edges are retained as cost subvectors: f = e) or standard mode (Least-squares edges are
 * multiplied with its transpose():  f = e^t * e, hereby f is scalar).
 * The selected mode influences the objective dimension. Note,
 * many Least-Squares solvers (e.g. LevenbergMarquardtDense) require the graph to be
 * in LS mode and only LS objectives are allowed.
 *
 * @see OptimizationProblemInterface StandardOptimizationProblem SolverInterface
 * @see VertexInterface EdgeInterface BaseEdge DiscretizationGrid
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class BaseHyperGraphOptimizationProblem : public OptimizationProblemInterface
{
 public:
    using Ptr = std::shared_ptr<BaseHyperGraphOptimizationProblem>;

    BaseHyperGraphOptimizationProblem() = default;
    BaseHyperGraphOptimizationProblem(OptimizationEdgeSet::Ptr edges, VertexSetInterface::Ptr vertices) : _graph(edges, vertices) {}

    virtual Ptr getInstance() const { return std::make_shared<BaseHyperGraphOptimizationProblem>(); }

    //! Get access to the accociated factory
    static Factory<BaseHyperGraphOptimizationProblem>& getFactory() { return Factory<BaseHyperGraphOptimizationProblem>::instance(); }

    void setGraph(OptimizationEdgeSet::Ptr edges, VertexSetInterface::Ptr vertices)
    {
        _graph.setEdgeSet(edges);
        _graph.setVertexSet(vertices);
    }

    const HyperGraph& getGraph() const { return _graph; }
    HyperGraph& getGraph() { return _graph; }

    virtual void precomputeGraphQuantities();
    virtual void precomputeVertexQuantities();
    virtual void precomputeEdgeQuantities();

    /** @name Specify the dimension of the optimization problem  */
    //@{

    //! Total dimension of objective function terms
    int getNonLsqObjectiveDimension() override
    {
        if (!_graph_precomputed) precomputeGraphQuantities();
        return _dim_non_lsq_obj;
    }
    int getLsqObjectiveDimension() override
    {
        if (!_graph_precomputed) precomputeGraphQuantities();
        return _dim_lsq_obj;
    }
    int getObjectiveDimension() override
    {
        if (!_graph_precomputed) precomputeGraphQuantities();
        return (int)(_dim_lsq_obj > 0 || _dim_non_lsq_obj > 0);
    }
    //! Total dimension of equality constraints
    int getEqualityDimension() override
    {
        if (!_graph_precomputed) precomputeGraphQuantities();
        return _dim_eq;
    }
    //! Total dimension of general inequality constraints
    int getInequalityDimension() override
    {
        if (!_graph_precomputed) precomputeGraphQuantities();
        return _dim_ineq;
    }

    //@}

    /** @name main equations of the optimization problem  */
    //@{

    // implmements interface method
    double computeValueNonLsqObjective() override;

    // implements interface method
    void computeValuesLsqObjective(Eigen::Ref<Eigen::VectorXd> values) override;

    // implements interface method
    double computeValueObjective() override;

    // implements interface method
    void computeValuesEquality(Eigen::Ref<Eigen::VectorXd> values) override;

    // implements interface method
    void computeValuesInequality(Eigen::Ref<Eigen::VectorXd> values) override;

    // implements interface method
    void computeValues(double& non_lsq_obj_value, Eigen::Ref<Eigen::VectorXd> lsq_obj_values, Eigen::Ref<Eigen::VectorXd> eq_values,
                       Eigen::Ref<Eigen::VectorXd> ineq_values) override;

    //@}

    /** @name Access parameter vector and bounds  */
    //@{

    // implements interface method
    double getParameterValue(int idx) override { return _graph.getVertexSetRaw()->getParameterValue(idx); }
    // implements interface method
    void setParameterValue(int idx, double x) override { _graph.getVertexSetRaw()->setParameterValue(idx, x); }
    // implements interface method
    void getParameterVector(Eigen::Ref<Eigen::VectorXd> x) override { _graph.getVertexSetRaw()->getParameterVector(x); }
    // implements interface method
    void setParameterVector(const Eigen::Ref<const Eigen::VectorXd>& x) override { _graph.getVertexSetRaw()->setParameterVector(x); }

    // implements interface method
    void getBounds(Eigen::Ref<Eigen::VectorXd> lb, Eigen::Ref<Eigen::VectorXd> ub) override { _graph.getVertexSetRaw()->getBounds(lb, ub); }
    // implements interface method
    void setBounds(const Eigen::Ref<const Eigen::VectorXd>& lb, const Eigen::Ref<const Eigen::VectorXd>& ub) override
    {
        _graph.getVertexSetRaw()->setBounds(lb, ub);
    }
    // implements interface method
    double getLowerBound(int idx) override { return _graph.getVertexSetRaw()->getLowerBound(idx); }
    // implements interface method
    double getUpperBound(int idx) override { return _graph.getVertexSetRaw()->getUpperBound(idx); }
    // implements interface method
    void setLowerBound(int idx, double lb) override { _graph.getVertexSetRaw()->setLowerBound(idx, lb); }
    // implements interface method
    void setUpperBound(int idx, double ub) override { _graph.getVertexSetRaw()->setUpperBound(idx, ub); }

    //@}

    /** @name Interface implementations  */
    //@{

    // implements interface method
    int getParameterDimension() override
    {
        if (!_graph_precomputed) precomputeGraphQuantities();
        return _dim_param;
    }

    // implements interface method
    void applyIncrement(const Eigen::Ref<const Eigen::VectorXd>& increment) override { _graph.getVertexSetRaw()->applyIncrementNonFixed(increment); }
    // implements interface method
    void applyIncrement(int idx, double increment) override { _graph.getVertexSetRaw()->applyIncrementNonFixed(idx, increment); }

    // implements interface method
    void backupParameters() override { _graph.getVertexSetRaw()->backupParametersActiveVertices(); }
    // implements interface method
    void restoreBackupParameters(bool keep_backup) override { _graph.getVertexSetRaw()->restoreBackupParametersActiveVertices(keep_backup); }
    // implements interface method
    void discardBackupParameters(bool all = false) override { _graph.getVertexSetRaw()->discardBackupParametersActiveVertices(all); }

    //@}

    /** @name Specify properties of the optimization problem  */
    //@{

    // implements interface method
    bool isLeastSquaresProblem() const override { return _graph.getEdgeSetRaw()->hasOnlyLeastSquaresObjectives(); }

    //@}

    /** @name Methods for dealing with bounds   */
    //@{

    // implements interface method
    int finiteCombinedBoundsDimension() override;
    // implements interface method
    int finiteBoundsDimension() override;
    // implements interface method
    void computeValuesActiveInequality(Eigen::Ref<Eigen::VectorXd> values, double weight = 1.0) override;
    // implements interface method
    void computeDistanceFiniteCombinedBounds(Eigen::Ref<Eigen::VectorXd> values) override;
    // implements interface method
    void computeLowerAndUpperBoundDiff(Eigen::Ref<Eigen::VectorXd> lb_minus_x, Eigen::Ref<Eigen::VectorXd> ub_minus_x) override;
    // implements interface method
    void getParametersAndBoundsFinite(Eigen::Ref<Eigen::VectorXd> lb_finite_bounds, Eigen::Ref<Eigen::VectorXd> ub_finite_bounds,
                                      Eigen::Ref<Eigen::VectorXd> x_finite_bounds) override;

    void computeDenseJacobianFiniteCombinedBounds(Eigen::Ref<Eigen::MatrixXd> jacobian, double weight = 1.0) override;
    int computeSparseJacobianFiniteCombinedBoundsNNZ() override;
    void computeSparseJacobianFiniteCombinedBoundsStructure(Eigen::Ref<Eigen::VectorXi> i_row, Eigen::Ref<Eigen::VectorXi> j_col) override;
    void computeSparseJacobianFiniteCombinedBoundsValues(Eigen::Ref<Eigen::VectorXd> values, double weight = 1.0) override;

    void computeDenseJacobianFiniteCombinedBoundsIdentity(Eigen::Ref<Eigen::MatrixXd> jacobian) override;

    //@}

    //! Check if a function taking the parameter value and unfixed-idx is true for all unfixed parameter values.
    bool checkIfAllUnfixedParam(std::function<bool(double, int)> fun);

    void clear() override;

#ifdef MESSAGE_SUPPORT
    // implements interface method
    virtual void toMessage(corbo::messages::HyperGraphOptimizationProblem& message) const {}
    // implements interface method
    virtual void fromMessage(const corbo::messages::HyperGraphOptimizationProblem& message, std::stringstream* issues = nullptr) {}
#endif

 protected:
    HyperGraph _graph;
    bool _graph_precomputed = false;

    int _dim_param       = 0;
    int _dim_non_lsq_obj = 0;
    int _dim_lsq_obj     = 0;
    int _dim_eq          = 0;
    int _dim_ineq        = 0;
};

using HyperGraphOptimizationProblemFactory = Factory<BaseHyperGraphOptimizationProblem>;
#define FACTORY_REGISTER_HYPER_GRAPH_OPTIMIZATION_PROBLEM(type) FACTORY_REGISTER_OBJECT_ID(type, BaseHyperGraphOptimizationProblem, 100)

}  // namespace corbo

#endif  // SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_HYPER_GRAPH_OPTIMIZATION_PROBLEM_BASE_H_
