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

#ifndef SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_GENERIC_EDGE_H_
#define SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_GENERIC_EDGE_H_

#include <corbo-optimization/hyper_graph/edge.h>

#include <corbo-optimization/hyper_graph/scalar_vertex.h>
#include <corbo-optimization/hyper_graph/vector_vertex.h>

#include <functional>

namespace corbo {

/**
 * @brief Generic edge for functions \f$ \mathbb{R}^n \to \mathbb{R} \f$.
 *
 * @ingroup optimization hyper-graph
 *
 * This edge can be used for passing a generic function (e.g. a lambda)
 * to the optimization problem. Jacobians and Hessians are calculated numerically
 * if necessary.
 *
 * This edge should be only used for Rapid-Prototyping puroposes, due to the
 * computational overhead that rises from copying and calling the external function resp. functor.
 *
 * You can pass a function that accepts a EdgeGenericScalarFun::VertexContainer as
 * argument. Each vertex stored in that container corresponds to the vertex passed using the
 * class constructor (in the same order). The corresponding Vertex types are specified as template
 * arguments of this class as well.
 *
 * Example usage for defining (x-1)^2 (Note, keep in mind the special definitions for squaring within different solvers):
 * @code
 *      VectorVertex<1> x;
 *      using MyEdgeT = EdgeGenericScalarFun<StateVertex<1>>;
 *      auto fun = [] (MyEdgeT::VertexContainer& vertices) {return vertices.at(0)->getData(0)-1;}; // Note, we define only (x-1) here
 *      MyEdgeT* my_edge = new MyEdgeT(fun, false, x);
 *      // now add to the hyper-graph
 * @endcode
 *
 * @see EdgeGenericVectorFun BaseEdge
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
template <class... VerticesT>
class EdgeGenericScalarFun : public Edge<VerticesT...>
{
    // Use dependent names and members from the parent template class.
 protected:
    using Edge<VerticesT...>::_vertices;

 public:
    using typename Edge<VerticesT...>::VertexContainer;
    using ExtFunT   = std::function<double(const VertexContainer&)>;
    using ExtJacobT = std::function<void(const VertexContainer&, int, Eigen::Ref<Eigen::MatrixXd>, const double*)>;

    using Ptr  = std::shared_ptr<EdgeGenericScalarFun>;
    using UPtr = std::unique_ptr<EdgeGenericScalarFun>;

    //! Construct generic scalar function edge by passing the function object and all vertices
    EdgeGenericScalarFun(const ExtFunT& fun, bool lsq_form, VerticesT&... vertices) : Edge<VerticesT...>(vertices...), _fun(fun), _lsq_form(lsq_form)
    {
    }

    //! Construct generic scalar function edge by passing the function object, a jacobian, and all vertices
    EdgeGenericScalarFun(const ExtFunT& fun, const ExtJacobT& jacobian, bool lsq_form, VerticesT&... vertices)
        : Edge<VerticesT...>(vertices...), _fun(fun), _jacob(jacobian), _lsq_form(lsq_form)
    {
    }

    int getDimension() const override { return 1; }

    // implements interface method
    bool isLinear() const override { return false; }
    // implements interface method
    bool isLeastSquaresForm() const override { return _lsq_form; }
    // implements interface method
    bool providesJacobian() const override { return (bool)_jacob; }

    // implements interface method
    void computeValues(Eigen::Ref<Eigen::VectorXd> values) override { values[0] = _fun(_vertices); }
    // implements interface method
    void computeJacobian(int vtx_idx, Eigen::Ref<Eigen::MatrixXd> block_jacobian, const double* multiplier = nullptr) override
    {
        if (_jacob)
            _jacob(_vertices, vtx_idx, block_jacobian, multiplier);
        else
            BaseEdge::computeJacobian(vtx_idx, block_jacobian, multiplier);
    }

    // implements interface method
    void setLinear(bool linear) { _linear = linear; }

 protected:
    ExtFunT _fun;  //!< Store the external function or functor
    ExtJacobT _jacob;
    bool _lsq_form;
    bool _linear = false;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

/**
 * @brief Generic edge for functions \f$ \mathbb{R}^n \to \mathbb{R}^D \f$.
 *
 * @ingroup optimization hyper-graph
 *
 * This edge can be used for passing a generic function (e.g. a lambda)
 * to the optimization problem. Jacobians and Hessians are calculated numerically
 * if necessary.
 *
 * This edge should be only used for Rapid-Prototyping puroposes, due to the
 * computational overhead that rises from copying and calling the external function resp. functor.
 *
 * You can pass a function that accepts a EdgeGenericScalarFun::VertexContainer as
 * argument. Each vertex stored in that container corresponds to the vertex passed using the
 * class constructor (in the same order). The corresponding Vertex types are specified as template
 * arguments of this class as well.
 *
 * The return type of the vector valued external function should be a Eigen::Matrix<double,D,1> type.
 *
 * @see EdgeGenericScalarFun BaseEdge
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
template <int D, class... VerticesT>
class EdgeGenericVectorFun : public Edge<VerticesT...>
{
    // Use dependent names and members from the parent template class.
 public:
    using typename Edge<VerticesT...>::VertexContainer;

    using Ptr  = std::shared_ptr<EdgeGenericVectorFun>;
    using UPtr = std::unique_ptr<EdgeGenericVectorFun>;

 protected:
    using Edge<VerticesT...>::_vertices;

 public:
    using ErrorVector = Eigen::Matrix<double, D, 1>;
    using ExtFunT     = std::function<void(const VertexContainer&, Eigen::Ref<ErrorVector>)>;
    using ExtJacobT   = std::function<void(const VertexContainer&, int, Eigen::Ref<Eigen::MatrixXd>, const double*)>;

    //! Construct generic vector function edge by passing the dimension D, the function object and all vertices
    EdgeGenericVectorFun(const ExtFunT& fun, bool lsq_form, VerticesT&... vertices) : Edge<VerticesT...>(vertices...), _fun(fun), _lsq_form(lsq_form)
    {
    }

    //! Construct generic vector function edge by passing the dimension D, the function object, a jacobian and all vertices
    EdgeGenericVectorFun(const ExtFunT& fun, const ExtJacobT& jacobian, bool lsq_form, VerticesT&... vertices)
        : Edge<VerticesT...>(vertices...), _fun(fun), _jacob(jacobian), _lsq_form(lsq_form)
    {
    }

    int getDimension() const override { return D; }
    // implements interface method
    bool isLinear() const override { return false; }
    // implements interface method
    bool isLeastSquaresForm() const override { return _lsq_form; }
    // implements interface method
    bool providesJacobian() const override { return (bool)_jacob; }

    // implements interface method
    void computeValues(Eigen::Ref<Eigen::VectorXd> values) override { _fun(_vertices, values); }
    // implements interface method
    void computeJacobian(int vtx_idx, Eigen::Ref<Eigen::MatrixXd> block_jacobian, const double* multiplier) override
    {
        if (_jacob)
            _jacob(_vertices, vtx_idx, block_jacobian, multiplier);
        else
            BaseEdge::computeJacobian(vtx_idx, block_jacobian, multiplier);
    }

    // implements interface method
    void setLinear(bool linear) { _linear = linear; }

 protected:
    ExtFunT _fun;  //!< Store the external function or functor
    ExtJacobT _jacob;
    bool _lsq_form;
    bool _linear = false;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template <class... VerticesT>
class MixedEdgeGenericVectorFun : public MixedEdge<VerticesT...>
{
    // Use dependent names and members from the parent template class.
 public:
    using typename MixedEdge<VerticesT...>::VertexContainer;

    using Ptr  = std::shared_ptr<MixedEdgeGenericVectorFun>;
    using UPtr = std::unique_ptr<MixedEdgeGenericVectorFun>;

 protected:
    using MixedEdge<VerticesT...>::_vertices;

 public:
    using ExtPrecomputeT = std::function<void(const VertexContainer&)>;
    using ExtFunT        = std::function<void(const VertexContainer&, Eigen::Ref<Eigen::VectorXd>)>;
    // using ExtJacobT   = std::function<void(const VertexContainer&, int, Eigen::Ref<Eigen::MatrixXd>, const double*)>;

    //! Construct generic vector function edge by passing the dimension D, the function object and all vertices
    MixedEdgeGenericVectorFun(int obj_dim, int eq_dim, int ineq_dim, const ExtPrecomputeT& precompute_fun, const ExtFunT& obj_fun,
                              const ExtFunT& eq_fun, const ExtFunT& ineq_fun, VerticesT&... vertices)
        : MixedEdge<VerticesT...>(vertices...),
          _preqcompute_fun(precompute_fun),
          _obj_fun(obj_fun),
          _eq_fun(eq_fun),
          _ineq_fun(ineq_fun),
          _obj_dim(obj_dim),
          _eq_dim(eq_dim),
          _ineq_dim(ineq_dim)
    {
    }

    //! Construct generic vector function edge by passing the dimension D, the function object, a jacobian and all vertices
    // MixedEdgeGenericVectorFun(const ExtFunT& fun, const ExtJacobT& jacobian, bool lsq_form, VerticesT&... vertices)
    //: Edge<VerticesT...>(vertices...), _fun(fun), _jacob(jacobian), _lsq_form(lsq_form)
    //{
    //}

    int getObjectiveDimension() const override { return _obj_dim; }
    int getEqualityDimension() const override { return _eq_dim; }
    int getInequalityDimension() const override { return _ineq_dim; }

    // implements interface method
    // bool isLinear() const override { return false; }
    // implements interface method
    bool isObjectiveLeastSquaresForm() const override { return _obj_lsq_form; }
    // implements interface method
    // bool providesJacobian() const override { return (bool)_jacob; }

    // implements interface method
    void precompute() override { _preqcompute_fun(_vertices); }
    void computeObjectiveValues(Eigen::Ref<Eigen::VectorXd> obj_values) override { _obj_fun(_vertices, obj_values); }
    void computeEqualityValues(Eigen::Ref<Eigen::VectorXd> eq_values) override { _eq_fun(_vertices, eq_values); }
    void computeInequalityValues(Eigen::Ref<Eigen::VectorXd> ineq_values) override { _ineq_fun(_vertices, ineq_values); }
    // implements interface method
    // void computeJacobian(int vtx_idx, Eigen::Ref<Eigen::MatrixXd> block_jacobian, const double* multiplier) override
    //{
    //   if (_jacob)
    //       _jacob(_vertices, vtx_idx, block_jacobian, multiplier);
    //   else
    //       BaseEdge::computeJacobian(vtx_idx, block_jacobian, multiplier);
    //}

    // implements interface method
    void setLinear(bool linear) { _linear = linear; }

    void setObjectiveLsqForm(bool obj_lsq) { _obj_lsq_form = obj_lsq; }

 protected:
    ExtPrecomputeT _preqcompute_fun;
    ExtFunT _obj_fun;  //!< Store the external function or functor
    ExtFunT _eq_fun;
    ExtFunT _ineq_fun;
    // ExtJacobT _obj_jacob;
    bool _obj_lsq_form = false;
    bool _linear       = false;
    int _obj_dim       = 0;
    int _eq_dim        = 0;
    int _ineq_dim      = 0;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template <typename T, void (T::*ComputeMethod)(int, const Eigen::Ref<const Eigen::VectorXd>&, Eigen::Ref<Eigen::VectorXd>) const>
class UnaryVectorVertexEdge : public Edge<VectorVertex>
{
 public:
    using Ptr = std::shared_ptr<UnaryVectorVertexEdge>;

    explicit UnaryVectorVertexEdge(int dim, int k, const T& fun_obj, VectorVertex& vertex, bool is_linear = false, bool is_lsq = false)
        : Edge<VectorVertex>(vertex), _dimension(dim), _k(k), _is_linear(is_linear), _is_lsq(is_lsq), _fun_obj(fun_obj)
    {
    }

    int getDimension() const override { return _dimension; }

    // implements interface method
    bool isLinear() const override { return _is_linear; }

    // implements interface method
    bool isLeastSquaresForm() const override { return _is_lsq; }

    // implements interface method
    bool providesJacobian() const override { return false; }

    // implements interface method
    void computeValues(Eigen::Ref<Eigen::VectorXd> values) override
    {
        const VectorVertex* vertex = static_cast<const VectorVertex*>(_vertices[0]);
        (_fun_obj.*ComputeMethod)(_k, vertex->values(), values);
    }

 private:
    int _dimension  = 0;
    int _k          = 0;
    bool _is_linear = false;
    bool _is_lsq    = false;

    const T& _fun_obj;
};

template <typename T, void (T::*ComputeMethod)(int, double, Eigen::Ref<Eigen::VectorXd>) const>
class UnaryScalarVertexEdge : public Edge<ScalarVertex>
{
 public:
    using Ptr = std::shared_ptr<UnaryScalarVertexEdge>;

    explicit UnaryScalarVertexEdge(int dim, int k, const T& fun_obj, ScalarVertex& vertex, bool is_linear, bool is_lsq)
        : Edge<ScalarVertex>(vertex), _dimension(dim), _k(k), _is_linear(is_linear), _is_lsq(is_lsq), _fun_obj(fun_obj)
    {
    }

    int getDimension() const override { return _dimension; }

    // implements interface method
    bool isLinear() const override { return _is_linear; }

    // implements interface method
    bool isLeastSquaresForm() const override { return _is_lsq; }

    // implements interface method
    bool providesJacobian() const override { return false; }

    // implements interface method
    void computeValues(Eigen::Ref<Eigen::VectorXd> values) override
    {
        const ScalarVertex* vertex = static_cast<const ScalarVertex*>(_vertices[0]);
        (_fun_obj.*ComputeMethod)(_k, vertex->value(), values);
    }

 private:
    int _dimension  = 0;
    int _k          = 0;
    bool _is_linear = false;
    bool _is_lsq    = false;

    const T& _fun_obj;
};

template <typename T, void (T::*ComputeMethod)(int, const Eigen::Ref<const Eigen::VectorXd>&, const Eigen::Ref<const Eigen::VectorXd>&,
                                               Eigen::Ref<Eigen::VectorXd>) const>
class BinaryVectorVertexEdge : public Edge<VectorVertex, VectorVertex>
{
 public:
    using Ptr = std::shared_ptr<BinaryVectorVertexEdge>;

    explicit BinaryVectorVertexEdge(int dim, int k, const T& fun_obj, VectorVertex& vertex1, VectorVertex& vertex2, bool is_linear, bool is_lsq)
        : Edge<VectorVertex, VectorVertex>(vertex1, vertex2), _dimension(dim), _k(k), _is_linear(is_linear), _is_lsq(is_lsq), _fun_obj(fun_obj)
    {
    }

    int getDimension() const override { return _dimension; }

    // implements interface method
    bool isLinear() const override { return _is_linear; }

    // implements interface method
    bool isLeastSquaresForm() const override { return _is_lsq; }

    // implements interface method
    bool providesJacobian() const override { return false; }

    // implements interface method
    void computeValues(Eigen::Ref<Eigen::VectorXd> values) override
    {
        const VectorVertex* vertex1 = static_cast<const VectorVertex*>(_vertices[0]);
        const VectorVertex* vertex2 = static_cast<const VectorVertex*>(_vertices[1]);
        (_fun_obj.*ComputeMethod)(_k, vertex1->values(), vertex2->values(), values);
    }

 private:
    int _dimension  = 0;
    int _k          = 0;
    bool _is_linear = false;
    bool _is_lsq    = false;

    const T& _fun_obj;
};

template <typename T, void (T::*ComputeMethod)(int, const Eigen::Ref<const Eigen::VectorXd>&, double, Eigen::Ref<Eigen::VectorXd>) const>
class BinaryVectorScalarVertexEdge : public Edge<VectorVertex, ScalarVertex>
{
 public:
    using Ptr = std::shared_ptr<BinaryVectorScalarVertexEdge>;

    explicit BinaryVectorScalarVertexEdge(int dim, int k, const T& fun_obj, VectorVertex& vertex1, ScalarVertex& vertex2, bool is_linear, bool is_lsq)
        : Edge<VectorVertex, ScalarVertex>(vertex1, vertex2), _dimension(dim), _k(k), _is_linear(is_linear), _is_lsq(is_lsq), _fun_obj(fun_obj)
    {
    }

    int getDimension() const override { return _dimension; }

    // implements interface method
    bool isLinear() const override { return _is_linear; }

    // implements interface method
    bool isLeastSquaresForm() const override { return _is_lsq; }

    // implements interface method
    bool providesJacobian() const override { return false; }

    // implements interface method
    void computeValues(Eigen::Ref<Eigen::VectorXd> values) override
    {
        const VectorVertex* vertex1 = static_cast<const VectorVertex*>(_vertices[0]);
        const ScalarVertex* vertex2 = static_cast<const ScalarVertex*>(_vertices[1]);
        (_fun_obj.*ComputeMethod)(_k, vertex1->values(), vertex2->value(), values);
    }

 private:
    int _dimension  = 0;
    int _k          = 0;
    bool _is_linear = false;
    bool _is_lsq    = false;

    const T& _fun_obj;
};

template <typename T, void (T::*ComputeMethod)(int, const Eigen::Ref<const Eigen::VectorXd>&, const Eigen::Ref<const Eigen::VectorXd>&, double,
                                               Eigen::Ref<Eigen::VectorXd>) const>
class TernaryVectorScalarVertexEdge : public Edge<VectorVertex, VectorVertex, ScalarVertex>
{
 public:
    using Ptr = std::shared_ptr<TernaryVectorScalarVertexEdge>;

    explicit TernaryVectorScalarVertexEdge(int dim, int k, const T& fun_obj, VectorVertex& vec_vtx1, VectorVertex& vec_vtx2, ScalarVertex& scalar_vtx,
                                           bool is_linear, bool is_lsq)
        : Edge<VectorVertex, VectorVertex, ScalarVertex>(vec_vtx1, vec_vtx2, scalar_vtx),
          _dimension(dim),
          _k(k),
          _is_linear(is_linear),
          _is_lsq(is_lsq),
          _fun_obj(fun_obj)
    {
    }

    int getDimension() const override { return _dimension; }

    // implements interface method
    bool isLinear() const override { return _is_linear; }

    // implements interface method
    bool isLeastSquaresForm() const override { return _is_lsq; }

    // implements interface method
    bool providesJacobian() const override { return false; }

    // implements interface method
    void computeValues(Eigen::Ref<Eigen::VectorXd> values) override
    {
        const VectorVertex* vec_vtx1   = static_cast<const VectorVertex*>(_vertices[0]);
        const VectorVertex* vec_vtx2   = static_cast<const VectorVertex*>(_vertices[1]);
        const ScalarVertex* scalar_vtx = static_cast<const ScalarVertex*>(_vertices[2]);
        (_fun_obj.*ComputeMethod)(_k, vec_vtx1->values(), vec_vtx2->values(), scalar_vtx->value(), values);
    }

 private:
    int _dimension  = 0;
    int _k          = 0;
    bool _is_linear = false;
    bool _is_lsq    = false;

    const T& _fun_obj;
};

}  // namespace corbo

#endif  // SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_GENERIC_EDGE_H_
