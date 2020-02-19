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

#ifndef SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_VECTOR_VERTEX_H_
#define SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_VECTOR_VERTEX_H_

#include <corbo-optimization/hyper_graph/vertex_interface.h>

#include <corbo-core/types.h>

#include <memory>
#include <vector>

namespace corbo {

/**
 * @brief Vertex implementation that stores an Eigen::VectorXd (dynamic dimension)
 *
 * @ingroup optimization hyper-graph
 *
 * The vertex can be either completely fixed or unfixed.
 * In order to partially fix components of the underlying vector
 * refer to class PartiallyFixedVectorVertex.
 *
 * @see VertexInterface PartiallyFixedVectorVertex HyperGraph EdgeInterface
 *      ScalarVertex
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class VectorVertex : public VertexInterface
{
 public:
    using Ptr  = std::shared_ptr<VectorVertex>;
    using UPtr = std::unique_ptr<VectorVertex>;

    //! Default constructor
    VectorVertex() = default;

    explicit VectorVertex(bool fixed) : _fixed(fixed) {}

    //! Construct and allocate memory for a given dimension
    explicit VectorVertex(int dimension, bool fixed = false)
        : _values(Eigen::VectorXd::Zero(dimension)),
          _lb(Eigen::VectorXd::Constant(dimension, -CORBO_INF_DBL)),
          _ub(Eigen::VectorXd::Constant(dimension, CORBO_INF_DBL)),
          _finite_lb_bounds(false),
          _finite_ub_bounds(false),
          _fixed(fixed)
    {
    }

    //! Construct vertex with given values
    explicit VectorVertex(const Eigen::Ref<const Eigen::VectorXd>& values, bool fixed = false)
        : _values(values),
          _lb(Eigen::VectorXd::Constant(values.size(), -CORBO_INF_DBL)),
          _ub(Eigen::VectorXd::Constant(values.size(), CORBO_INF_DBL)),
          _finite_lb_bounds(false),
          _finite_ub_bounds(false),
          _fixed(fixed)
    {
    }

    //! Construct vertex with given values, lower and upper bounds
    explicit VectorVertex(const Eigen::Ref<const Eigen::VectorXd>& values, const Eigen::Ref<const Eigen::VectorXd>& lb,
                          const Eigen::Ref<const Eigen::VectorXd>& ub, bool fixed = false)
        : _values(values), _fixed(fixed)
    {
        assert(values.size() == lb.size());
        assert(values.size() == ub.size());
        setLowerBounds(lb);
        setUpperBounds(ub);
    }

    // implements interface method
    int getDimension() const override { return _values.rows(); }
    // implements interface method
    int getDimensionUnfixed() const override { return _fixed ? 0 : _values.rows(); }

    //! Change the dimension of the vertex (lower and upper bounds needs to be redefined)
    virtual void setDimension(int dim)
    {
        _values           = Eigen::VectorXd::Zero(dim);
        _lb               = Eigen::VectorXd::Constant(dim, -CORBO_INF_DBL);
        _ub               = Eigen::VectorXd::Constant(dim, CORBO_INF_DBL);
        _finite_lb_bounds = _finite_ub_bounds = false;
    }

    // implements interface method
    void plus(int idx, double inc) override { _values[idx] += inc; }
    // implements interface method
    void plus(const double* inc) override { _values += Eigen::Map<const Eigen::VectorXd>(inc, getDimension()); }
    // implements interface method
    void plusUnfixed(const double* inc) override { plus(inc); }

    // implements interface method
    const double* getData() const override { return _values.data(); }
    // implements interface method

    double* getDataRaw() override { return _values.data(); }

    // implements interface method
    void setData(int idx, double data) override { _values[idx] = data; }

    //! Set values and bounds at once
    virtual void set(const Eigen::Ref<const Eigen::VectorXd>& values, const Eigen::Ref<const Eigen::VectorXd>& lb,
                     const Eigen::Ref<const Eigen::VectorXd>& ub, bool fixed = false)
    {
        assert(values.size() == lb.size());
        assert(values.size() == ub.size());
        _values = values;
        setLowerBounds(lb);
        setUpperBounds(ub);

        setFixed(false);
    }

    // implements interface method
    bool hasFixedComponents() const override { return _fixed; }
    // implements interface method
    bool isFixedComponent(int /*idx*/) const override { return _fixed; }

    // implements interface method
    void setLowerBounds(const Eigen::Ref<const Eigen::VectorXd>& lb) override
    {
        _lb               = lb;
        _finite_lb_bounds = (_lb.array() > -CORBO_INF_DBL).any();
    }
    // implements interface method
    void setLowerBound(int idx, double lb) override
    {
        _lb[idx]          = lb;
        _finite_lb_bounds = (_lb.array() > -CORBO_INF_DBL).any();
    }
    // implements interface method
    void setUpperBounds(const Eigen::Ref<const Eigen::VectorXd>& ub) override
    {
        _ub               = ub;
        _finite_ub_bounds = (_ub.array() < CORBO_INF_DBL).any();
    }
    // implements interface method
    void setUpperBound(int idx, double ub) override
    {
        _ub[idx]          = ub;
        _finite_ub_bounds = (_ub.array() < CORBO_INF_DBL).any();
    }
    // implements interface method
    bool hasFiniteBounds() const override { return _finite_lb_bounds || _finite_ub_bounds; }
    // implements interface method
    bool hasFiniteLowerBounds() const override { return _finite_lb_bounds; }
    // implements interface method
    bool hasFiniteUpperBounds() const override { return _finite_ub_bounds; }
    // implements interface method
    bool hasFiniteLowerBound(int idx) const override
    {
        assert(idx < _lb.size());
        return _lb[idx] > -CORBO_INF_DBL;
    }
    // implements interface method
    bool hasFiniteUpperBound(int idx) const override
    {
        assert(idx < _lb.size());
        return _ub[idx] < CORBO_INF_DBL;
    }
    // implements interface method
    int getNumberFiniteLowerBounds(bool unfixed_only) const override
    {
        if (unfixed_only && _fixed)
            return 0;
        else
            return (_lb.array() > -CORBO_INF_DBL).count();
    }
    // implements interface method
    int getNumberFiniteUpperBounds(bool unfixed_only) const override
    {
        if (unfixed_only && _fixed)
            return 0;
        else
            return (_ub.array() < CORBO_INF_DBL).count();
    }
    int getNumberFiniteBounds(bool unfixed_only) const override
    {
        if (unfixed_only && _fixed)
            return 0;
        else
            return (_ub.array() < CORBO_INF_DBL || _lb.array() > -CORBO_INF_DBL).count();
    }

    //! Set complete vertex to fixed (and hence skip during optimization)
    virtual void setFixed(bool fixed) { _fixed = fixed; }

    // implements interface method
    const double* getLowerBounds() const override { return _lb.data(); }
    // implements interface method
    const double* getUpperBounds() const override { return _ub.data(); }

    // Backup values
    // implements interface method
    void push() override { _backup.push_back(_values); }
    // implements interface method
    void pop() override
    {
        top();
        _backup.pop_back();
    }
    // implements interface method
    void top() override
    {
        assert(!_backup.empty());
        _values = _backup.back();
    }
    // implements interface method
    void discardTop() override { _backup.pop_back(); }
    // implements interface method
    void clear() override { _backup.clear(); }
    // implements interface method
    int getNumBackups() const override { return (int)_backup.size(); }

    //! Read-access to the underlying value vector
    const Eigen::VectorXd& values() const { return _values; }
    //! Write-access to the underlying value vector
    Eigen::VectorXd& values() { return _values; }

    //! Read-access to the underlying lower bound vector
    const Eigen::VectorXd& lowerBound() const { return _lb; }
    //! Read-access to the underlying upper bound vector
    const Eigen::VectorXd& upperBound() const { return _ub; }

 protected:
    Eigen::VectorXd _values;
    Eigen::VectorXd _lb;
    Eigen::VectorXd _ub;
    bool _finite_lb_bounds = false;
    bool _finite_ub_bounds = false;

    bool _fixed = false;

    std::vector<Eigen::VectorXd> _backup;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

/**
 * @brief Vector based vertex with support for partially fixed components
 *
 * @ingroup optimization hyper-graph
 *
 * The vertex extends VectorVertex by allowing the user to
 * partially fix components of the underlying vector.
 *
 * @see VertexInterface VectorVertex HyperGraph EdgeInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class PartiallyFixedVectorVertex : public VectorVertex
{
 public:
    using Ptr  = std::shared_ptr<PartiallyFixedVectorVertex>;
    using UPtr = std::unique_ptr<PartiallyFixedVectorVertex>;

    //! Default constructor
    PartiallyFixedVectorVertex() = default;

    explicit PartiallyFixedVectorVertex(int dimension)
        : VectorVertex(dimension), _fixed(Eigen::Array<bool, -1, 1>::Constant(dimension, false)), _num_unfixed(dimension)
    {
    }

    //! Construct and allocate memory for a given dimension
    explicit PartiallyFixedVectorVertex(int dimension, const Eigen::Ref<const Eigen::Array<bool, -1, 1>>& fixed)
        : VectorVertex(dimension), _fixed(fixed), _num_unfixed(dimension)
    {
    }

    //! Construct vertex with given values
    explicit PartiallyFixedVectorVertex(const Eigen::Ref<const Eigen::VectorXd>& values)
        : VectorVertex(values), _fixed(Eigen::Array<bool, -1, 1>::Constant(values.size(), false)), _num_unfixed(values.size())
    {
    }

    //! Construct vertex with given values and fixed components
    explicit PartiallyFixedVectorVertex(const Eigen::Ref<const Eigen::VectorXd>& values, const Eigen::Ref<const Eigen::Array<bool, -1, 1>>& fixed)
        : VectorVertex(values), _fixed(fixed), _num_unfixed(fixed.size() - fixed.count())
    {
    }

    //! Construct vertex with given values, lower and upper bounds
    explicit PartiallyFixedVectorVertex(const Eigen::Ref<const Eigen::VectorXd>& values, const Eigen::Ref<const Eigen::VectorXd>& lb,
                                        const Eigen::Ref<const Eigen::VectorXd>& ub)
        : VectorVertex(values, lb, ub), _fixed(Eigen::Array<bool, -1, 1>::Constant(values.size(), false)), _num_unfixed(values.size())
    {
    }

    // implements interface method
    int getDimensionUnfixed() const override { return _num_unfixed; }

    // implements parent method
    void setDimension(int dim) override
    {
        VectorVertex::setDimension(dim);
        _fixed.setConstant(dim, false);
        _num_unfixed = dim;
    }

    //! Set values and bounds at once
    void set(const Eigen::Ref<const Eigen::VectorXd>& values, const Eigen::Ref<const Eigen::VectorXd>& lb,
             const Eigen::Ref<const Eigen::VectorXd>& ub, bool fixed = false) override
    {
        assert(values.size() == lb.size());
        assert(values.size() == ub.size());
        _values = values;
        setLowerBounds(lb);
        setUpperBounds(ub);

        setFixed(fixed);
    }

    //! Set values and bounds at once (overload with fixed vector)
    void set(const Eigen::Ref<const Eigen::VectorXd>& values, const Eigen::Ref<const Eigen::VectorXd>& lb,
             const Eigen::Ref<const Eigen::VectorXd>& ub, const Eigen::Ref<const Eigen::Array<bool, -1, 1>>& fixed)
    {
        assert(values.size() == lb.size());
        assert(values.size() == ub.size());
        _values = values;
        setLowerBounds(lb);
        setUpperBounds(ub);

        setFixed(fixed);
    }

    //! Set component with idx (0 <= idx < dimension()) to (un)fixed
    void setFixed(int idx, bool fixed)
    {
        _fixed[idx]  = fixed;
        _num_unfixed = getDimension() - _fixed.count();
    }

    //! Set logical array [dimension() x 1] in order to fix selected components
    void setFixed(const Eigen::Ref<const Eigen::Array<bool, -1, 1>>& fixed)
    {
        _fixed       = fixed;
        _num_unfixed = getDimension() - _fixed.count();
    }

    // implements interface method
    void setFixed(bool fixed) override
    {
        _fixed.setConstant(_values.size(), fixed);
        _num_unfixed = fixed ? 0 : getDimension();
    }

    // implements interface method
    void plusUnfixed(const double* inc) override
    {
        int idx = 0;
        for (int i = 0; i < getDimension(); ++i)
        {
            if (!_fixed(i))
            {
                plus(i, inc[idx]);
                ++idx;
            }
        }
    }

    // implements interface method
    bool hasFixedComponents() const override { return _num_unfixed < getDimension(); }
    // implements interface method
    bool isFixedComponent(int idx) const override { return _fixed[idx]; }

    // implements interface method
    int getNumberFiniteLowerBounds(bool unfixed_only) const override
    {
        if (unfixed_only && _num_unfixed > 0)
        {
            int num = 0;
            for (int i = 0; i < getDimension(); ++i)
            {
                if (!_fixed[i] && _lb[i] > -CORBO_INF_DBL) num += 1;
            }
            return num;
        }
        else
            return (_lb.array() > -CORBO_INF_DBL).count();
    }

    // implements interface method
    int getNumberFiniteUpperBounds(bool unfixed_only) const override
    {
        if (unfixed_only && _num_unfixed > 0)
        {
            int num = 0;
            for (int i = 0; i < getDimension(); ++i)
            {
                if (!_fixed[i] && _ub[i] < CORBO_INF_DBL) num += 1;
            }
            return num;
        }
        else
            return (_ub.array() < CORBO_INF_DBL).count();
    }

    // implements interface method
    int getNumberFiniteBounds(bool unfixed_only) const override
    {
        if (unfixed_only && _num_unfixed > 0)
        {
            int num = 0;
            for (int i = 0; i < getDimension(); ++i)
            {
                if (!_fixed[i] && (_ub[i] < CORBO_INF_DBL || _lb[i] > -CORBO_INF_DBL)) num += 1;
            }
            return num;
        }
        else
            return (_ub.array() < CORBO_INF_DBL || _lb.array() > -CORBO_INF_DBL).count();
    }

    //! Read-only access to the underlying logical array for fixed components
    const Eigen::Array<bool, -1, 1> fixedArray() const { return _fixed; }

 protected:
    Eigen::Array<bool, -1, 1> _fixed;
    int _num_unfixed;
};

}  // namespace corbo

#endif  // SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_VECTOR_VERTEX_H_
