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

#ifndef SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_SCALAR_VERTEX_H_
#define SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_SCALAR_VERTEX_H_

#include <corbo-optimization/hyper_graph/vertex_interface.h>

#include <corbo-core/types.h>

#include <memory>
#include <vector>

namespace corbo {

/**
 * @brief Vertex implementation for scalar values
 *
 * @ingroup optimization hyper-graph
 *
 * This vertex is optimized for scalar (1D) values in
 * contrast to VectorVertex which stores multi-dimensional
 * vector.
 *
 * @see VertexInterface VectorVertex HyperGraph EdgeInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class ScalarVertex : public VertexInterface
{
 public:
    using Ptr  = std::shared_ptr<ScalarVertex>;
    using UPtr = std::unique_ptr<ScalarVertex>;

    //! Default constructor
    ScalarVertex() {}
    //! Construct vertex with given value
    explicit ScalarVertex(double value) : _value(value), _lb(-CORBO_INF_DBL), _ub(CORBO_INF_DBL) {}

    //! Construct vertex with given value and fixed flag
    explicit ScalarVertex(double value, bool fixed) : _value(value), _lb(-CORBO_INF_DBL), _ub(CORBO_INF_DBL), _fixed(fixed) {}

    //! Construct vertex with given value, lower and upper bound
    explicit ScalarVertex(double value, double lb, double ub) : _value(value), _lb(lb), _ub(ub) {}

    //! Construct vertex with given value, lower and upper bound and fixed flag
    explicit ScalarVertex(double value, double lb, double ub, bool fixed) : _value(value), _lb(lb), _ub(ub), _fixed(fixed) {}

    // implements interface method
    int getDimension() const override { return 1; }
    // implements interface method
    int getDimensionUnfixed() const override { return _fixed ? 0 : 1; }

    // implements interface method
    void plus(int /*idx*/, double inc) override { _value += inc; }
    // implements interface method
    void plus(const double* inc) override { _value += *inc; }
    // implements interface method
    void plusUnfixed(const double* inc) override { _value += *inc; }

    // implements interface method
    const double* getData() const override { return &_value; }
    double* getDataRaw() override { return &_value; }

    // implements interface method
    void setData(int /*idx*/, double data) override { _value = data; }

    // directly set all properties
    void set(double value, double lb, double ub, bool fixed)
    {
        _value = value;
        _lb    = lb;
        _ub    = ub;
        _fixed = fixed;
    }

    // implements interface method
    bool hasFixedComponents() const override { return _fixed; }
    // implements interface method
    bool isFixedComponent(int /*idx*/) const override { return _fixed; }

    //! Set lower bound
    void setLowerBound(double lb) { _lb = lb; }
    //! Set upper bound
    void setUpperBound(double ub) { _ub = ub; }
    // implements interface method
    void setLowerBounds(const Eigen::Ref<const Eigen::VectorXd>& lb) override { _lb = lb[0]; }
    // implements interface method
    void setUpperBounds(const Eigen::Ref<const Eigen::VectorXd>& ub) override { _ub = ub[0]; }
    // implements interface method
    void setLowerBound(int /*idx*/, double lb) override { _lb = lb; }
    // implements interface method
    void setUpperBound(int /*idx*/, double ub) override { _ub = ub; }
    // implements interface method
    bool hasFiniteBounds() const override { return _lb > -CORBO_INF_DBL || _ub < CORBO_INF_DBL; }
    // implements interface method
    bool hasFiniteLowerBounds() const override { return _lb > -CORBO_INF_DBL; }
    // implements interface method
    bool hasFiniteUpperBounds() const override { return _ub < CORBO_INF_DBL; }
    // implements interface method
    bool hasFiniteLowerBound(int /*idx*/) const override { return _lb > -CORBO_INF_DBL; }
    // implements interface method
    bool hasFiniteUpperBound(int /*idx*/) const override { return _ub < CORBO_INF_DBL; }
    // implements interface method
    int getNumberFiniteLowerBounds(bool unfixed_only) const override
    {
        if (unfixed_only && _fixed)
            return 0;
        else
            return (int)hasFiniteLowerBounds();
    }
    // implements interface method
    int getNumberFiniteUpperBounds(bool unfixed_only) const override
    {
        if (unfixed_only && _fixed)
            return 0;
        else
            return (int)hasFiniteUpperBounds();
    }
    // implements interface method
    int getNumberFiniteBounds(bool unfixed_only) const override
    {
        if (unfixed_only && _fixed)
            return 0;
        else
            return (int)hasFiniteBounds();
    }

    //! Set vertex (un)fixed
    void setFixed(bool fixed) { _fixed = fixed; }

    // implements interface method
    const double* getLowerBounds() const override { return &_lb; }
    // implements interface method
    const double* getUpperBounds() const override { return &_ub; }

    // Backup values
    // implements interface method
    void push() override { _backup.push_back(_value); }
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
        _value = _backup.back();
    }
    // implements interface method
    void discardTop() override { _backup.pop_back(); }
    // implements interface method
    void clear() override { _backup.clear(); }
    // implements interface method
    int getNumBackups() const override { return (int)_backup.size(); }

    //! Get underlying value
    const double& value() const { return _value; }
    //! Raw access to the underlying value
    double& value() { return _value; }

    //! Get underlying value (this method is for compatibility purposes)
    const double& values() const { return _value; }
    //! Raw access to the underlying value (this method is for compatibility purposes)
    double& values() { return _value; }

 protected:
    double _value;
    double _lb;
    double _ub;

    bool _fixed = false;

    std::vector<double> _backup;
};

}  // namespace corbo

#endif  // SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_SCALAR_VERTEX_H_
