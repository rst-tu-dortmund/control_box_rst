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

#ifndef SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_EDGE_CACHE_H_
#define SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_EDGE_CACHE_H_

#include <corbo-core/types.h>

#include <Eigen/Core>

#include <corbo-core/console.h>

#include <vector>

namespace corbo {

// TODO(roesmann): we might change the cache to default construct the eigen vectors to avoid memory allocation during caching.
//                 we could think of using a container with id access rather than a stack as underlying datatype
class EdgeCache
{
 public:
    enum Recent : int { Current = 0, Previous = 1 };

    void reserveMemoryValues(int num_value_vectors) { _values.reserve(num_value_vectors); }
    void reserveMemoryJacobians(int num_jacobians) { _jacobians.reserve(num_jacobians); }

    Eigen::VectorXd& pushValues(int value_dim)
    {
        PRINT_DEBUG_COND_ONCE(_values.size() >= _values.capacity(),
                              "EdgeCache::pushValues(): cache capacity reached; you might better reserve more space in advance.");
#if __cplusplus > 201402L
        return _values.emplace_back();
#else
        _values.emplace_back(value_dim);
        return _values.back();
#endif
    }

    void popValues() { _values.pop_back(); }
    Eigen::VectorXd& topValues()
    {
        assert(!_values.empty());
        return _values.back();
    }
    Eigen::VectorXd& recentValues(int reverse_idx)
    {
        assert(reverse_idx < _values.size());
        return *(&_values.back() - reverse_idx);
    }
    int sizeValues() { return (int)_values.size(); }
    void clearValues() { _values.clear(); }
    Eigen::MatrixXd& pushJacobian(int value_dim, int param_dim)
    {
        PRINT_DEBUG_COND_ONCE(_jacobians.size() >= _values.capacity(),
                              "EdgeCache::pushJacobian(): cache capacity reached; you might better reserve more space in advance.");
#if __cplusplus > 201402L
        return _values.emplace_back();
#else
        _jacobians.emplace_back(value_dim, param_dim);
        return _jacobians.back();
#endif
    }

    void popJacobians() { _jacobians.pop_back(); }
    Eigen::MatrixXd& topJacobians()
    {
        assert(!_jacobians.empty());
        return _jacobians.back();
    }
    Eigen::MatrixXd& recentJacobians(int reverse_idx)
    {
        assert(reverse_idx < _jacobians.size());
        return *(&_jacobians.back() - reverse_idx);
    }
    int sizeJacobians() { return (int)_jacobians.size(); }
    void clearJacobians() { _jacobians.clear(); }

    void clear()
    {
        clearValues();
        clearJacobians();
    }

    bool getCustomFlag() const { return _custom_flag; }
    void setCustomFlag(bool flag) { _custom_flag = flag; }

 protected:
    std::vector<Eigen::VectorXd> _values;
    std::vector<Eigen::MatrixXd> _jacobians;

    bool _custom_flag = false;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

}  // namespace corbo

#endif  // SRC_OPTIMIZATION_INCLUDE_CORBO_OPTIMIZATION_HYPER_GRAPH_EDGE_CACHE_H_
