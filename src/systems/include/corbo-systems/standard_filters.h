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

#ifndef SRC_SYSTEMS_INCLUDE_CORBO_SYSTEMS_STANDARD_FILTERS_H_
#define SRC_SYSTEMS_INCLUDE_CORBO_SYSTEMS_STANDARD_FILTERS_H_

#include <corbo-systems/filter_interface.h>

#ifdef MESSAGE_SUPPORT
#include <corbo-communication/messages/systems/filters.pb.h>
#endif

#include <deque>
#include <memory>

namespace corbo {

/**
 * @brief Moving Average Filter
 *
 * @ingroup systems
 *
 * This class implements a moving average filter based on a specified window size.
 *
 * @remark This interface is provided with factory support (FilterFactory).
 *
 * @see FilterInterface MovingMedianFilter MovingLeastSquaresFilter
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class MovingAverageFilter : public FilterInterface
{
 public:
    using Ptr = std::shared_ptr<MovingAverageFilter>;

    // implements interface method
    FilterInterface::Ptr getInstance() const override { return std::make_shared<MovingAverageFilter>(); }

    // implements interface method
    double filter(double t, double value) override;

    // implements interface method
    void reset() override;

    void setWindowSize(int window_size) { _window_size = window_size; }

#ifdef MESSAGE_SUPPORT
    //! Export to message
    virtual void toMessage(corbo::messages::MovingAverageFilter& message) const;
    //! Import from message
    virtual void fromMessage(const corbo::messages::MovingAverageFilter& message, std::stringstream* issues = nullptr);

    // implements interface method
    void toMessage(corbo::messages::Filter& message) const override { toMessage(*message.mutable_moving_average()); }
    // implements interface method
    void fromMessage(const corbo::messages::Filter& message, std::stringstream* issues = nullptr) override
    {
        fromMessage(message.moving_average(), issues);
    }
#endif

 private:
    int _window_size = 5;

    std::deque<double> _values;
    bool _first_value = true;
};

FACTORY_REGISTER_FILTER(MovingAverageFilter)

/**
 * @brief Moving Median Filter
 *
 * @ingroup systems
 *
 * This class implements a moving median filter based on a specified window size.
 *
 * @remark This interface is provided with factory support (FilterFactory).
 *
 * @see FilterInterface MovingAverageFilter MovingLeastSquaresFilter
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class MovingMedianFilter : public FilterInterface
{
 public:
    using Ptr = std::shared_ptr<MovingMedianFilter>;

    // implements interface method
    FilterInterface::Ptr getInstance() const override { return std::make_shared<MovingMedianFilter>(); }

    // implements interface method
    double filter(double t, double value) override;

    // implements interface method
    void reset() override;

    void setWindowSize(int window_size) { _window_size = window_size; }

#ifdef MESSAGE_SUPPORT
    //! Export to message
    virtual void toMessage(corbo::messages::MovingMedianFilter& message) const;
    //! Import from message
    virtual void fromMessage(const corbo::messages::MovingMedianFilter& message, std::stringstream* issues = nullptr);

    // implements interface method
    void toMessage(corbo::messages::Filter& message) const override { toMessage(*message.mutable_moving_median()); }
    // implements interface method
    void fromMessage(const corbo::messages::Filter& message, std::stringstream* issues = nullptr) override
    {
        fromMessage(message.moving_median(), issues);
    }
#endif

 private:
    int _window_size = 5;

    std::deque<double> _values;
    bool _first_value = true;
};

FACTORY_REGISTER_FILTER(MovingMedianFilter)

/**
 * @brief Moving LeastSquares Filter
 *
 * @ingroup systems
 *
 * This class implements a moving filter which solves an unconstrained
 * least squares problem along the specified window size.
 * It returns the last point of the regression result.
 *
 * @remark This interface is provided with factory support (FilterFactory).
 *
 * @see FilterInterface MovingAverageFilter MovingMedianFilter
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class MovingLeastSquaresFilter : public FilterInterface
{
 public:
    using Ptr = std::shared_ptr<MovingMedianFilter>;

    // implements interface method
    FilterInterface::Ptr getInstance() const override { return std::make_shared<MovingLeastSquaresFilter>(); }

    // implements interface method
    double filter(double t, double value) override;

    // implements interface method
    void reset() override;

    void setWindowSize(int window_size) { _window_size = window_size; }
    void setPolynomialOrder(int order) { _order = order; }

#ifdef MESSAGE_SUPPORT
    //! Export to message
    virtual void toMessage(corbo::messages::MovingLeastSquaresFilter& message) const;
    //! Import from message
    virtual void fromMessage(const corbo::messages::MovingLeastSquaresFilter& message, std::stringstream* issues = nullptr);

    // implements interface method
    void toMessage(corbo::messages::Filter& message) const override { toMessage(*message.mutable_moving_least_squares()); }
    // implements interface method
    void fromMessage(const corbo::messages::Filter& message, std::stringstream* issues = nullptr) override
    {
        fromMessage(message.moving_least_squares(), issues);
    }
#endif

 protected:
    void allocateMemory();
    void buildRegressionMatrix();

 private:
    // parameters
    int _window_size = 5;
    int _order       = 3;

    // internal states
    Eigen::VectorXd _values;
    Eigen::VectorXd _time;
    Eigen::MatrixXd _regression_matrix;
    int _num_values = 0;
};

FACTORY_REGISTER_FILTER(MovingLeastSquaresFilter)

}  // namespace corbo

#endif  // SRC_SYSTEMS_INCLUDE_CORBO_SYSTEMS_STANDARD_FILTERS_H_
