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

#ifndef SRC_SYSTEMS_INCLUDE_CORBO_SYSTEMS_TIME_VALUE_BUFFER_H_
#define SRC_SYSTEMS_INCLUDE_CORBO_SYSTEMS_TIME_VALUE_BUFFER_H_

#include <Eigen/Core>

#include <corbo-core/console.h>

#include <limits>
#include <memory>
#include <vector>

namespace corbo {

/**
 * @brief Time Delay Object for Piecewise-Constant Signals
 *
 * @ingroup systems
 *
 * This class realzies a time delay for piecewise-constant signals
 * (discretization can be arbitrarily and must not be uniform).
 * Internally, the object caches the last few controls in combination
 * with their time durations.
 * This time is then mapped to the delayed time basis.
 * In case no previous signal value is present, a previously defined
 * initial value is returned (default: zero vector).
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class TimeValueBuffer
{
 public:
    //! Default constructor
    TimeValueBuffer() = default;

    /**
     * @brief Compute the delayed values
     *
     * Adds the input u to the internal cache and returns all vectors u
     * that fall into the time interval [t, dt].
     * If no previously cache vector is available (e.g. at the beginning),
     * an initial value is returned (see setInitialValue()).
     *
     * @see setInitialValue(), setTimeDelay(), reset()
     *
     * @param[in]   t           Current time stamp
     * @param[in]   u           Value which is expected at \c t in case of ignoring the time-delay
     * @param[in]   dt          Duration for \c u (note, we assume a piecewise constant signal)
     * @param[out]  useq_out    vector containing pairs of (dt and controls) (from the delayed time-base) [useq_out is cleared]
     */
    void getValues(double ts, double dt, std::vector<std::pair<double, Eigen::VectorXd>>& useq_out);

    void appendValues(double t, const Eigen::Ref<const Eigen::VectorXd>& u);

    //! Specify initial value
    void setInitialValue(const Eigen::Ref<const Eigen::VectorXd>& uinit) { _uinit = uinit; }

    bool isEmpty() const { return _ucache.empty(); }

    //! Reset internal buffer
    void reset() { _ucache.clear(); }

 protected:
 private:
    Eigen::VectorXd _uinit;

    std::vector<std::pair<double, Eigen::VectorXd>> _ucache;
};

}  // namespace corbo

#endif  // SRC_SYSTEMS_INCLUDE_CORBO_SYSTEMS_TIME_VALUE_BUFFER_H_
