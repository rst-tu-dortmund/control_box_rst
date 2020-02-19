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

#ifndef SRC_CORE_INCLUDE_CORBO_CORE_SIGNAL_TARGET_INTERFACE_H_
#define SRC_CORE_INCLUDE_CORBO_CORE_SIGNAL_TARGET_INTERFACE_H_

#include <corbo-core/signals.h>

#include <memory>
#include <string>
#include <vector>

namespace corbo {

/**
 * @brief Interface class for signal targets
 *
 * @ingroup core signals
 *
 * Signals are generated in certain modules like
 * controllers, obervers or tasks in order to provide internal data
 * (e.g. measurements, states, statistics, ...).
 * Signals are usually derived from SignalInterface which implements
 * an adapter to the actually data including general signal information
 * (refer to SignalHeader).
 *
 * Signals are identified using unique names. Signalnamespaces
 * can be set using a "/" as delimiter, e.g.:
 * /c this/is/a_signal/with/namespaces.
 * Names can only be composed of alphanumeric characters and underscores (_).
 *
 * In general signals are provided as data stream and hence multiple signals with
 * different time stamps can share the same identifier (name).
 * Targets need make sure to store or combine multiple signals of the same identifier
 * appropriately.
 *
 * @see SignalTargetInterface SignalTargetRPC YamlExporter
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class SignalTargetInterface
{
 public:
    using Ptr = std::shared_ptr<SignalTargetInterface>;

    //! Virtual destructor
    virtual ~SignalTargetInterface() {}

    /**
     * @brief Register a measurement type at current target
     *
     * A measurement is meant to represent a vector of floating point numbers at a specified
     * point in time.
     * Hence, targets might store measurements of the same type (unique_name) into a TimeSeries object
     * in order to record the whole sequence for batch-processing or visualization.
     *
     * @remarks Signal registration is optional (at least for the current set of targets), but it is
     *          recommended for pre-memory-allocation purposes and e.g. for online plotting support in the GUI.
     *
     * @see Measurement sendMeasurement
     *
     * @param[in]  unique_name       Identifier for the current streams of measurements (same type, but different time-stamps allowed).
     * @param[in]  value_dimension   Dimension of the measurement vector (number of floating point values).
     * @param[in]  value_labels      Labels for individual components of the measurement vector, e.g. for data export or plotting (axes
     * labels).
     * @param[in]  zero_order_hold   Specify whether the signal is intended to be holded (zero order hold) between two consecutive values.
     */
    virtual void registerMeasurement(const std::string& unique_name, int value_dimension, const std::vector<std::string>& value_labels = {},
                                     bool zero_order_hold = false) = 0;

    /**
     * @brief Register a time series type at current target
     *
     * A TimeSeries is meant to represent a sequence of value vectors or measurements.
     *
     * The channel associated with \c unique_name might be used to send a single TimeSeriesSignal
     * or a stream of TimeSeriesSignal objects.
     * Hence, targets might store signals of the same type (unique_name) into a collection of TimeSeries objects
     * in order to record the whole sequence for batch-processing or visualization.
     *
     * @remark Make sure to set TimeSeries::setTimeFromStart() in order to specify a valid sequence of TimeSeries objects
     *         w.r.t. time. The SignalHeader time stamp is currently not in use for this purpose.
     *
     * @remark Signal registration is optional (at least for the current set of targets), but it is
     *         recommended for pre-memory-allocation purposes and e.g. for online plotting support in the GUI.
     *
     * @see TimeSeriesSignal TimeSeries sendTimeSeries
     *
     * @param[in]  unique_name       Identifier for the current streams of measurements (same type, but different time-stamps allowed).
     * @param[in]  value_dimension   Dimension of the measurement vector (number of floating point values).
     * @param[in]  zero_order_hold   Specify whether the signal is intended to be holded (zero order hold) between two consecutive values.
     */
    virtual void registerTimeSeries(const std::string& unique_name, int value_dimension, bool zero_order_hold = false) = 0;

    /**
     * @brief Send a measurement to the target
     *
     * Refer to registerMeasurement() for a detailed description.
     *
     * @warning Make sure to provide a valid SignalHeader including a unique name
     *          or use on of the method overloads for a more convenient usage.
     *
     * @see Measurement registerMeasurement
     *
     * @param[in]  measurement     Measurement object with valid header to be sent.
     */
    virtual void sendMeasurement(Measurement::ConstPtr measurement) = 0;

    /**
     * @brief Send a measurement to the target
     *
     * Refer to registerMeasurement() for a detailed description.
     * This method overloads sendMeasurement(Measurement::Ptr measurement)
     * and constructs a Measurement object from a given unique_name,
     * time and value vector (provided as std::vector<double> type).
     *
     * Note, the time corresponds to the value vector and not to the SignalHeader time-stamp.
     * The SignalHeader time-stamp is set to the current time.
     *
     * @see Measurement registerMeasurement
     *
     * @param[in]  unique_name     Identifier for the current streams of measurements (same type, but different time-stamps allowed).
     * @param[in]  time            Time at which the values are recorded (usually in seconds).
     * @param[in]  values          Value vector containing the actual recorded data.
     */
    virtual void sendMeasurement(const std::string& unique_name, double time, const std::vector<double>& values,
                                 const std::vector<std::string>& value_labels = {})
    {
        Measurement::Ptr measurement        = std::make_shared<Measurement>(time, values);
        measurement->header.time            = Time::now();
        measurement->header.name            = unique_name;
        measurement->header.value_dimension = values.size();
        measurement->getValueLabelsRef()    = value_labels;
        sendMeasurement(measurement);
    }

    /**
     * @brief Send a measurement to the target
     *
     * Refer to registerMeasurement() for a detailed description.
     * This method overloads sendMeasurement(Measurement::Ptr measurement)
     * and constructs a Measurement object from a given unique_name,
     * time and value vector (provided as Eigen::VectorXd type).
     *
     * Note, the time corresponds to the value vector and not to the SignalHeader time-stamp.
     * The SignalHeader time-stamp is set to the current time.
     *
     * @see Measurement registerMeasurement
     *
     * @param[in]  unique_name     Identifier for the current streams of measurements (same type, but different time-stamps allowed).
     * @param[in]  time            Time at which the values are recorded (usually in seconds).
     * @param[in]  values          Value vector containing the actual recorded data.
     */
    virtual void sendMeasurement(const std::string& unique_name, double time, const Eigen::Ref<const Eigen::VectorXd>& values,
                                 const std::vector<std::string>& value_labels = {})
    {
        Measurement::Ptr measurement        = std::make_shared<Measurement>(time, values);
        measurement->header.time            = Time::now();
        measurement->header.name            = unique_name;
        measurement->header.value_dimension = values.rows();
        measurement->getValueLabelsRef()    = value_labels;
        sendMeasurement(measurement);
    }

    /**
     * @brief Send a time series to the target
     *
     * Refer to registerTimeSeries() for a detailed description.
     *
     * @warning Make sure to provide a valid SignalHeader including a unique name
     *          or use on of the method overloads for a more convenient usage.
     *
     * @see TimeSeriesSignal TimeSeries registerTimeSeries
     *
     * @param[in]  time_series   TimeSeriesSignal object with valid header to be sent.
     */
    virtual void sendTimeSeries(TimeSeriesSignal::Ptr time_series) = 0;

    /**
     * @brief Send a time series to the target
     *
     * Refer to registerTimeSeries() for a detailed description.
     * This method overloads sendTimeSeries(TimeSeriesSignal::Ptr time_series)
     * and constructs a TimeSeriesSignal object from a given unique_name,
     * and TimeSeries object.
     *
     * @see TimeSeriesSignal TimeSeries registerTimeSeries
     *
     * @todo Make sure that everything is thread-safe if we forward a shared_ptr
     *
     * @param[in]  unique_name     Identifier for the current stream of TimeSeriesSignal objects
     * @param[in]  time_series     Shared TimeSeries object
     */
    virtual void sendTimeSeries(const std::string& unique_name, TimeSeries::Ptr time_series)
    {
        TimeSeriesSignal::Ptr ts_signal   = std::make_shared<TimeSeriesSignal>();
        ts_signal->header.time            = Time::now();
        ts_signal->header.name            = unique_name;
        ts_signal->header.value_dimension = time_series->getValueDimension();
        ts_signal->set(time_series);
        sendTimeSeries(ts_signal);
    }

    /**
     * @brief Send signal containing values indexed by a single integer
     *
     * @warning Make sure to provide a valid SignalHeader including a unique name
     *          or use on of the method overloads for a more convenient usage.
     *
     * @see IndexedValuesSignal IndexedValuesSetSignal
     *
     * @param[in]  indexed_values   IndexedValuesSignal object with valid header to be sent.
     */
    virtual void sendIndexedValues(IndexedValuesSignal::Ptr indexed_values) = 0;

    /**
     * @brief Send signal containing values indexed by a single integer
     *
     * This method overloads sendIndexedValues(IndexedValuesSignal::Ptr indexed_values)
     * and constructs a IndexedValuesSignal object from a given unique_name,
     * index and value vector (provided as std::vector<double> type).
     *
     * The SignalHeader time-stamp is set to the current time.
     *
     * @see IndexedValuesSignal IndexedValuesSetSignal
     *
     * @param[in]  unique_name     Identifier for the current streams of values (same type, but different time-stamps allowed).
     * @param[in]  index           Index as (non-unique) identifier for the current value vector (values with same identifier will be stacked)
     * @param[in]  values          Value vector (STL version)
     */
    virtual void sendIndexedValues(const std::string& unique_name, int index, const std::vector<double>& values)
    {
        IndexedValuesSignal::Ptr indexed_values = std::make_shared<IndexedValuesSignal>(index, values);
        indexed_values->header.time             = Time::now();
        indexed_values->header.name             = unique_name;
        indexed_values->header.value_dimension  = 1;  // for indexed values just one, since we won't have a fixed dimension, just statistical data
        sendIndexedValues(indexed_values);
    }

    /**
     * @brief Send signal containing values indexed by a single integer
     *
     * This method overloads sendIndexedValues(IndexedValuesSignal::Ptr indexed_values)
     * and constructs a IndexedValuesSignal object from a given unique_name,
     * index and value vector (provided as Eigen::VectorXd type).
     *
     * The SignalHeader time-stamp is set to the current time.
     *
     * @see IndexedValuesSignal IndexedValuesSetSignal
     *
     * @param[in]  unique_name     Identifier for the current streams of values (same type, but different time-stamps allowed).
     * @param[in]  index           Index as (non-unique) identifier for the current value vector (values with same identifier will be stacked)
     * @param[in]  values          Value vector (Eigen::VectorXd version)
     */
    virtual void sendIndexedValues(const std::string& unique_name, int index, const Eigen::Ref<const Eigen::VectorXd>& values)
    {
        IndexedValuesSignal::Ptr indexed_values = std::make_shared<IndexedValuesSignal>(index, values);
        indexed_values->header.time             = Time::now();
        indexed_values->header.name             = unique_name;
        indexed_values->header.value_dimension  = 1;  // for indexed values just one, since we won't have a fixed dimension, just statistical data
        sendIndexedValues(indexed_values);
    }

    /**
     * @brief Send signal containing a set of values indexed by integers (int to double[] map)
     *
     * @warning Make sure to provide a valid SignalHeader including a unique name
     *          or use on of the method overloads for a more convenient usage.
     *
     * @see IndexedValuesSetSignal IndexedValuesSignal
     *
     * @param[in]  indexed_values_set  IndexedValuesSetSignal object with valid header to be sent.
     */
    virtual void sendIndexedValuesSet(IndexedValuesSetSignal::Ptr indexed_values_set) = 0;

    /**
     * @brief Send a matrix to the target
     *
     * @warning Make sure to provide a valid SignalHeader including a unique name
     *          or use on of the method overloads for a more convenient usage.
     *
     * @see TimeSeriesSignal TimeSeries registerTimeSeries
     *
     * @param[in]  time_series   TimeSeriesSignal object with valid header to be sent.
     */
    virtual void sendMatrix(MatrixSignal::Ptr matrix) = 0;

    /**
     * @brief Send a matrix to the target
     *
     * This method overloads sendMatrix(MatrixSignal::Ptr matrix)
     * and constructs a MatrixSignal object from a given unique_name and matrix (provided as Eigen::MatrixXd type).
     * Each matrix can be associated with an additional label.
     *
     * The SignalHeader time-stamp is set to the current time.
     *
     * @see MatrixSignal MatrixSetSignal
     *
     * @param[in]  unique_name     Identifier for the current streams of matrices (same type, but different time-stamps allowed).
     * @param[in]  matrix          Matrix with actual values (Eigen::MatrixXd version)
     * @param[in]  label           Optional label which describes the content of the matrix [optional]
     */
    virtual void sendMatrix(const std::string& unique_name, const Eigen::Ref<const Eigen::MatrixXd>& matrix, const std::string& label = "")
    {
        MatrixSignal::Ptr matrix_sig       = std::make_shared<MatrixSignal>(matrix, label);
        matrix_sig->header.time            = Time::now();
        matrix_sig->header.name            = unique_name;
        matrix_sig->header.value_dimension = 1;  // not relevant for this type yet (TODO(roesmann): unify)
        sendMatrix(matrix_sig);
    }
};

}  // namespace corbo

#endif  // SRC_CORE_INCLUDE_CORBO_CORE_SIGNAL_TARGET_INTERFACE_H_
