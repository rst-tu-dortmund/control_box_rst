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

#ifndef SRC_CORE_INCLUDE_CORBO_CORE_SIGNALS_H_
#define SRC_CORE_INCLUDE_CORBO_CORE_SIGNALS_H_

#include <corbo-core/console.h>
#include <corbo-core/time.h>
#include <corbo-core/time_series.h>
#include <corbo-core/types.h>

#ifdef MESSAGE_SUPPORT
#include <corbo-communication/messages/core/signals.pb.h>
#endif

#include <map>
#include <memory>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

namespace corbo {

constexpr const char SIGNAL_NAMESPACE_DELIMITER = '/';

//! Available signal types (must match messages::SignalType enumeration)
enum class SignalType {
    TimeSeries,
    Measurement,
    Values,
    TimeSeriesSequence,
    IndexedValues,
    IndexedValuesSet,
    LabeledValues,
    LabeledValuesSet,
    Matrix,
    MatrixSet
};

/**
 * @brief Signal header
 *
 * @ingroup core signals
 *
 * The signal header stores general information about the underlying signal
 * e.g. its unique identifier, the time-stamp and the dimension of the underlying value
 * (the actual interpretation depends on the specific signal type).
 *
 * @see SignalInterface SignalTargetInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
struct SignalHeader
{
    std::string name = "unknown";  //!< Unique signal identifier
    Time time;                     //!< Time-stamp
    int value_dimension  = 0;      //!< Dimension of the underlying signal (type-dependent)
    bool zero_order_hold = false;  //!< If true, two consecutive in the signal are intended to be connected via 0-order hold

    std::string getShortName() const
    {
        std::size_t found = name.find_last_of("/");
        return name.substr(found + 1);
    }

#ifdef MESSAGE_SUPPORT
    //! Export header to message
    void toMessage(corbo::messages::SignalHeader& message) const
    {
        message.set_time(time.toSec());
        message.set_name(name);
        message.set_value_dimension(value_dimension);
        message.set_zero_order_hold(zero_order_hold);
    }
    //! Import header from message
    void fromMessage(const corbo::messages::SignalHeader& message, std::stringstream* issues = nullptr)
    {
        time.fromSec(message.time());
        name            = message.name();
        value_dimension = message.value_dimension();
        zero_order_hold = message.zero_order_hold();
    }
#endif
};

/**
 * @brief Interface class for signals
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
 * Signals are usually sent/streamed to a SignalTargetInterface.
 *
 * @see SignalHeader SignalTargetInterface
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class SignalInterface
{
 public:
    using Ptr      = std::shared_ptr<SignalInterface>;
    using ConstPtr = std::shared_ptr<const SignalInterface>;

    //! Virtual destructor
    virtual ~SignalInterface() {}

    //! Get the signal type according to enumeration SignalType
    virtual SignalType getType() const = 0;

    //! Return labels for the underlying components of the signal (e.g. axes labels)
    virtual void getValueLabels(std::vector<std::string>& sublabels) const {}

    //! The header of the signal
    SignalHeader header;

#ifdef MESSAGE_SUPPORT
    //! Export signal to message
    virtual void toMessage(corbo::messages::Signal& message) const {}
    //! Import signal from message
    virtual void fromMessage(const corbo::messages::Signal& message, std::stringstream* issues = nullptr) {}
#endif
};

/**
 * @brief Measurement signal (value vector recorded at a specific time)
 *
 * @ingroup core signals
 *
 * A measurement is meant to represent a vector of floating point numbers at a specified
 * point in time.
 *
 * @see SignalInterface TimeSeriesSignal
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class Measurement : public SignalInterface
{
 public:
    using Ptr      = std::shared_ptr<Measurement>;
    using ConstPtr = std::shared_ptr<const Measurement>;

    //! Default constructor
    Measurement() {}
    //! Construct measurement with a given time and value vector (Std version)
    Measurement(double time, const std::vector<double>& values) : _time(time), _values(values) {}
    //! Construct measurement with a given time and value vector (Eigen version)
    Measurement(double time, const Eigen::Ref<const Eigen::VectorXd>& values) { set(time, values); }

    // Implements interface method
    SignalType getType() const override { return SignalType::Measurement; }

    //! Set a time and value vector of arbitrary dimension (overrides any existing values)
    void set(double time, const std::vector<double>& values)
    {
        _time   = time;
        _values = values;
    }

    //! Set a time and value vector of arbitrary dimension (overrides any existing values)
    void set(double time, const Eigen::Ref<const Eigen::VectorXd>& values)
    {
        _time = time;
        _values.assign(values.data(), values.data() + values.size());
    }

    //! Access value vector (read-only)
    const std::vector<double>& getValues() const { return _values; }
    //! Access value vector
    std::vector<double>& getValuesRef() { return _values; }

    //! Retrieve recorded time of the measurement
    double getTime() const { return _time; }

    //! Access labels of signal components (might be empty if not provided) [read-only]
    const std::vector<std::string>& getValueLabels() const { return _value_labels; }
    //! Access labels of signal components (might be empty if not provided, allowed to modify directly)
    std::vector<std::string>& getValueLabelsRef() { return _value_labels; }

    // Implements interface method
    void getValueLabels(std::vector<std::string>& sublabels) const override { sublabels = _value_labels; }

#ifdef MESSAGE_SUPPORT
    // Implements interface method
    void toMessage(corbo::messages::Signal& message) const override;
    // Implements interface method
    void fromMessage(const corbo::messages::Signal& message, std::stringstream* issues = nullptr) override;
#endif

 private:
    double _time;                            //!< measurement time
    std::vector<double> _values;             //!< corresponding value vector
    std::vector<std::string> _value_labels;  //!< labels for value vector components (optional)
};

/**
 * @brief Time Series signal (trajectory resp. sequence of values w.r.t. time)
 *
 * @ingroup core signals
 *
 * A time series is meant to represent a sequence of vectors of floating point numbers
 * with individual time stamps. A time series might also be interpreted as a
 * discrete-time trajectory.
 * This class wraps a TimeSeries object into a signal along with a SignalHeader.
 *
 * @see TimeSeries SignalInterface Measurement
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class TimeSeriesSignal : public SignalInterface
{
 public:
    using Ptr      = std::shared_ptr<TimeSeriesSignal>;
    using ConstPtr = std::shared_ptr<const TimeSeriesSignal>;

    //! Default constructor
    TimeSeriesSignal() : _time_series(new TimeSeries) {}
    //! Construct empty time series signal with a dresired value vector dimension
    explicit TimeSeriesSignal(int value_dim) : _time_series(new TimeSeries(value_dim)) {}
    //! Construct time series signal from a time_series object (avoids copying)
    explicit TimeSeriesSignal(TimeSeries::Ptr time_series) : _time_series(time_series) {}
    //! Constructs and initializes time series signal and adds a single measurement
    explicit TimeSeriesSignal(const Measurement& measurment) { add(measurment); }

    //! Check if the time series is empty
    bool isEmpty() const { return _time_series ? _time_series->isEmpty() : true; }

    // Implements interface method
    SignalType getType() const override { return SignalType::TimeSeries; }

    // Implements interface method
    void getValueLabels(std::vector<std::string>& sublabels) const override
    {
        if (_time_series) sublabels = _time_series->getValueLabels();
    }

    //! Add measurement (time and value pair)
    void add(const Measurement& measurement);
    //! Add time and value vector pair (Eigen version)
    void add(double time, const Eigen::Ref<const Eigen::VectorXd>& values);
    //! Add time and value vector pair (STL version)
    void add(double time, const std::vector<double>& values);

    //! Set time series (and override any existing)
    void set(TimeSeries::Ptr time_series) { _time_series = time_series; }

    /**
     * @brief Set time vector and values matrix
     *
     * Any existing values are being replaced.
     * Columns of the values_matrix correspond to vector values at a specific time instance,
     * and hence it is values_matrix: [TimeSeries::getValueDimension() x TimeSeries::getTimeDimension()]
     * @param[in]  time              Vector of all time instances [timeDimension() x 1]
     * @param[in]  values_matrix     Matrix containing all values w.r.t. time [TimeSeries::getValueDimension() x TimeSeries::getTimeDimension()]
     * @param[in]  time_from_start   Specify an offset for all time values for post-processing
     * @return true if setting was successful, otherwise false.
     */
    void set(const Eigen::Ref<const Eigen::VectorXd>& time, const Eigen::Ref<const Eigen::MatrixXd>& values_matrix, double time_from_start = 0.0);

    //! Set time from start (offset to all time stamps in time())
    void setTimeFromStart(double time_from_start);

    //! Read access to the underlying time series object (returns null if not initialized)
    const TimeSeries* getTimeSeries() const { return _time_series.get(); }
    //! Raw access to the underlying time series object (returns null if not initialized)
    TimeSeries* getTimeSeriesRaw() { return _time_series.get(); }

    //! Return shared pointer of the underlying time series (can be empty if not initialized)
    TimeSeries::Ptr getTimeSeriesPtr() const { return _time_series; }

#ifdef MESSAGE_SUPPORT
    // Implements interface method
    void toMessage(corbo::messages::Signal& message) const override;
    // Implements interface method
    void fromMessage(const corbo::messages::Signal& message, std::stringstream* issues = nullptr) override;
#endif

 private:
    TimeSeries::Ptr _time_series;
};

/**
 * @brief Signal for a sequence of time series objects
 *
 * @ingroup core signals
 *
 * This signal stores a sequence of time series objects.
 * Individual time series might be sorted w.r.t different
 * TimeSeries::timeFromStart() values e.g. in order to
 * indicate planned or predicted trajectories at each individual
 * sampling interval.
 * This class wraps a TimeSeriesSequence object into a signal
 * along with a SignalHeader.
 *
 * @see TimeSeriesSequence TimeSeriesSignal TimeSeries SignalInterface Measurement
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class TimeSeriesSequenceSignal : public SignalInterface
{
 public:
    using Ptr      = std::shared_ptr<TimeSeriesSequenceSignal>;
    using ConstPtr = std::shared_ptr<const TimeSeriesSequenceSignal>;

    //! Default constructor
    TimeSeriesSequenceSignal() {}
    //! Construct empty signal with a desired value vector dimension
    explicit TimeSeriesSequenceSignal(int value_dim) : _ts_sequence(new TimeSeriesSequence(value_dim)) {}

    //! Determine if no time series is available
    bool isEmpty() const { return _ts_sequence ? _ts_sequence->isEmpty() : true; }
    // Implements interface method
    SignalType getType() const override { return SignalType::TimeSeriesSequence; }
    // Implements interface method
    void getValueLabels(std::vector<std::string>& sublabels) const override
    {
        // copy values from very first object in the sequence
        if (_ts_sequence && !_ts_sequence->isEmpty() && _ts_sequence->getSequence().front())
            sublabels = _ts_sequence->getSequence().front()->getValueLabels();
    }

    //! Add a new time serie to the sequence
    void add(TimeSeries::Ptr ts);

    //! Set time series sequence (and override any existing)
    void set(TimeSeriesSequence::Ptr ts_sequence) { _ts_sequence = ts_sequence; }

    //! Read access to the underlying time series sequence (returns null if not initialized)
    const TimeSeriesSequence* getSequence() const { return _ts_sequence.get(); }
    //! Raw access to the underlying time series sequence (returns null if not initialized)
    TimeSeriesSequence* getSequenceRaw() { return _ts_sequence.get(); }

    //! Return shared pointer of the underlying time series sequence (can be empty if not initialized)
    TimeSeriesSequence::Ptr getSequencePtr() const { return _ts_sequence; }

#ifdef MESSAGE_SUPPORT
    // Implements interface method
    void toMessage(corbo::messages::Signal& message) const override;
    // Implements interface method
    void fromMessage(const corbo::messages::Signal& message, std::stringstream* issues = nullptr) override;
#endif

 private:
    TimeSeriesSequence::Ptr _ts_sequence;
};

/**
 * @brief Signal containing values indexed by a single integer
 *
 * @ingroup core signals
 *
 * @see IndexedValuesSetSignal SignalInterface Measurement
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class IndexedValuesSignal : public SignalInterface
{
 public:
    using Ptr      = std::shared_ptr<IndexedValuesSignal>;
    using ConstPtr = std::shared_ptr<const IndexedValuesSignal>;

    //! Default constructor
    IndexedValuesSignal() = default;
    //! Construct with desired index
    explicit IndexedValuesSignal(int index) : _index(index) {}
    //! Construct with desired index and single value
    IndexedValuesSignal(int index, double value) : _index(index) { add(value); }
    //! Construct with desired index and value vector
    IndexedValuesSignal(int index, const std::vector<double>& values) : _index(index), _values(values) {}
    //! Construct with desired index and value vector (STL version)
    IndexedValuesSignal(int index, const Eigen::Ref<const Eigen::VectorXd>& values) : _index(index) { add(values); }

    //! Check if the underlying map is empty
    bool isEmpty() const { return _values.empty(); }

    // Implements interface method
    SignalType getType() const override { return SignalType::IndexedValues; }

    // Get internal value dimension
    int getValueDimension() const { return _values.size(); }

    //! Set desired index
    void setIndex(int index) { _index = index; }

    //! Add value
    void add(double value) { _values.push_back(value); }
    //! Add several values
    void add(const Eigen::Ref<const Eigen::VectorXd>& values);
    //! Add several values (STL version)
    void add(const std::vector<double>& values);

    //! Set index value pair (wipes off previous values)
    void set(int index, double value);
    //! Set several values to the desired index (wipes off previous values)
    void set(int index, const Eigen::Ref<const Eigen::VectorXd>& values);
    //! Set several values to the desired index (STL version) (wipes off previous values)
    void set(int index, const std::vector<double>& values);

    //! Return current index
    int getIndex() const { return _index; }
    //! Read access to the underlying values object
    const std::vector<double>& getValues() const { return _values; }
    //! Write access to the underlying values data (use with care)
    std::vector<double>& getValuesRef() { return _values; }

    // Clear value vector
    void clear() { _values.clear(); }

#ifdef MESSAGE_SUPPORT
    void toMessage(corbo::messages::IndexedValues& message) const;
    void fromMessage(const corbo::messages::IndexedValues& message, std::stringstream* issues = nullptr);
    // Implements interface method
    void toMessage(corbo::messages::Signal& message) const override
    {
        header.toMessage(*message.mutable_header());
        toMessage(*message.mutable_indexed_values());
    }
    // Implements interface method
    void fromMessage(const corbo::messages::Signal& message, std::stringstream* issues = nullptr) override
    {
        header.fromMessage(message.header());
        fromMessage(message.indexed_values(), issues);
    }
#endif

 private:
    int _index = 0;
    std::vector<double> _values;
};

/**
 * @brief Signal containing values indexed by an integer (int to double[] map)
 *
 * @ingroup core signals
 *
 * @see IndexedValuesSignal SignalInterface Measurement
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class IndexedValuesSetSignal : public SignalInterface
{
 public:
    using Ptr      = std::shared_ptr<IndexedValuesSetSignal>;
    using ConstPtr = std::shared_ptr<const IndexedValuesSetSignal>;
    using Map      = std::map<int, std::vector<double>>;

    //! Default constructor
    IndexedValuesSetSignal() = default;

    //! Check if the underlying map is empty
    bool isEmpty() const { return _values_map.empty(); }

    // Implements interface method
    SignalType getType() const override { return SignalType::IndexedValuesSet; }

    //! Add index value pair
    void add(int index, double value);
    //! Add from IndexedValuesSignal
    void add(const IndexedValuesSignal& indexed_values);
    //! Add several values to the desired index
    void add(int index, const Eigen::Ref<const Eigen::VectorXd>& values);
    //! Add several values to the desired index (STL version)
    void add(int index, const std::vector<double>& values);

    //! Read access to the underlying map object
    const Map& getData() const { return _values_map; }
    //! Write access to the underlying map (use with care)
    Map& getDataRef() { return _values_map; }

    //! Iterate internal map to find the largest value vector dimension
    int getMaxValueDimension() const;

    void clear() { _values_map.clear(); }

#ifdef MESSAGE_SUPPORT
    void toMessage(corbo::messages::IndexedValuesSet& message) const;
    void fromMessage(const corbo::messages::IndexedValuesSet& message, std::stringstream* issues = nullptr);
    // Implements interface method
    void toMessage(corbo::messages::Signal& message) const override
    {
        header.toMessage(*message.mutable_header());
        toMessage(*message.mutable_indexed_values_set());
    }
    // Implements interface method
    void fromMessage(const corbo::messages::Signal& message, std::stringstream* issues = nullptr) override
    {
        header.fromMessage(message.header());
        fromMessage(message.indexed_values_set(), issues);
    }
#endif

 private:
    Map _values_map;
};

/**
 * @brief Signal containing a simple matrix
 *
 * @ingroup core signals
 *
 * @see MatrixSetSignal IndexedValuesSetSignal SignalInterface Measurement
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class MatrixSignal : public SignalInterface
{
 public:
    using Ptr      = std::shared_ptr<MatrixSignal>;
    using ConstPtr = std::shared_ptr<const MatrixSignal>;

    //! Default constructor
    MatrixSignal() = default;
    //! Construct with matrix
    explicit MatrixSignal(const Eigen::Ref<const Eigen::MatrixXd>& matrix, const std::string& label = "") : _matrix(matrix), _label(label) {}

    //! Check if the underlying map is empty
    bool isEmpty() const { return _matrix.rows() == 0 && _matrix.cols() == 0; }

    // Implements interface method
    SignalType getType() const override { return SignalType::Matrix; }

    int getRowDimension() const { return _matrix.rows(); }
    int getColDimension() const { return _matrix.cols(); }

    //! Set matrix (overwrites internal matrix)
    void set(const Eigen::Ref<const Eigen::MatrixXd>& matrix, const std::string& label = "")
    {
        _matrix = matrix;
        _label  = label;
    }

    //! Read access to the underlying matrix object
    const Eigen::MatrixXd& getMatrix() const { return _matrix; }
    //! Write access to the underlying matrix (use with care)
    Eigen::MatrixXd& getMatrixRef() { return _matrix; }

    //! Read access to the label
    const std::string& getLabel() const { return _label; }
    //! Write access to the label
    std::string& getLabelRef() { return _label; }

#ifdef MESSAGE_SUPPORT
    void toMessage(corbo::messages::Matrix& message) const;
    void fromMessage(const corbo::messages::Matrix& message, std::stringstream* issues = nullptr);
    // Implements interface method
    void toMessage(corbo::messages::Signal& message) const override
    {
        header.toMessage(*message.mutable_header());
        toMessage(*message.mutable_matrix());
    }
    // Implements interface method
    void fromMessage(const corbo::messages::Signal& message, std::stringstream* issues = nullptr) override
    {
        header.fromMessage(message.header());
        fromMessage(message.matrix(), issues);
    }
#endif

 private:
    Eigen::MatrixXd _matrix;
    std::string _label;
};

/**
 * @brief Signal containing a set of matrices
 *
 * @ingroup core signals
 *
 * @see MatrixSignal IndexedValuesSignal SignalInterface Measurement
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class MatrixSetSignal : public SignalInterface
{
 public:
    using Ptr      = std::shared_ptr<MatrixSetSignal>;
    using ConstPtr = std::shared_ptr<const MatrixSetSignal>;
    using Map      = std::map<int, std::vector<double>>;

    //! Default constructor
    MatrixSetSignal() = default;

    //! Check if the underlying map is empty
    bool isEmpty() const { return _matrix_set.empty(); }

    // Implements interface method
    SignalType getType() const override { return SignalType::MatrixSet; }

    //! Add matrix signal
    void add(MatrixSignal::Ptr& matrix_signal);
    //! Add matrix from Eigen type
    void add(const Eigen::Ref<const Eigen::MatrixXd>& matrix, const std::string& label = "");

    //! Read access to the underlying map object
    const std::vector<MatrixSignal::Ptr>& getData() const { return _matrix_set; }
    //! Write access to the underlying map (use with care)
    std::vector<MatrixSignal::Ptr>& getDataRef() { return _matrix_set; }

    void clear() { _matrix_set.clear(); }

#ifdef MESSAGE_SUPPORT
    void toMessage(corbo::messages::MatrixSet& message) const;
    void fromMessage(const corbo::messages::MatrixSet& message, std::stringstream* issues = nullptr);
    // Implements interface method
    void toMessage(corbo::messages::Signal& message) const override
    {
        header.toMessage(*message.mutable_header());
        toMessage(*message.mutable_matrix_set());
    }
    // Implements interface method
    void fromMessage(const corbo::messages::Signal& message, std::stringstream* issues = nullptr) override
    {
        header.fromMessage(message.header());
        fromMessage(message.matrix_set(), issues);
    }
#endif

 private:
    std::vector<MatrixSignal::Ptr> _matrix_set;
};

}  // namespace corbo

#endif  // SRC_CORE_INCLUDE_CORBO_CORE_SIGNALS_H_
