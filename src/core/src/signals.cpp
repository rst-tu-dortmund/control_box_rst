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

#include <corbo-core/signals.h>

namespace corbo {

// void ValuesSignal::toMessage(corbo::messages::Signal& message) const
//{
//    header.toMessage(*message.mutable_header());
//    google::protobuf::RepeatedField<double> values(_values.begin(), _values.end());
//    message.mutable_values()->mutable_values()->Swap(&values);
//}
// void ValuesSignal::fromMessage(const corbo::messages::Signal& message, std::stringstream* issues)
//{
//    header.fromMessage(message.header(), issues);
//    _values.assign(message.values().values().begin(), message.values().values().end());
//}

#ifdef MESSAGE_SUPPORT
void Measurement::toMessage(corbo::messages::Signal& message) const
{
    header.toMessage(*message.mutable_header());
    message.mutable_measurement()->set_time(_time);
    google::protobuf::RepeatedField<double> values(_values.begin(), _values.end());
    message.mutable_measurement()->mutable_values()->Swap(&values);
    // labels
    message.mutable_measurement()->clear_value_labels();
    for (const std::string& label : _value_labels) message.mutable_measurement()->add_value_labels(label);
}
void Measurement::fromMessage(const corbo::messages::Signal& message, std::stringstream* issues)
{
    header.fromMessage(message.header(), issues);
    _time = message.measurement().time();
    _values.assign(message.measurement().values().begin(), message.measurement().values().end());
    // labels
    _value_labels.clear();
    for (int i = 0; i < message.measurement().value_labels_size(); ++i) _value_labels.push_back(message.measurement().value_labels(i));
}
#endif

void TimeSeriesSignal::add(const Measurement& measurement)
{
    if (!_time_series) _time_series = std::make_shared<TimeSeries>();
    _time_series->add(measurement.getTime(), measurement.getValues());

    if (!measurement.getValueLabels().empty() && _time_series->getValueLabels().empty())
    {
        _time_series->getValueLabelsRef() = measurement.getValueLabels();
    }
}

void TimeSeriesSignal::add(double time, const Eigen::Ref<const Eigen::VectorXd>& values)
{
    if (!_time_series) _time_series = std::make_shared<TimeSeries>();
    _time_series->add(time, values);
}

void TimeSeriesSignal::add(double time, const std::vector<double>& values)
{
    if (!_time_series) _time_series = std::make_shared<TimeSeries>();
    _time_series->add(time, values);
}

void TimeSeriesSignal::set(const Eigen::Ref<const Eigen::VectorXd>& time, const Eigen::Ref<const Eigen::MatrixXd>& values_matrix,
                           double time_from_start)
{
    if (!_time_series) _time_series = std::make_shared<TimeSeries>();
    _time_series->set(time, values_matrix, time_from_start);
}

void TimeSeriesSignal::setTimeFromStart(double time_from_start)
{
    if (!_time_series) _time_series = std::make_shared<TimeSeries>();
    _time_series->setTimeFromStart(time_from_start);
}

#ifdef MESSAGE_SUPPORT
void TimeSeriesSignal::toMessage(corbo::messages::Signal& message) const
{
    header.toMessage(*message.mutable_header());
    if (_time_series) _time_series->toMessage(*message.mutable_time_series());
}

void TimeSeriesSignal::fromMessage(const corbo::messages::Signal& message, std::stringstream* issues)
{
    header.fromMessage(message.header());
    if (!_time_series) _time_series = std::make_shared<TimeSeries>();
    _time_series->fromMessage(message.time_series(), issues);
}
#endif

void TimeSeriesSequenceSignal::add(TimeSeries::Ptr ts)
{
    if (!_ts_sequence) _ts_sequence = std::make_shared<TimeSeriesSequence>();
    _ts_sequence->add(ts);
}

#ifdef MESSAGE_SUPPORT
void TimeSeriesSequenceSignal::toMessage(corbo::messages::Signal& message) const
{
    header.toMessage(*message.mutable_header());
    if (_ts_sequence) _ts_sequence->toMessage(*message.mutable_ts_sequence());
}
void TimeSeriesSequenceSignal::fromMessage(const corbo::messages::Signal& message, std::stringstream* issues)
{
    header.fromMessage(message.header());
    if (!_ts_sequence) _ts_sequence = std::make_shared<TimeSeriesSequence>();
    _ts_sequence->fromMessage(message.ts_sequence(), issues);
}
#endif

void IndexedValuesSignal::add(const Eigen::Ref<const Eigen::VectorXd>& values)
{
    std::copy(values.data(), values.data() + values.size(), std::back_inserter(_values));
}
void IndexedValuesSignal::add(const std::vector<double>& values) { std::copy(values.begin(), values.end(), std::back_inserter(_values)); }

void IndexedValuesSignal::set(int index, double value)
{
    clear();
    setIndex(index);
    add(value);
}
void IndexedValuesSignal::set(int index, const Eigen::Ref<const Eigen::VectorXd>& values)
{
    clear();
    setIndex(index);
    add(values);
}
void IndexedValuesSignal::set(int index, const std::vector<double>& values)
{
    clear();
    setIndex(index);
    add(values);
}

#ifdef MESSAGE_SUPPORT
void IndexedValuesSignal::toMessage(messages::IndexedValues& message) const
{
    message.set_index(_index);
    google::protobuf::RepeatedField<double> values(_values.begin(), _values.end());
    message.mutable_values()->Swap(&values);
}

void IndexedValuesSignal::fromMessage(const messages::IndexedValues& message, std::stringstream* issues)
{
    _index = message.index();
    _values.assign(message.values().begin(), message.values().end());
}
#endif

void IndexedValuesSetSignal::add(int index, double value) { _values_map[index].push_back(value); }

void IndexedValuesSetSignal::add(const IndexedValuesSignal& indexed_values)
{
    std::copy(indexed_values.getValues().begin(), indexed_values.getValues().end(), std::back_inserter(_values_map[indexed_values.getIndex()]));
}

void IndexedValuesSetSignal::add(int index, const Eigen::Ref<const Eigen::VectorXd>& values)
{
    std::copy(values.data(), values.data() + values.size(), std::back_inserter(_values_map[index]));
}

void IndexedValuesSetSignal::add(int index, const std::vector<double>& values)
{
    std::copy(values.begin(), values.end(), std::back_inserter(_values_map[index]));
}

int IndexedValuesSetSignal::getMaxValueDimension() const
{
    int max_dim = 0;
    for (const auto& elem : _values_map)
    {
        if ((int)elem.second.size() > max_dim) max_dim = (int)elem.second.size();
    }
    return max_dim;
}

#ifdef MESSAGE_SUPPORT
void IndexedValuesSetSignal::toMessage(messages::IndexedValuesSet& message) const
{
    for (auto& item : _values_map)
    {
        messages::Vector& value_vec = (*message.mutable_indexed_values())[item.first];
        for (const double& value : item.second) value_vec.add_values(value);
    }
}

void IndexedValuesSetSignal::fromMessage(const messages::IndexedValuesSet& message, std::stringstream* issues)
{
    clear();

    for (const auto& item : message.indexed_values())
    {
        std::vector<double>& value_vec = _values_map[item.first];
        value_vec.assign(item.second.values().begin(), item.second.values().end());
    }
}
#endif

#ifdef MESSAGE_SUPPORT
void MatrixSignal::toMessage(messages::Matrix& message) const
{
    message.set_rows(_matrix.rows());
    message.set_cols(_matrix.cols());
    message.set_row_major(true);

    message.mutable_data()->Resize(_matrix.rows() * _matrix.cols(), 0);  // TODO(roesmann): inefficient, to first fill with zeros!
    Eigen::Map<Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(message.mutable_data()->mutable_data(), _matrix.rows(), _matrix.cols()) = _matrix;

    message.set_label(_label);
}

void MatrixSignal::fromMessage(const messages::Matrix& message, std::stringstream* isses)
{
    if (message.row_major())
    {
        _matrix = Eigen::Map<const Eigen::Matrix<double, -1, -1, Eigen::RowMajor>>(message.data().data(), message.rows(), message.cols());
    }
    else
    {
        _matrix = Eigen::Map<const Eigen::Matrix<double, -1, -1, Eigen::ColMajor>>(message.data().data(), message.rows(), message.cols());
    }

    _label = message.label();
}

#endif

void MatrixSetSignal::add(MatrixSignal::Ptr& matrix_signal) { _matrix_set.push_back(matrix_signal); }

void MatrixSetSignal::add(const Eigen::Ref<const Eigen::MatrixXd>& matrix, const std::string& label)
{
    MatrixSignal::Ptr matrix_signal = std::make_shared<MatrixSignal>(matrix, label);
    _matrix_set.push_back(matrix_signal);
}

#ifdef MESSAGE_SUPPORT
void MatrixSetSignal::toMessage(messages::MatrixSet& message) const
{
    for (const MatrixSignal::Ptr& mat : _matrix_set)
    {
        messages::Matrix* mat_msg = message.add_matrices();
        mat->toMessage(*mat_msg);
    }
}

void MatrixSetSignal::fromMessage(const messages::MatrixSet& message, std::stringstream* issues)
{
    clear();
    for (int i = 0; i < message.matrices_size(); ++i)
    {
        MatrixSignal::Ptr matrix_signal = std::make_shared<MatrixSignal>();
        matrix_signal->fromMessage(message.matrices(i));
    }
}

#endif

}  // namespace corbo
