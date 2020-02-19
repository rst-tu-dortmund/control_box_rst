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

#ifdef RPC_SUPPORT

#include <corbo-communication/signal_target_rpc.h>

#include <corbo-core/console.h>

#include <memory>

namespace corbo {

void SignalTargetRPC::registerMeasurement(const std::string& unique_name, int value_dimension, const std::vector<std::string>& value_labels,
                                          bool zero_order_hold)
{
    // default construct measurement with appropriate dimension
    Measurement::Ptr measurement        = std::make_shared<Measurement>();
    measurement->header.name            = unique_name;
    measurement->header.time            = Time::now();
    measurement->header.value_dimension = value_dimension;
    measurement->header.zero_order_hold = zero_order_hold;

    measurement->getValueLabelsRef() = value_labels;

    sendMeasurement(measurement);
}

void SignalTargetRPC::registerTimeSeries(const std::string& unique_name, int value_dimension, bool zero_order_hold)
{
    // default construct time_series with appropriate dimension
    TimeSeriesSignal::Ptr time_series_signal   = std::make_shared<TimeSeriesSignal>(value_dimension);
    time_series_signal->header.name            = unique_name;
    time_series_signal->header.time            = Time::now();
    time_series_signal->header.value_dimension = value_dimension;
    time_series_signal->header.zero_order_hold = zero_order_hold;
    sendTimeSeries(time_series_signal);
}

void SignalTargetRPC::sendMeasurement(Measurement::ConstPtr measurement)
{
    if (!measurement) return;

    messages::Signal signal_msg;
    measurement->toMessage(signal_msg);
    _signal_writer->Write(signal_msg);
}

void SignalTargetRPC::sendTimeSeries(TimeSeriesSignal::Ptr time_series)
{
    if (!time_series) return;

    messages::Signal signal_msg;
    time_series->toMessage(signal_msg);
    _signal_writer->Write(signal_msg);
}

void SignalTargetRPC::sendIndexedValues(IndexedValuesSignal::Ptr indexed_values)
{
    if (!indexed_values) return;
    messages::Signal signal_msg;
    indexed_values->toMessage(signal_msg);
    _signal_writer->Write(signal_msg);
}

void SignalTargetRPC::sendIndexedValuesSet(IndexedValuesSetSignal::Ptr indexed_values_set)
{
    if (!indexed_values_set) return;

    messages::Signal signal_msg;
    indexed_values_set->toMessage(signal_msg);
    _signal_writer->Write(signal_msg);
}

void SignalTargetRPC::sendMatrix(MatrixSignal::Ptr matrix)
{
    if (!matrix) return;

    messages::Signal signal_msg;
    matrix->toMessage(signal_msg);
    _signal_writer->Write(signal_msg);
}

}  // namespace corbo

// RPC_SUPPORT
#endif
