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

#ifndef SRC_COMMUNICATION_INCLUDE_CORBO_COMMUNICATION_SIGNAL_TARGET_RPC_H_
#define SRC_COMMUNICATION_INCLUDE_CORBO_COMMUNICATION_SIGNAL_TARGET_RPC_H_

#ifdef RPC_SUPPORT

#include <corbo-communication/services/master_service.grpc.pb.h>
#include <corbo-core/signal_target_interface.h>

#include <memory>
#include <string>
#include <vector>

namespace corbo {

/**
 * @brief Signal target for RPC communication
 *
 * @ingroup communication signals
 *
 * This signal target implements a SignalTargetInterface
 * for RPC communication, in particular gRPC is chosed
 * as actual implementation.
 * Signals are directly passed to grpc::ServerWriter which
 * realizes an output stream to a desired target.
 * gRPC already takes care of buffering and multiple threading.
 *
 * @see SignalTargetInterface CommonSignalTarget
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class SignalTargetRPC : public SignalTargetInterface
{
 public:
    using Ptr = std::shared_ptr<SignalTargetRPC>;

    /**
     * @brief Construct signal target with a given grpc::ServerWriter instance
     * @param[in] signal_writer (refer to the gRPC doc)
     */
    explicit SignalTargetRPC(grpc::ServerWriter<messages::Signal>* signal_writer) : _signal_writer(signal_writer) {}

    // implements interface method
    void registerMeasurement(const std::string& unique_name, int value_dimension, const std::vector<std::string>& value_labels = {},
                             bool zero_order_hold = false) override;
    // implements interface method
    void registerTimeSeries(const std::string& unique_name, int value_dimension, bool zero_order_hold = false) override;

    // implements interface method
    void sendMeasurement(Measurement::ConstPtr measurement) override;
    // implements interface method
    void sendTimeSeries(TimeSeriesSignal::Ptr time_series) override;

    // implements interface method
    void sendIndexedValues(IndexedValuesSignal::Ptr indexed_values) override;

    // implements interface method
    void sendIndexedValuesSet(IndexedValuesSetSignal::Ptr indexed_values_set) override;

    // implements interface method
    void sendMatrix(MatrixSignal::Ptr matrix) override;

 private:
    grpc::ServerWriter<messages::Signal>* _signal_writer;
};

}  // namespace corbo

#endif  // RPC_SUPPORT

#endif  // SRC_COMMUNICATION_INCLUDE_CORBO_COMMUNICATION_SIGNAL_TARGET_RPC_H_
