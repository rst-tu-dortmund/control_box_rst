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

#include <corbo-communication/main_service_client.h>

#include <corbo-communication/messages/corbo_parameters.pb.h>

#include <cmath>
#include <memory>

namespace corbo {

std::unique_ptr<google::protobuf::Message> MasterServiceClient::createParameterMsg() const
{
    return std::unique_ptr<messages::corboParameters>(new messages::corboParameters);
}

bool MasterServiceClient::setParameters(const google::protobuf::Message& parameters, messages::Status& status_msg, int timeout_ms)
{
    const messages::corboParameters* params = dynamic_cast<const messages::corboParameters*>(&parameters);
    if (!params)
    {
        status_msg.set_ok(false);
        status_msg.set_text("Parameter message is not of type corbo::messages::corboParameters.");
        return false;
    }
    bool ret_val = true;

    int timeout_individual = std::floor((double)timeout_ms / 4);

    messages::Status sub_status;
    ret_val = setPlant(params->plant(), sub_status, timeout_individual);
    status_msg.set_ok(status_msg.ok() && sub_status.ok());
    if (!sub_status.text().empty()) status_msg.set_text(status_msg.text() + "\n" + sub_status.text());

    ret_val = setController(params->controller(), sub_status, timeout_individual);
    status_msg.set_ok(status_msg.ok() && sub_status.ok());
    if (!sub_status.text().empty()) status_msg.set_text(status_msg.text() + "\n" + sub_status.text());

    ret_val = setObserver(params->observer(), status_msg, timeout_individual);
    status_msg.set_ok(status_msg.ok() && sub_status.ok());
    if (!sub_status.text().empty()) status_msg.set_text(status_msg.text() + "\n" + sub_status.text());

    ret_val = setTask(params->task(), status_msg, timeout_individual);
    status_msg.set_ok(status_msg.ok() && sub_status.ok());
    if (!sub_status.text().empty()) status_msg.set_text(status_msg.text() + "\n" + sub_status.text());

    return ret_val;
}

bool MasterServiceClient::getParameters(google::protobuf::Message& parameters, int timeout_ms)
{
    parameters.Clear();

    messages::corboParameters* params = dynamic_cast<messages::corboParameters*>(&parameters);
    if (!params)
    {
        PRINT_ERROR("MasterServiceClient::getParameters(): Message type is not compatible!");
        return false;
    }

    int timeout_individual = std::floor((double)timeout_ms / 4);

    bool ret_val = true;

    ret_val = getPlant(*params->mutable_plant(), timeout_individual);

    ret_val = getController(*params->mutable_controller(), timeout_individual);

    ret_val = getObserver(*params->mutable_observer(), timeout_individual);

    ret_val = getTask(*params->mutable_task(), timeout_individual);

    return ret_val;
}

bool MasterServiceClient::setPlant(const messages::Plant& plant_msg, messages::Status& status_msg, int timeout_ms)
{
    if (!_stub) return false;
    grpc::ClientContext context;
    std::chrono::system_clock::time_point deadline = std::chrono::system_clock::now() + std::chrono::milliseconds(timeout_ms);
    context.set_deadline(deadline);
    grpc::Status status = _stub->setPlant(&context, plant_msg, &status_msg);
    if (!status.ok())
    {
        return false;
    }
    return true;
}

bool MasterServiceClient::getPlant(messages::Plant& plant_msg, int timeout_ms)
{
    if (!_stub) return false;
    grpc::ClientContext context;
    std::chrono::system_clock::time_point deadline = std::chrono::system_clock::now() + std::chrono::milliseconds(timeout_ms);
    context.set_deadline(deadline);
    messages::Void void_msg;
    grpc::Status status = _stub->getPlant(&context, void_msg, &plant_msg);
    if (!status.ok())
    {
        return false;
    }
    return true;
}

bool MasterServiceClient::setController(const messages::Controller& ctrl_msg, messages::Status& status_msg, int timeout_ms)
{
    if (!_stub) return false;
    grpc::ClientContext context;
    std::chrono::system_clock::time_point deadline = std::chrono::system_clock::now() + std::chrono::milliseconds(timeout_ms);
    context.set_deadline(deadline);
    grpc::Status status = _stub->setController(&context, ctrl_msg, &status_msg);
    if (!status.ok())
    {
        return false;
    }
    return true;
}

bool MasterServiceClient::getController(messages::Controller& ctrl_msg, int timeout_ms)
{
    if (!_stub) return false;
    grpc::ClientContext context;
    std::chrono::system_clock::time_point deadline = std::chrono::system_clock::now() + std::chrono::milliseconds(timeout_ms);
    context.set_deadline(deadline);
    messages::Void void_msg;
    grpc::Status status = _stub->getController(&context, void_msg, &ctrl_msg);
    if (!status.ok())
    {
        return false;
    }
    return true;
}

bool MasterServiceClient::setObserver(const messages::Observer& obs_msg, messages::Status& status_msg, int timeout_ms)
{
    if (!_stub) return false;
    grpc::ClientContext context;
    std::chrono::system_clock::time_point deadline = std::chrono::system_clock::now() + std::chrono::milliseconds(timeout_ms);
    context.set_deadline(deadline);
    grpc::Status status = _stub->setObserver(&context, obs_msg, &status_msg);
    if (!status.ok())
    {
        return false;
    }
    return true;
}

bool MasterServiceClient::getObserver(messages::Observer& obs_msg, int timeout_ms)
{
    if (!_stub) return false;
    grpc::ClientContext context;
    std::chrono::system_clock::time_point deadline = std::chrono::system_clock::now() + std::chrono::milliseconds(timeout_ms);
    context.set_deadline(deadline);
    messages::Void void_msg;
    grpc::Status status = _stub->getObserver(&context, void_msg, &obs_msg);
    if (!status.ok())
    {
        return false;
    }
    return true;
}

bool MasterServiceClient::setTask(const messages::Task& task_msg, messages::Status& status_msg, int timeout_ms)
{
    if (!_stub) return false;
    grpc::ClientContext context;
    std::chrono::system_clock::time_point deadline = std::chrono::system_clock::now() + std::chrono::milliseconds(timeout_ms);
    context.set_deadline(deadline);
    grpc::Status status = _stub->setTask(&context, task_msg, &status_msg);
    if (!status.ok())
    {
        return false;
    }
    return true;
}

bool MasterServiceClient::getTask(messages::Task& task_msg, int timeout_ms)
{
    if (!_stub) return false;
    grpc::ClientContext context;
    std::chrono::system_clock::time_point deadline = std::chrono::system_clock::now() + std::chrono::milliseconds(timeout_ms);
    context.set_deadline(deadline);
    messages::Void void_msg;
    grpc::Status status = _stub->getTask(&context, void_msg, &task_msg);
    if (!status.ok())
    {
        return false;
    }
    return true;
}

bool MasterServiceClient::getAvailableSignals(std::function<void(const messages::Signal& signal)> feedback, int timeout_ms)
{
    if (!_stub) return false;
    grpc::ClientContext context;
    std::chrono::system_clock::time_point deadline = std::chrono::system_clock::now() + std::chrono::milliseconds(timeout_ms);
    context.set_deadline(deadline);
    messages::Void void_msg;
    messages::Signal signal;
    std::unique_ptr<grpc::ClientReader<messages::Signal>> reader(_stub->getAvailableSignals(&context, void_msg));
    while (reader->Read(&signal))
    {
        feedback(signal);
    }
    grpc::Status status = reader->Finish();
    if (!status.ok()) return false;
    return true;
}

bool MasterServiceClient::performTask(std::function<void(const messages::Signal& signal)> feedback, std::string* msg, int timeout_ms)
{
    if (!_stub) return false;
    grpc::ClientContext context;
    std::chrono::system_clock::time_point deadline = std::chrono::system_clock::now() + std::chrono::milliseconds(timeout_ms);
    context.set_deadline(deadline);

    messages::Void void_msg;
    messages::Signal signal;

    std::unique_ptr<grpc::ClientReader<messages::Signal>> reader(_stub->performTask(&context, void_msg));
    while (reader->Read(&signal))
    {
        feedback(signal);
    }
    grpc::Status status = reader->Finish();

    if (!status.ok())
    {
        if (msg) *msg = status.error_message();
        return false;
    }
    return true;
}

bool MasterServiceClient::verifyConfig(messages::Status& status_msg, int timeout_ms)
{
    if (!_stub) return false;
    grpc::ClientContext context;
    std::chrono::system_clock::time_point deadline = std::chrono::system_clock::now() + std::chrono::milliseconds(timeout_ms);
    context.set_deadline(deadline);
    messages::Void void_msg;
    grpc::Status status = _stub->verifyConfig(&context, void_msg, &status_msg);
    if (!status.ok())
    {
        return false;
    }
    return true;
}

bool MasterServiceClient::ping(int timeout_ms)
{
    if (!_stub) return false;
    grpc::ClientContext context;
    // context.set_fail_fast(true);
    std::chrono::system_clock::time_point deadline = std::chrono::system_clock::now() + std::chrono::milliseconds(timeout_ms);
    context.set_deadline(deadline);
    messages::Void void_msg;
    messages::Status reponse;
    grpc::Status status = _stub->ping(&context, void_msg, &reponse);
    if (!status.ok())
    {
        return false;
    }
    return true;
}

bool MasterServiceClient::stopTask(int timeout_ms)
{
    if (!_stub) return false;
    grpc::ClientContext context;
    // context.set_fail_fast(true);
    std::chrono::system_clock::time_point deadline = std::chrono::system_clock::now() + std::chrono::milliseconds(timeout_ms);
    context.set_deadline(deadline);
    messages::Void void_msg;
    messages::Void reponse;
    grpc::Status status = _stub->stop(&context, void_msg, &reponse);
    if (!status.ok())
    {
        return false;
    }
    return true;
}

}  // namespace corbo

// RPC_SUPPORT
#endif
