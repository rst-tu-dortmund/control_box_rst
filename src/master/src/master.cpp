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

#include <corbo-master/master.h>

#include <corbo-communication/signal_target_rpc.h>
#include <corbo-communication/utilities.h>
#include <corbo-core/global.h>

#include <fstream>
#include <memory>
#include <string>

namespace corbo {

void Master::start(const std::string& server_address, bool blocking)
{
    if (!_environment.hasController() && !_environment.hasObserver() && !_environment.hasPlant()) setDefault();

    grpc::ServerBuilder builder;

    // Listen on the given address without any authentication mechanism.
    builder.AddListeningPort(server_address, grpc::InsecureServerCredentials());

    // Create main service
    builder.RegisterService(this);

    _server = builder.BuildAndStart();
    if (_server)
    {
        PRINT_INFO("RpcServer listening on " << server_address);
        if (blocking) _server->Wait();
    }
    else
    {
        PRINT_ERROR("RpcServer failed");
    }
}

void Master::setDefault()
{
    PRINT_DEBUG_NAMED("Initializing default environment.");

    // controller
    ControllerFactory::instance().printObjects();
    ControllerInterface::Ptr controller = ControllerFactory::instance().create("PidController");
    if (controller) _environment.setController(controller);

    // observer
    ObserverFactory::instance().printObjects();
    ObserverInterface::Ptr observer = ObserverFactory::instance().create("NoObserver");
    if (observer) _environment.setObserver(observer);

    corbo::SerialIntegratorSystem::Ptr system = std::make_shared<corbo::SerialIntegratorSystem>();
    corbo::FullStateSystemOutput::Ptr output  = std::make_shared<corbo::FirstStateSystemOutput>();
    corbo::PlantInterface::Ptr plant          = std::make_shared<corbo::SimulatedPlant>(system, output);
    if (plant) _environment.setPlant(plant);

    // task
    TaskFactory::instance().printObjects();
    _task = TaskFactory::instance().create("ClosedLoopControlTask");
}

bool Master::loadFromFile(const std::string& filename)
{
    if (filename.empty()) return false;

    corbo::messages::corboParameters msg;

    std::ifstream file(filename);

    if (!file.is_open())
    {
        return false;
    }

    if (!msg.ParsePartialFromIstream(&file))
    {
        return false;
    }

    std::stringstream issues;
    if (!setcorboMessage(msg, &issues))
    {
        PRINT_ERROR_NAMED(issues.str());
        return false;
    }

    return true;
}

bool Master::setcorboMessage(const messages::corboParameters& msg, std::stringstream* issues)
{
    bool result = setPlant(msg.plant(), issues);
    PRINT_ERROR_COND_NAMED(!result, "setPlant failed.");

    result &= setController(msg.controller(), issues);
    PRINT_ERROR_COND_NAMED(!result, "setController failed.");

    result &= setObserver(msg.observer(), issues);
    PRINT_ERROR_COND_NAMED(!result, "setObserver failed.");

    result &= setTask(msg.task(), issues);
    PRINT_ERROR_COND_NAMED(!result, "setTask failed.");

    return result;
}

bool Master::setPlant(const messages::Plant& msg, std::stringstream* issues)
{
    // if (!_environment.getPlant()) return false;  // grpc::Status::CANCELLED;

    if (msg.plant_case() == corbo::messages::Plant::PLANT_NOT_SET)
    {
        if (issues) *issues << "No plant selected." << std::endl;
        return false;
    }
    std::string plant_type;
    if (util::get_oneof_field_type(msg, "plant", plant_type, false))
    {
        PlantInterface::Ptr plant = PlantFactory::instance().create(plant_type, false);
        if (!plant)
        {
            // also check if we have nested, isolated one-of messages
            util::get_oneof_field_type_expand_isolated(msg, "plant", plant_type, false);
            plant = PlantFactory::instance().create(plant_type);
        }

        if (plant)
        {
            plant->fromMessage(msg, issues);
            _environment.setPlant(plant);
            return issues ? issues->str().empty() : true;  // TODO(roesmann) correct return value
        }
        if (issues) *issues << "Plant type '" + plant_type + "' not found." << std::endl;
    }
    else if (issues)
        *issues << "Cannot determine type of specified plant. Protobuf error." << std::endl;

    return false;
}

bool Master::setController(const messages::Controller& msg, std::stringstream* issues)
{
    // if (!_environment.getController()) return false;  // grpc::Status::CANCELLED;

    if (msg.controller_case() == corbo::messages::Controller::CONTROLLER_NOT_SET)
    {
        if (issues) *issues << "No controller selected." << std::endl;
        return false;
    }
    std::string controller_type;
    if (util::get_oneof_field_type(msg, "controller", controller_type, false))
    {
        ControllerInterface::Ptr controller = ControllerFactory::instance().create(controller_type, false);
        if (!controller)
        {
            // also check if we have nested, isolated one-of messages
            util::get_oneof_field_type_expand_isolated(msg, "controller", controller_type, false);
            controller = ControllerFactory::instance().create(controller_type);
        }

        if (controller)
        {
            controller->fromMessage(msg, issues);
            _environment.setController(controller);
            return issues ? issues->str().empty() : true;  // TODO(roesmann) correct return value
        }
        if (issues) *issues << "Controller type '" + controller_type + "' not found." << std::endl;
    }
    else if (issues)
        *issues << "Cannot determine type of specified controller. Protobuf error." << std::endl;

    return false;
}

bool Master::setObserver(const messages::Observer& msg, std::stringstream* issues)
{
    // if (!_environment.getObserver()) return false;  // grpc::Status::CANCELLED;

    if (msg.observer_case() == corbo::messages::Observer::OBSERVER_NOT_SET)
    {
        if (issues) *issues << "No observer selected." << std::endl;
        return false;
    }
    std::string observer_type;
    if (util::get_oneof_field_type(msg, "observer", observer_type, false))
    {
        ObserverInterface::Ptr observer = ObserverFactory::instance().create(observer_type, false);
        if (!observer)
        {
            // also check if we have nested, isolated one-of messages
            util::get_oneof_field_type_expand_isolated(msg, "observer", observer_type, false);
            observer = ObserverFactory::instance().create(observer_type);
        }

        if (observer)
        {
            observer->fromMessage(msg, issues);
            _environment.setObserver(observer);
            return issues ? issues->str().empty() : true;  // TODO(roesmann) correct return value
        }
        if (issues) *issues << "Observer type '" + observer_type + "' not found." << std::endl;
    }
    else if (issues)
        *issues << "Cannot determine type of specified observer. Protobuf error." << std::endl;

    return false;
}

bool Master::setTask(const messages::Task& msg, std::stringstream* issues)
{
    if (msg.task_case() == corbo::messages::Task::TASK_NOT_SET)
    {
        if (issues) *issues << "No task selected." << std::endl;
        return false;
    }
    std::string task_type;
    if (util::get_oneof_field_type(msg, "task", task_type, false))
    {
        _task = TaskFactory::instance().create(task_type, false);
        if (!_task)
        {
            // also check if we have nested, isolated one-of messages
            util::get_oneof_field_type_expand_isolated(msg, "task", task_type, false);
            _task = TaskFactory::instance().create(task_type);
        }

        if (_task)
        {
            _task->fromMessage(msg, issues);
            return issues ? issues->str().empty() : true;  // TODO(roesmann) correct return value
        }
        if (issues) *issues << "Task type '" + task_type + "' not found." << std::endl;
    }
    else if (issues)
        *issues << "Cannot determine type of specified task. Protobuf error." << std::endl;

    return false;
}

grpc::Status Master::setPlant(grpc::ServerContext* context, const corbo::messages::Plant* request, corbo::messages::Status* response)
{
    std::stringstream issues;
    if (setPlant(*request, &issues))
    {
        response->set_ok(true);
    }
    else
    {
        response->set_ok(false);
        response->set_text(issues.str());
    }
    return grpc::Status::OK;
}

grpc::Status Master::getPlant(grpc::ServerContext* context, const corbo::messages::Void* request, corbo::messages::Plant* response)
{
    if (!_environment.getPlant()) return grpc::Status::CANCELLED;
    _environment.getPlant()->toMessage(*response);
    return grpc::Status::OK;
}

grpc::Status Master::setController(grpc::ServerContext* context, const corbo::messages::Controller* request, corbo::messages::Status* response)
{
    std::stringstream issues;
    if (setController(*request, &issues))
    {
        response->set_ok(true);
    }
    else
    {
        response->set_ok(false);
        response->set_text(issues.str());
    }
    return grpc::Status::OK;
}

grpc::Status Master::getController(grpc::ServerContext* context, const corbo::messages::Void* request, corbo::messages::Controller* response)
{
    if (!_environment.getController()) return grpc::Status::CANCELLED;
    _environment.getController()->toMessage(*response);
    return grpc::Status::OK;
}

grpc::Status Master::setObserver(grpc::ServerContext* context, const corbo::messages::Observer* request, corbo::messages::Status* response)
{
    std::stringstream issues;
    if (setObserver(*request, &issues))
    {
        response->set_ok(true);
    }
    else
    {
        response->set_ok(false);
        response->set_text(issues.str());
    }
    return grpc::Status::OK;
}

grpc::Status Master::getObserver(grpc::ServerContext* context, const corbo::messages::Void* request, corbo::messages::Observer* response)
{
    if (!_environment.getObserver()) return grpc::Status::CANCELLED;
    _environment.getObserver()->toMessage(*response);
    return grpc::Status::OK;
}

grpc::Status Master::setTask(grpc::ServerContext* context, const corbo::messages::Task* request, corbo::messages::Status* response)
{
    std::stringstream issues;
    if (setTask(*request, &issues))
    {
        response->set_ok(true);
    }
    else
    {
        response->set_ok(false);
        response->set_text(issues.str());
    }
    return grpc::Status::OK;
}

grpc::Status Master::getTask(grpc::ServerContext* context, const corbo::messages::Void* request, corbo::messages::Task* response)
{
    if (!_task) return grpc::Status::CANCELLED;
    _task->toMessage(*response);
    return grpc::Status::OK;
}

grpc::Status Master::getAvailableSignals(grpc::ServerContext* context, const corbo::messages::Void* request,
                                         grpc::ServerWriter<messages::Signal>* response_stream)
{
    if (_task)
    {
        SignalTargetRPC::Ptr signal_target = std::make_shared<SignalTargetRPC>(response_stream);
        _task->getAvailableSignals(_environment, *signal_target);
    }
    return grpc::Status::OK;
}

grpc::Status Master::performTask(grpc::ServerContext* context, const corbo::messages::Void* request,
                                 grpc::ServerWriter<messages::Signal>* response_stream)
{
    setOk(true);  // reset global ok-state
    std::string err_msg;
    if (_task)
    {
        SignalTargetRPC::Ptr signal_target = std::make_shared<SignalTargetRPC>(response_stream);
        _task->performTask(_environment, signal_target.get(), &err_msg);
    }

    if (!err_msg.empty())
    {
        return grpc::Status(grpc::StatusCode::ABORTED, err_msg);
    }

    return grpc::Status::OK;
}

grpc::Status Master::verifyConfig(grpc::ServerContext* context, const corbo::messages::Void* request, corbo::messages::Status* response)
{
    if (_task)
    {
        std::string result;
        response->set_ok(_task->verify(_environment, &result));
        response->set_text(result);
    }
    else
    {
        response->set_text("No task selected!");
    }
    return grpc::Status::OK;
}

grpc::Status Master::ping(grpc::ServerContext* context, const corbo::messages::Void* request, corbo::messages::Status* response)
{
    response->set_text("got it");
    response->set_ok(true);
    return grpc::Status::OK;
}

grpc::Status Master::stop(grpc::ServerContext* context, const corbo::messages::Void* request, corbo::messages::Void* response)
{
    setOk(false);  // set global ok flag which should be checked frequently by the current task
    return grpc::Status::OK;
}

}  // namespace corbo
