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

#ifndef SRC_MASTER_INCLUDE_CORBO_MASTER_MASTER_H_
#define SRC_MASTER_INCLUDE_CORBO_MASTER_MASTER_H_

#include <corbo-communication/messages/corbo_parameters.pb.h>
#include <corbo-communication/services/master_service.grpc.pb.h>
#include <corbo-controllers/corbo_controllers.h>
#include <corbo-core/console.h>
#include <corbo-numerics/corbo_numerics.h>
#include <corbo-observers/corbo_observers.h>
#include <corbo-optimal-control/corbo_optimal_control.h>
#include <corbo-optimization/corbo_optimization.h>
#include <corbo-plants/corbo_plants.h>
#include <corbo-systems/corbo_systems.h>
#include <corbo-tasks/corbo_tasks.h>

#include <grpc++/grpc++.h>

#include <iostream>
#include <memory>
#include <string>

namespace corbo {

/**
 * @brief General service client for rpc communication
 *
 * @ingroup master communication rpc
 *
 * This class implements a master for (g)RPC communication
 * with the main service client (MasterServiceClient, defined in the communication module).
 *
 * @remark Make sure that any executable that instanciates a Master also includes
 *         relevant module headers in order to register classes at their corresponding
 *         factories.
 *
 * @see MasterServiceClient
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class Master : public services::MasterService::Service
{
 public:
    //! Default constructor
    Master() {}

    /**
     * @brief Start master server (blocking call)
     *
     * Start master with a given server address,
     * e.g. "localhost:50051".
     * This call blocks further program execution.
     *
     * @param[in] server_address  Server address including port.
     * @param[in] blocking        If true, the method blocks thread until shutdown
     */
    void start(const std::string& server_address, bool blocking = true);

    //! Restore default master settings
    void setDefault();

    bool loadFromFile(const std::string& filename);

    bool setcorboMessage(const corbo::messages::corboParameters& msg, std::stringstream* issues);

    bool setPlant(const messages::Plant& msg, std::stringstream* issues = nullptr);
    bool setController(const corbo::messages::Controller& msg, std::stringstream* issues = nullptr);
    bool setObserver(const corbo::messages::Observer& msg, std::stringstream* issues = nullptr);
    bool setTask(const corbo::messages::Task& msg, std::stringstream* issues = nullptr);

 protected:
    //! Set and configure plant
    grpc::Status setPlant(grpc::ServerContext* context, const corbo::messages::Plant* request, corbo::messages::Status* response) override;
    //! Get current plant configuration
    grpc::Status getPlant(grpc::ServerContext* context, const corbo::messages::Void* request, corbo::messages::Plant* response) override;

    //! Set and configure controller
    grpc::Status setController(grpc::ServerContext* context, const corbo::messages::Controller* request, corbo::messages::Status* response) override;
    //! Get current controller configuration
    grpc::Status getController(grpc::ServerContext* context, const corbo::messages::Void* request, corbo::messages::Controller* response) override;

    //! Set and configure observer
    grpc::Status setObserver(grpc::ServerContext* context, const corbo::messages::Observer* request, corbo::messages::Status* response) override;
    //! Get current observer configuration
    grpc::Status getObserver(grpc::ServerContext* context, const corbo::messages::Void* request, corbo::messages::Observer* response) override;

    //! Set and configure task
    grpc::Status setTask(grpc::ServerContext* context, const corbo::messages::Task* request, corbo::messages::Status* response) override;
    //! Get current task configuration
    grpc::Status getTask(grpc::ServerContext* context, const corbo::messages::Void* request, corbo::messages::Task* response) override;

    //! Retrieve available signals from task (and included modules)
    grpc::Status getAvailableSignals(grpc::ServerContext* context, const corbo::messages::Void* request,
                                     grpc::ServerWriter<corbo::messages::Signal>* response_stream) override;

    //! Perform the current task and broadcast stream of signals.
    grpc::Status performTask(grpc::ServerContext* context, const corbo::messages::Void* request,
                             grpc::ServerWriter<corbo::messages::Signal>* response_stream) override;

    //! Check if the current configuration seems to be valid (refer to TaskInterface::verify())
    grpc::Status verifyConfig(grpc::ServerContext* context, const corbo::messages::Void* request, corbo::messages::Status* response) override;

    //! Send and receive dummy message in order to check the connectivity between master and client
    grpc::Status ping(grpc::ServerContext* context, const corbo::messages::Void* request, corbo::messages::Status* response) override;
    //! Execution stop of the current task requested (setting global ok() variable to false, needs to be supported by the selected task)
    grpc::Status stop(grpc::ServerContext* context, const corbo::messages::Void* request, corbo::messages::Void* response) override;

 private:
    Environment _environment;
    TaskInterface::Ptr _task;

    std::unique_ptr<grpc::Server> _server;
};

}  // namespace corbo

#endif  // SRC_MASTER_INCLUDE_CORBO_MASTER_MASTER_H_
