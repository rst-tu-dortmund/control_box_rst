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

#ifndef SRC_COMMUNICATION_INCLUDE_CORBO_COMMUNICATION_MAIN_SERVICE_CLIENT_H_
#define SRC_COMMUNICATION_INCLUDE_CORBO_COMMUNICATION_MAIN_SERVICE_CLIENT_H_

#ifdef RPC_SUPPORT

#include <corbo-communication/services/master_service.grpc.pb.h>
#include <corbo-core/console.h>
#include <grpc++/grpc++.h>
#include <functional>
#include <memory>

namespace corbo {

/**
 * @brief General service client for rpc communication
 *
 * @ingroup communication rpc master
 *
 * This class implements a client for (g)RPC communication
 * with the corbo-master (defined in the master module).
 * The client might be utiltized for a GUI implementation
 * or command line based network communication.
 *
 * @see Master
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class MasterServiceClient
{
 public:
    //! Default constructor
    MasterServiceClient() = default;

    //! Default destructor
    virtual ~MasterServiceClient() = default;

    /**
     * @brief Construct an rpc service client with a given channel
     *
     * @code
     *  std::shared_ptr<grpc::Channel> channel = grpc::CreateChannel("localhost:50051", grpc::InsecureChannelCredentials());
     *  MasterServiceClient client(channel);
     * @endcode
     * @param channel   rpc channel of type grpc::Channel (refer to grpc docs)
     */
    explicit MasterServiceClient(std::shared_ptr<grpc::Channel> channel) : _stub(services::MasterService::NewStub(channel)) {}

    //! Set channel to server
    void setChannel(std::shared_ptr<grpc::Channel> channel) { _stub = services::MasterService::NewStub(channel); }

    /**
     * @brief Get the a new instance of the complete parameter message type
     *
     * This function might be overriden in subclasses to provide customized
     * parameter protobuf messages.
     * However, setParameters() and getParamters() must be updated accordingly in order to
     * do not break the general API.
     *
     * The default parameter message is corbo::messages::corboParameters.
     *
     * @see setParameters getParameters getInstance
     *
     * @return Empty message of the complete parameter message type
     */
    virtual std::unique_ptr<google::protobuf::Message> createParameterMsg() const;

    /**
     * @brief Set modules according to combined parameter message
     *
     * @see createParameterMsg getParameters
     *
     * @param[in]  parameters    Parameter message according to createParameterMsg()
     * @param[out] status_msg    Response (status message)
     * @param[in] timeout_ms     Maximum duration in ms to wait for a server response
     * @return true if service call was successful, false if not (e.g. connection lost after timeout)
     */
    virtual bool setParameters(const google::protobuf::Message& parameters, messages::Status& status_msg, int timeout_ms = 1000);

    /**
     * @brief Get the currently selected modules via combined parameter message
     * @see createParameterMsg setParameters
     * @param[out] parameters     Response (combined parameter message according to createParameterMsg())
     * @param[in] timeout_ms     Maximum duration in ms to wait for a server response
     * @return true if service call was successful, false if not (e.g. connection lost after timeout)
     */
    virtual bool getParameters(google::protobuf::Message& parameters, int timeout_ms = 1000);

    /**
     * @brief Set a plant via message at server
     * @param[in]  plant_msg     Plant message
     * @param[out] status_msg    Response (status message)
     * @param[in] timeout_ms     Maximum duration in ms to wait for a server response
     * @return true if service call was successful, false if not (e.g. connection lost after timeout)
     */
    bool setPlant(const messages::Plant& plant_msg, messages::Status& status_msg, int timeout_ms = 500);

    /**
     * @brief Get the currently selected plant via message
     * @param[out] plant_msg     Response (plant message)
     * @param[in] timeout_ms     Maximum duration in ms to wait for a server response
     * @return true if service call was successful, false if not (e.g. connection lost after timeout)
     */
    bool getPlant(messages::Plant& plant_msg, int timeout_ms = 500);

    /**
     * @brief Set a controller via message at server
     * @param[in]  ctrl_msg      Controller message
     * @param[out] status_msg    Response (status message)
     * @param[in] timeout_ms     Maximum duration in ms to wait for a server response
     * @return true if service call was successful, false if not (e.g. connection lost after timeout)
     */
    bool setController(const messages::Controller& ctrl_msg, messages::Status& status_msg, int timeout_ms = 500);

    /**
     * @brief Get the currently selected controller via message
     * @param[out] ctrl_msg      Response (controller message)
     * @param[in] timeout_ms     Maximum duration in ms to wait for a server response
     * @return true if service call was successful, false if not (e.g. connection lost after timeout)
     */
    bool getController(messages::Controller& ctrl_msg, int timeout_ms = 500);

    /**
     * @brief Set an observer via message at server
     * @param[in]  obs_msg       Observer message
     * @param[out] status_msg    Response (status message)
     * @param[in] timeout_ms     Maximum duration in ms to wait for a server response
     * @return true if service call was successful, false if not (e.g. connection lost after timeout)
     */
    bool setObserver(const messages::Observer& obs_msg, messages::Status& status_msg, int timeout_ms = 500);

    /**
     * @brief Get the currently selected observer via message
     * @param[out] obs_msg      Response (observer message)
     * @param[in] timeout_ms    Maximum duration in ms to wait for a server response
     * @return true if service call was successful, false if not (e.g. connection lost after timeout)
     */
    bool getObserver(messages::Observer& obs_msg, int timeout_ms = 500);

    /**
     * @brief Set a task via message at server
     * @param[in]  task_msg      Task message
     * @param[out] status_msg    Response (status message)
     * @param[in] timeout_ms     Maximum duration in ms to wait for a server response
     * @return true if service call was successful, false if not (e.g. connection lost after timeout)
     */
    bool setTask(const messages::Task& task_msg, messages::Status& status_msg, int timeout_ms = 500);

    /**
     * @brief Get the currently selected task via message
     * @param[out] task_msg     Response (task message)
     * @param[in] timeout_ms    Maximum duration in ms to wait for a server response
     * @return true if service call was successful, false if not (e.g. connection lost after timeout)
     */
    bool getTask(messages::Task& task_msg, int timeout_ms = 500);

    /**
     * @brief Retrieve available signals
     *
     * The service call accepts a callback with prototype
     * \code void(const messages::Signal&) \endcode
     * The callback is invoked for every registered signal.
     * Note, not all signals are registered in advance,
     * some might occur first at task execution stage.
     *
     * @param[in] feedback     Callback for available signals
     * @param[in] timeout_ms   Maximum duration in ms to wait for a server response
     * @return true if service call was successful, false if not (e.g. connection lost after timeout)
     */
    bool getAvailableSignals(std::function<void(const messages::Signal&)> feedback, int timeout_ms = 500);

    /**
     * @brief Perform the currently selected task
     *
     * Note, the master needs to be configured approriately.
     * You might need to invoke verifyConfig() first and
     * investigate the resulting status message.
     *
     * The service call accepts a callback with prototype
     * \code void(const messages::Signal&) \endcode
     * The callback is invoked for every occuring signal
     * (buffered stream stream during exeuction).
     *
     * Note, this methods blocks until the task exeuction is completed.
     * In case the tasks execution loop validates the global variable
     * ok() frequently (declared in corbo-core/global.h),
     * you can possibly interrupt execution (e.g. from another thread)
     * by invoking the service call stopTask().
     *
     * @param[in] feedback     Callback for available signals
     * @return true if service call was successful, false if not (e.g. connection lost after timeout)
     */
    bool performTask(std::function<void(const messages::Signal&)> feedback, std::string* msg = nullptr,
                     int timeout_ms = 99999999);  // TODO(roesmann): timeout != deadline!!!

    /**
     * @brief Verify if the current master configuration seems to be valid.
     *
     * Check if the master detects some issues and faults in the current setup
     * including choice for plant, controller, observer and task.
     * The resulting status message should contain a non-empty text field
     * if any issues occured.
     * @param[out] status_msg  Response containing possible issues
     * @return true if service call was successful, false if not (e.g. connection lost after timeout)
     */
    bool verifyConfig(messages::Status& status_msg, int timeout_ms = 2000);

    /**
     * @brief Ping master (check connection state)
     * @param[in] timeout_ms  maximum duration in ms to wait for a server response
     * @return true if service call was successful, false if not (e.g. connection lost after timeout)
     */
    bool ping(int timeout_ms);

    /**
     * @brief Stop execution of the current task if supported by the task
     *
     * In case the tasks execution loop validates the global variable
     * ok() frequently (declared in corbo-core/global.h),
     * you can possibly interrupt execution (e.g. from another thread)
     * by invoking this service call.
     * @return true if service call was successful, false if not (e.g. connection lost after timeout)
     */
    bool stopTask(int timeout_ms = 2000);

 private:
    std::unique_ptr<services::MasterService::Stub> _stub;
};

}  // namespace corbo

#endif  // RPC_SUPPORT

#endif  // SRC_COMMUNICATION_INCLUDE_CORBO_COMMUNICATION_MAIN_SERVICE_CLIENT_H_
