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

#ifndef SRC_GUI_INCLUDE_CORBO_GUI_RPC_TASK_WORKER_H_
#define SRC_GUI_INCLUDE_CORBO_GUI_RPC_TASK_WORKER_H_

#include <corbo-communication/main_service_client.h>

#include <QObject>
#include <QThread>
#include <QTimer>

#include <memory>

namespace corbo {
namespace gui {

class RpcTaskWorker : public QObject
{
    Q_OBJECT

 public:
    explicit RpcTaskWorker(QObject* parent = nullptr) : QObject(parent) {}

    ~RpcTaskWorker() {}

    void setRpcClient(std::shared_ptr<MasterServiceClient> client) { _rpc_client = client; }

 public slots:

    void performTask()
    {
        if (!_rpc_client) emit taskFinished(false, "RPC Client not connected");
        std::string msg;

        auto signal_feedback = [this](const messages::Signal& signal) { emit feedback(signal); };
        bool success         = _rpc_client->performTask(signal_feedback, &msg);

        emit taskFinished(success, QString::fromStdString(msg));
    }

 signals:
    void feedback(const messages::Signal& signal);
    void taskFinished(bool success, const QString& msg);

 private:
    std::shared_ptr<MasterServiceClient> _rpc_client;
};

}  // namespace gui
}  // namespace corbo

#endif  // SRC_GUI_INCLUDE_CORBO_GUI_RPC_TASK_WORKER_H_
