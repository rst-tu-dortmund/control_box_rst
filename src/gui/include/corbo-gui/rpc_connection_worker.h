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

#ifndef SRC_GUI_INCLUDE_CORBO_GUI_RPC_CONNECTION_WORKER_H_
#define SRC_GUI_INCLUDE_CORBO_GUI_RPC_CONNECTION_WORKER_H_

#include <corbo-communication/main_service_client.h>

#include <QObject>
#include <QThread>
#include <QTimer>

#include <memory>

namespace corbo {
namespace gui {

class RpcConnectionWorker : public QObject
{
    Q_OBJECT

 public:
    explicit RpcConnectionWorker(std::shared_ptr<MasterServiceClient> rpc_client, QObject* parent = nullptr) : QObject(parent)
    {
        if (rpc_client)
            _rpc_client = rpc_client;
        else
            _rpc_client = std::make_shared<MasterServiceClient>();
    }

    ~RpcConnectionWorker()
    {
        if (_timer) _timer->stop();
    }

 public slots:

    void connectRpc(const QString& address)
    {
        _server_address = address;

        if (!_timer)
        {
            _timer = new QTimer(this);
            connect(_timer, SIGNAL(timeout()), this, SLOT(tryConnect()));
            _timer->setInterval(1000);
        }

        _timer->start();  // this should also stop any previous timer
    }

 signals:
    void connectionResult(std::shared_ptr<MasterServiceClient> client);

 private slots:
    void tryConnect()
    {
        std::shared_ptr<grpc::Channel> channel = grpc::CreateChannel(_server_address.toStdString(), grpc::InsecureChannelCredentials());
        _rpc_client->setChannel(channel);

        if (_rpc_client->ping(500))
        {
            _timer->stop();
            emit connectionResult(_rpc_client);
        }
    }

 private:
    QTimer* _timer = nullptr;

    QString _server_address;

    std::shared_ptr<MasterServiceClient> _rpc_client;
};

}  // namespace gui
}  // namespace corbo

#endif  // SRC_GUI_INCLUDE_CORBO_GUI_RPC_CONNECTION_WORKER_H_
