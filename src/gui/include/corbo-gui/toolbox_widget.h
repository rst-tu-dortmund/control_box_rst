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

#ifndef SRC_GUI_INCLUDE_CORBO_GUI_TOOLBOX_WIDGET_H_
#define SRC_GUI_INCLUDE_CORBO_GUI_TOOLBOX_WIDGET_H_

#include <corbo-communication/main_service_client.h>
#include <corbo-communication/messages/corbo_parameters.pb.h>
#include <corbo-gui/parameter_cache.h>
#include <corbo-gui/parameter_widget.h>
#include <corbo-gui/signal_helper.h>

#include <QDateTime>
#include <QPushButton>
#include <QString>
#include <QThread>
#include <QToolBox>
#include <QVector>

#include <memory>

namespace corbo {
namespace gui {

class ToolboxWidget : public QWidget
{
    Q_OBJECT

 public:
    explicit ToolboxWidget(SignalHelper::Ptr signal_helper, std::shared_ptr<MasterServiceClient> rpc, QWidget* parent = 0);
    virtual ~ToolboxWidget();

    QSize sizeHint() const override;

    QVector<ParameterWidget*> getParameterWidgets() { return _param_widgets; }

    void setSignalAutoRequest(bool active) { _signal_auto_update = active; }

    void setRpcClient(std::shared_ptr<MasterServiceClient> rpc_client);

    void fromMessage(const google::protobuf::Message& parameters);
    void toMessage(google::protobuf::Message& parameters) const;

    void clearParameters();

 signals:
    void taskCompleted(const QDateTime& time);
    void parameterLoadedFromFile(const QString& filename);
    void parameterSavedToFile(const QString& filename);

 public slots:
    void getAvailableSignals();
    void requestParametersFromService();
    void restoreTaskParameters(int task_id);
    void removeFromTaskCache(int task_id);
    void saveCurrentParametersToFile();
    void saveParametersToFile(int task_id);
    void loadParametersFromFile();
    void loadParametersFromFile(const QString& filepath);

 protected:
    void createToolboxContent();
    void initializeToolboxes();
    void updateParameterWidgets();
    bool sendParametersToService(bool suppress_warning = false);
    void saveParametersToFile(const google::protobuf::Message& params);

    void setTaskRunningFlag(bool running);

 protected slots:
    void performTask();
    void rpcTaskFinished(bool success, const QString& msg);

 private:
    QToolBox* _toolbox;
    QPushButton* _btn_perform_task;

    std::unique_ptr<google::protobuf::Message> _param_message;  //!< Internal cache of the current parameter setting/message
    QVector<ParameterWidget*> _param_widgets;

    ParameterCache _param_cache;
    ParameterCache _task_cache;  // caches paramters that belongs to a single task run

    SignalHelper::Ptr _signal_helper;
    std::shared_ptr<MasterServiceClient> _rpc_client;

    bool _signal_auto_update = false;  // this should be false for the initial messages parsing, but enabled afterwards

    QThread _rpc_task_thread;
    bool _task_running = false;  //!< Keep track if a task is currently running

    const QString _task_button_label_active   = "Perform Task";
    const QString _task_button_label_inactive = "Stop Task";
};

}  // namespace gui
}  // namespace corbo

#endif  // SRC_GUI_INCLUDE_CORBO_GUI_TOOLBOX_WIDGET_H_
