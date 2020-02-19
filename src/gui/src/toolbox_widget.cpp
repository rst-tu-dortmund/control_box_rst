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

#include <corbo-gui/toolbox_widget.h>

#include <corbo-gui/rpc_task_worker.h>

#include <QFileDialog>
#include <QMessageBox>

#include <fstream>
#include <memory>

namespace corbo {
namespace gui {

ToolboxWidget::ToolboxWidget(SignalHelper::Ptr signal_helper, std::shared_ptr<MasterServiceClient> rpc, QWidget* parent)
    : QWidget(parent), _signal_helper(signal_helper), _rpc_client(rpc)
{
    setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::Maximum);

    if (rpc) _param_message = rpc->createParameterMsg();

    createToolboxContent();

    // register meta types for queued signal-slot connections
    qRegisterMetaType<messages::Signal>("messages::Signal");
}

ToolboxWidget::~ToolboxWidget()
{
    _rpc_task_thread.quit();
    _rpc_task_thread.wait();
}

QSize ToolboxWidget::sizeHint() const { return QSize(350, 800); }

void ToolboxWidget::createToolboxContent()
{
    QVBoxLayout* layout = new QVBoxLayout(this);
    // layout->setContentsMargins(0, 0, 0, 0);
    layout->setAlignment(Qt::AlignTop);

    _toolbox = new QToolBox;
    // _toolbox->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::Preferred);

    // create toolbox items according to corboParameters message;
    initializeToolboxes();

    layout->addWidget(_toolbox);

    // add hline and button:
    layout->addSpacing(20);
    QFrame* hline = new QFrame;
    hline->setFrameShape(QFrame::HLine);
    hline->setFrameShadow(QFrame::Sunken);
    hline->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Maximum);
    layout->addWidget(hline);

    layout->addSpacing(20);

    QHBoxLayout* btn_layout = new QHBoxLayout;
    _btn_perform_task       = new QPushButton(_task_button_label_active);
    _btn_perform_task->setMaximumWidth(200);
    _btn_perform_task->setMaximumHeight(40);
    connect(_btn_perform_task, &QPushButton::clicked, [this](bool) { performTask(); });
    btn_layout->addWidget(_btn_perform_task);
    layout->addLayout(btn_layout);

    //  layout->addStretch();
}

void ToolboxWidget::initializeToolboxes()
{
    if (!_param_message) return;

    if (!_param_widgets.empty())
    {
        for (ParameterWidget* widget : _param_widgets) widget->deleteLater();
    }

    // create toolbox items according to corboParameters message;
    namespace p               = google::protobuf;
    const p::Descriptor* desc = _param_message->GetDescriptor();
    int num_fields            = desc->field_count();
    for (int i = 0; i < num_fields; ++i)
    {
        const p::FieldDescriptor* field = desc->field(i);
        if (field->cpp_type() != p::FieldDescriptor::CPPTYPE_MESSAGE)
        {
            PRINT_ERROR(
                "ToolboxWidget::createToolboxContent(): non-message-field found in parameter message. Only message fields are allowed in top-level");
            continue;
        }

        QString name = QString::fromStdString(field->name());
        name[0]      = name[0].toUpper();

        ParameterWidget* param_widget = new ParameterWidget(&_param_cache, this);
        _toolbox->addItem(param_widget, name);
        connect(param_widget, &ParameterWidget::updatedOneOfField, [this](const QString& /*text*/) {
            if (_signal_auto_update) getAvailableSignals();
        });
        connect(param_widget, &ParameterWidget::signalUpdateRequested, this, &ToolboxWidget::getAvailableSignals);

        // TODO(roesmann): We can directly modify _param_message field instead of copying each time in fromParameter and toParameters
        _param_widgets.push_back(param_widget);
    }
}

void ToolboxWidget::updateParameterWidgets()
{
    namespace p               = google::protobuf;
    const p::Descriptor* desc = _param_message->GetDescriptor();
    int num_fields            = desc->field_count();

    if (num_fields != _param_widgets.size())
    {
        PRINT_ERROR("ToolboxWidget::updateParameterWidgets(): Number of fields does not match number of ParameterWidgets.");
        return;
    }

    _signal_auto_update = false;  // deactivate signal auto update since we do not want intermediate callbacks to trigger new signals

    for (int i = 0; i < num_fields; ++i)
    {
        const p::FieldDescriptor* field = desc->field(i);

        if (field->cpp_type() != p::FieldDescriptor::CPPTYPE_MESSAGE)
        {
            PRINT_ERROR("ToolboxWidget::fromMessage(): non-message-field found in parameter message. Only message fields are allowed in top-level");
            return;
        }
        _param_widgets[i]->generateFromAllocatedField(_param_message.get(), field);
    }

    _signal_auto_update = true;  // reactivate signal auto update

    // now also get available signals
    getAvailableSignals();
}

void ToolboxWidget::requestParametersFromService()
{
    if (!_rpc_client)
    {
        clearParameters();
        return;
    }

    if (_rpc_client->getParameters(*_param_message))
    {
        updateParameterWidgets();
    }
}

bool ToolboxWidget::sendParametersToService(bool suppress_warning)
{
    bool ret_val = true;
    QString status;

    if (!_rpc_client)
    {
        status += "Cannot connect to RPC client / connection lost\n";
        if (!suppress_warning) QMessageBox::warning(this, tr("RPC Connection Error"), status);
        return false;
    }

    messages::Status status_msg;
    if (_rpc_client->setParameters(*_param_message, status_msg))
    {
        if (!status_msg.ok())
        {
            status += QString::fromStdString(status_msg.text()) + "\n";
            ret_val = false;
        }
    }

    if (ret_val == false && !suppress_warning) QMessageBox::warning(this, tr("Invalid Configuration"), status);
    return ret_val;
}

void ToolboxWidget::performTask()
{
    if (!_rpc_client) return;

    if (_task_running)
    {
        if (_rpc_client->stopTask()) setTaskRunningFlag(false);
        return;
    }

    // first update parameters
    if (!sendParametersToService()) return;

    // perform actual task
    setTaskRunningFlag(true);

    RpcTaskWorker* rpc_worker = new RpcTaskWorker;
    rpc_worker->setRpcClient(_rpc_client);
    rpc_worker->moveToThread(&_rpc_task_thread);
    connect(&_rpc_task_thread, SIGNAL(started()), rpc_worker, SLOT(performTask()));
    connect(rpc_worker, &RpcTaskWorker::feedback, this, [this](const messages::Signal& signal) { _signal_helper->addSignal(signal); });
    connect(rpc_worker, SIGNAL(taskFinished(bool, const QString&)), this, SLOT(rpcTaskFinished(bool, const QString&)));
    connect(&_rpc_task_thread, SIGNAL(finished()), rpc_worker, SLOT(deleteLater()));

    _rpc_task_thread.start(QThread::NormalPriority);
}

void ToolboxWidget::rpcTaskFinished(bool success, const QString& msg)
{
    setTaskRunningFlag(false);

    _rpc_task_thread.quit();  // stop thread's event-loop

    if (success)
    {
        emit taskCompleted(QDateTime::currentDateTime());

        // add parameters to cache
        // TODO(roesmann) avoid copying by directly writing form ParameterWidgets to _param_message
        // toMessage(*_param_message);
        _task_cache.toCache(std::to_string(_signal_helper->currentSeriesId()), *_param_message);
        // TODO(roesmann) we create a copy in both toMessage() and toCache() again: place for improving efficiency

        // increase the id for the next run
        _signal_helper->startNewSeries();
    }
    else
    {
        QMessageBox::warning(this, tr("Something went wrong..."), msg);
        return;
    }

    // request parameters for next task
    if (_signal_auto_update) getAvailableSignals();
}

void ToolboxWidget::getAvailableSignals()
{
    if (!_signal_auto_update) return;  // TODO(roesmann) quick workaround: we do not want to update signals at all if signal update is deactivated

    if (!_rpc_client) return;

    // first update parameters
    sendParametersToService(true);

    // get current signals
    _signal_helper->clearCurrentSeries();  // TODO(roesmann) instead of clearing everything, check which signals are updated and erase only
                                           // the other ones!
    auto signal_feedback = [this](const messages::Signal& signal) { _signal_helper->addSignal(signal); };
    _rpc_client->getAvailableSignals(signal_feedback);
}

void ToolboxWidget::restoreTaskParameters(int task_id)
{
    std::unique_ptr<google::protobuf::Message> cached_msg = _task_cache.fromCache(std::to_string(task_id));
    if (cached_msg)
    {
        fromMessage(*cached_msg);
    }
}

void ToolboxWidget::removeFromTaskCache(int task_id) { _task_cache.erase(std::to_string(task_id)); }

void ToolboxWidget::saveCurrentParametersToFile()
{
    if (!_param_message) return;
    saveParametersToFile(*_param_message);
}

void ToolboxWidget::saveParametersToFile(int task_id)
{
    std::unique_ptr<google::protobuf::Message> cached_msg = _task_cache.fromCache(std::to_string(task_id));
    if (cached_msg)
        saveParametersToFile(*cached_msg);
    else
        QMessageBox::warning(this, "Warning", "Cannot save parameters to file, since not found in cache.");
}

void ToolboxWidget::saveParametersToFile(const google::protobuf::Message& params)
{
    QString default_filter = tr("Parameter file (*.cparams)");
#ifdef __linux__
    QString filename = QFileDialog::getSaveFileName(this, tr("Export Parameters"), "parameters.cparams",
                                                    tr("Parameter file (*.cparams);;All files (*.*)"), &default_filter);
#else
    QString filename =
        QFileDialog::getSaveFileName(this, tr("Export Parameters"), "parameters", tr("Parameter file (*.cparams);;All files (*.*)"), &default_filter);
#endif
    if (filename.isEmpty()) return;

    // QFile file(filename);
    std::ofstream file(filename.toStdString(), std::ofstream::binary | std::ofstream::trunc);

    if (!file.is_open())
    {
        QMessageBox::warning(this, "Warning", "Cannot open file with write-permissions for exporting parameters.");
        return;
    }
    params.SerializeToOstream(&file);

    emit parameterSavedToFile(filename);
}

void ToolboxWidget::loadParametersFromFile()
{
    QString default_filter = tr("Parameter file (*.cparams)");
    QString filename =
        QFileDialog::getOpenFileName(this, tr("Import Parameters"), "parameters", tr("Parameter file (*.cparams);;All files (*.*)"), &default_filter);
    loadParametersFromFile(filename);
}

void ToolboxWidget::loadParametersFromFile(const QString& filepath)
{
    if (!_param_message)
    {
        QMessageBox::warning(this, "Warning", "Cannot load parameter before RPC client is configured.");
        return;
    }

    if (filepath.isEmpty()) return;

    std::ifstream file(filepath.toStdString());

    if (!file.is_open())
    {
        QMessageBox::warning(this, "Warning", "Cannot open file with read-permissions for importing parameters.");
        return;
    }

    _param_message->Clear();
    if (_param_message->ParsePartialFromIstream(&file))
    {
        updateParameterWidgets();
        emit parameterLoadedFromFile(filepath);
    }
}

void ToolboxWidget::clearParameters()
{
    for (ParameterWidget* param_widget : _param_widgets) param_widget->clearElements();
}

void ToolboxWidget::setRpcClient(std::shared_ptr<MasterServiceClient> rpc_client)
{
    _rpc_client = rpc_client;
    if (_rpc_client)
    {
        _param_message = rpc_client->createParameterMsg();
        initializeToolboxes();
        requestParametersFromService();
    }
}

void ToolboxWidget::fromMessage(const google::protobuf::Message& parameters)
{
    if (!_param_message) return;

    _param_message->CopyFrom(parameters);
    updateParameterWidgets();
}

void ToolboxWidget::toMessage(google::protobuf::Message& parameters) const
{
    if (!_param_message) return;

    parameters.CopyFrom(*_param_message);
}

void ToolboxWidget::setTaskRunningFlag(bool running)
{
    _task_running = running;
    if (running)
        _btn_perform_task->setText(_task_button_label_inactive);
    else
        _btn_perform_task->setText(_task_button_label_active);
}

}  // namespace gui
}  // namespace corbo
