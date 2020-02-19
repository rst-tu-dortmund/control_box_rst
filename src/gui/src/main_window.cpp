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

#include <corbo-gui/main_window.h>

#include <QCoreApplication>
#include <QDir>
#include <QInputDialog>
#include <QMenuBar>
#include <QMessageBox>
#include <QScrollArea>
#include <QSettings>
#include <QStatusBar>
#include <QTimer>
#include <QVBoxLayout>

#include <memory>
#include <utility>

namespace corbo {
namespace gui {

corboMainWindow::corboMainWindow(QWidget* parent)
{
    // set working path of qdir to path of the executable:
    QDir::setCurrent(QCoreApplication::applicationDirPath());

    _central_widget = new QWidget;
    setCentralWidget(_central_widget);

    setWindowTitle(tr("Control-Box RST"));

    // resize(800,500);
    setWindowState(Qt::WindowMaximized);

    // TODO(roesmann): this should be handled in a separate thread with low priority
    // create signal helper which stores all available signals
    _signal_helper = std::make_shared<SignalHelper>();

    _rpc_client = std::make_shared<MasterServiceClient>();

    createDockWindows();
    createActions();
    createMenus();
    createStatusBar();

    addMenuActions();

    setDefaultLayout();

    // read window settings from previous runs if available
    readWindowSettings();

    configureRpcWorker();

    // register meta types for queued signal-slot connections
    qRegisterMetaType<std::shared_ptr<MasterServiceClient>>("std::shared_ptr<MasterServiceClient>");
}

corboMainWindow::corboMainWindow(std::shared_ptr<MasterServiceClient> rpc_client, QWidget* parent) : corboMainWindow(parent)
{
    if (_rpc_client) _rpc_client = rpc_client;
}

corboMainWindow::~corboMainWindow()
{
    _rpc_worker_thread.quit();
    _rpc_worker_thread.wait();
}

void corboMainWindow::closeEvent(QCloseEvent* event) { writeWindowSettings(); }

void corboMainWindow::createActions()
{
    // save parameters to file
    _action_save_parameters_file = new QAction(tr("Save Parameters..."), this);
    _action_save_parameters_file->setStatusTip(tr("Save current parameters to file."));
    connect(_action_save_parameters_file, &QAction::triggered, _toolbox_widget,
            [this](bool /*checked*/) { _toolbox_widget->saveCurrentParametersToFile(); });

    // load parameters from file
    _action_load_parameters_file = new QAction(tr("Load Parameters..."), this);
    _action_load_parameters_file->setStatusTip(tr("Load parameters from file."));
    connect(_action_load_parameters_file, &QAction::triggered, _toolbox_widget,
            [this](bool /*checked*/) { _toolbox_widget->loadParametersFromFile(); });

    // Clear scopes
    _action_remove_signals = new QAction(tr("Clear Signals"), this);
    _action_remove_signals->setStatusTip(tr("Clear all signals"));
    connect(_action_remove_signals, &QAction::triggered, [this](bool /*checked*/) {
        _signal_helper->clearAll();
        _toolbox_widget->getAvailableSignals();  // TODO(roesmann) should we better clear only recent and not current signals?
    });

    // Close all scopes
    _action_close_scopes = new QAction(tr("Close Scopes"), this);
    _action_close_scopes->setStatusTip(tr("Clear all scopes"));
    connect(_action_close_scopes, SIGNAL(triggered()), _scope_collection_widget, SLOT(closeAllScopes()));

    // rescale scope axes
    _action_rescale_scope_axes = new QAction(tr("Rescale Scope Axes"), this);
    connect(_action_rescale_scope_axes, &QAction::triggered, _scope_collection_widget, &ScopeCollectionWidget::resizeScopeAxes);

    // Reset layout to default
    _action_default_layout = new QAction(tr("Default Layout"), this);
    _action_default_layout->setStatusTip(tr("Set docks to its default position"));
    connect(_action_default_layout, SIGNAL(triggered()), this, SLOT(setDefaultLayout()));

    // rpc server address
    _action_rpc_server_address = new QAction(tr("Change Server Address..."), this);
    connect(_action_rpc_server_address, &QAction::triggered, this, [this](bool /*checked*/) {
        bool ok;
        QString text = QInputDialog::getText(this, tr("QInputDialog::getText()"), tr("Server Address:"), QLineEdit::Normal, _rpc_server_address, &ok);
        if (ok && !text.isEmpty())
        {
            _rpc_server_address = text;
            emit requestRpcServerAddress(_rpc_server_address);
        }
    });

    // About dialog
    _action_about = new QAction(tr("About"), this);
    _action_about->setStatusTip(tr("Display general information about the program"));
    connect(_action_about, SIGNAL(triggered()), this, SLOT(aboutPopup()));

    // Exit action
    _action_exit_program = new QAction(tr("Exit"), this);
    connect(_action_exit_program, SIGNAL(triggered()), this, SLOT(close()));

    // Recent files
    for (std::size_t i = 0; i < _action_recent_files.size(); ++i)
    {
        _action_recent_files[i] = new QAction(this);
        _action_recent_files[i]->setVisible(false);
        QAction* cur_action = _action_recent_files[i];
        connect(_action_recent_files[i], &QAction::triggered, _toolbox_widget, [this, cur_action]() {
            if (cur_action) _toolbox_widget->loadParametersFromFile(cur_action->data().toString());
        });
    }
}

void corboMainWindow::createMenus()
{
    // File menu
    _menu_file = menuBar()->addMenu(tr("File"));

    // Param menu
    _menu_parameters = menuBar()->addMenu(tr("Parameters"));

    // Scopes menu
    _menu_scopes = menuBar()->addMenu(tr("Scopes"));

    // Signals menu
    _menu_signals = menuBar()->addMenu(tr("Signals"));

    // View menu
    _menu_view = menuBar()->addMenu(tr("View"));

    // Help menu
    _menu_help = menuBar()->addMenu(tr("&Help"));

#ifdef __linux__
    menuBar()->setNativeMenuBar(false);  // Otherwise we do not see a menu bar on Ubuntu 14.04
#endif
}

void corboMainWindow::createStatusBar() { displayStatusReady(); }

void corboMainWindow::createDockWindows()
{
    // Toolbox dock
    QScrollArea* toolbox_scroll_area = new QScrollArea;
    _toolbox_widget                  = new ToolboxWidget(_signal_helper, _rpc_client);
    toolbox_scroll_area->setWidget(_toolbox_widget);
    connect(_toolbox_widget, &ToolboxWidget::parameterSavedToFile, this, &corboMainWindow::addToRecentFileAction);
    connect(_toolbox_widget, &ToolboxWidget::parameterLoadedFromFile, this, &corboMainWindow::addToRecentFileAction);

    _dock_toolbox_widget = new QDockWidget(tr("Toolbox"), this);
    _dock_toolbox_widget->setFloating(false);
    _dock_toolbox_widget->setObjectName("ToolboxDock");
    _dock_toolbox_widget->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);
    _dock_toolbox_widget->setWidget(_toolbox_widget);
    _dock_toolbox_widget->setFeatures(QDockWidget::AllDockWidgetFeatures);
    addDockWidget(Qt::LeftDockWidgetArea, _dock_toolbox_widget);

    // Scope dock
    // QScrollArea* scope_scroll_area = new QScrollArea;
    _scope_collection_widget = new ScopeCollectionWidget(_signal_helper);
    // scope_scroll_area->setWidget(_scope_collection_widget);
    connect(_signal_helper.get(), &SignalHelper::newMeasurement, _scope_collection_widget, &ScopeCollectionWidget::addMeasurement);
    connect(_signal_helper.get(), &SignalHelper::signalRemoved, _scope_collection_widget, &ScopeCollectionWidget::removeSignal);
    connect(_signal_helper.get(), &SignalHelper::newSeries, _scope_collection_widget, &ScopeCollectionWidget::initializeTask);

    _dock_scope_widget = new QDockWidget(tr("Scopes"), this);
    _dock_scope_widget->setFloating(false);
    _dock_scope_widget->setObjectName("ScopesDock");
    _dock_scope_widget->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);
    // _scope_collection_widget->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    // _dock_scope_widget->setWidget(scope_scroll_area);
    _dock_scope_widget->setWidget(_scope_collection_widget);
    _dock_scope_widget->setFeatures(QDockWidget::AllDockWidgetFeatures);
    addDockWidget(Qt::RightDockWidgetArea, _dock_scope_widget);

    // Signals dock
    QScrollArea* signal_scroll_area = new QScrollArea;

    // pass the signal map to the scope widget such that it can lookup signal keys and receive the corresponding object
    // note, we cannot pass signal objects directly via drag-and-drop from the signal-widget, but we can drop the key.
    _signal_collection_widget = new SignalCollectionWidget(_signal_helper);
    connect(_signal_helper.get(), &SignalHelper::newSignal, _signal_collection_widget, &SignalCollectionWidget::addSignal);
    connect(_signal_helper.get(), &SignalHelper::signalRemoved, _signal_collection_widget, &SignalCollectionWidget::removeTreeItem);
    connect(_toolbox_widget, &ToolboxWidget::taskCompleted, _signal_collection_widget, &SignalCollectionWidget::moveToRecent);
    connect(_signal_collection_widget, &SignalCollectionWidget::requestTaskParameterBackup, _toolbox_widget, &ToolboxWidget::restoreTaskParameters);
    connect(_signal_collection_widget, &SignalCollectionWidget::taskRemoved, _toolbox_widget, &ToolboxWidget::removeFromTaskCache);

    signal_scroll_area->setWidget(_signal_collection_widget);
    signal_scroll_area->setWidgetResizable(true);

    _dock_signal_widget = new QDockWidget(tr("Signals"), this);
    _dock_signal_widget->setFloating(false);
    _dock_signal_widget->setObjectName("SignalsDock");
    _dock_signal_widget->setAllowedAreas(Qt::AllDockWidgetAreas);
    _dock_signal_widget->setWidget(signal_scroll_area);
    _dock_signal_widget->setFeatures(QDockWidget::AllDockWidgetFeatures);

    // update available signals
    _toolbox_widget->getAvailableSignals();
    _toolbox_widget->setSignalAutoRequest(true);
}

void corboMainWindow::addMenuActions()
{
    _menu_file->addAction(_action_rpc_server_address);
    _menu_file->addSeparator();
    _menu_file->addAction(_action_exit_program);

    _menu_view->addAction(_action_default_layout);

    _menu_parameters->addAction(_action_save_parameters_file);
    _menu_parameters->addSeparator();
    _menu_parameters->addAction(_action_load_parameters_file);
    _menu_recent_parameters = _menu_parameters->addMenu(tr("Recent Parameters"));
    for (std::size_t i = 0; i < _action_recent_files.size(); ++i) _menu_recent_parameters->addAction(_action_recent_files[i]);
    updateRecentFileActions();

    _menu_scopes->addAction(_action_rescale_scope_axes);
    _menu_scopes->addAction(_action_close_scopes);

    _menu_signals->addAction(_action_remove_signals);

    _menu_help->addAction(_action_about);
}

void corboMainWindow::setDefaultLayout()
{
    _dock_toolbox_widget->hide();
    _dock_scope_widget->hide();
    _dock_signal_widget->hide();
    addDockWidget(Qt::LeftDockWidgetArea, _dock_toolbox_widget);
    _dock_toolbox_widget->setFloating(false);
    _dock_toolbox_widget->show();

    addDockWidget(Qt::RightDockWidgetArea, _dock_scope_widget);
    _dock_scope_widget->setFloating(false);
    _dock_scope_widget->show();

    addDockWidget(Qt::RightDockWidgetArea, _dock_signal_widget);
    splitDockWidget(_dock_scope_widget, _dock_signal_widget, Qt::Horizontal);
    _dock_signal_widget->setFloating(false);
    _dock_signal_widget->show();

    // bring scope window to front
    // _dock_scope_widget->raise();
}

void corboMainWindow::readWindowSettings()
{
    QSettings settings(OrganizationName, ApplicationName);
    QVariant pos_var                = settings.value("pos");
    QVariant size_var               = settings.value("size");
    QVariant dock_var               = settings.value("dock_locations");
    QVariant rpc_server_address_var = settings.value("rpc_server_address");

    if (pos_var.canConvert<QPoint>()) move(pos_var.toPoint());
    if (size_var.canConvert<QSize>()) resize(size_var.toSize());
    if (dock_var.canConvert<QByteArray>()) restoreState(dock_var.toByteArray());
    if (rpc_server_address_var.canConvert<QString>()) _rpc_server_address = rpc_server_address_var.toString();
}

void corboMainWindow::writeWindowSettings()
{
    QSettings settings(OrganizationName, ApplicationName);
    settings.setValue("pos", pos());
    settings.setValue("size", size());
    settings.setValue("dock_locations", saveState());
    settings.setValue("rpc_server_address", _rpc_server_address);
}

void corboMainWindow::displayStatus(const QString& status, double duration_ms)
{
    if (status.isEmpty())
    {
        displayStatusReady();
        return;
    }
    statusBar()->showMessage(status);

    if (duration_ms > 0)
    {
        // start timer to reset to ready after the specified number of miliseconds
        QTimer::singleShot(duration_ms, this, SLOT(displayStatusReady()));
    }
}

void corboMainWindow::updateRecentFileActions()
{
    QSettings settings(OrganizationName, ApplicationName);
    QStringList files = settings.value("recentFileList").toStringList();

    int num_recent_files = qMin(files.size(), (int)_action_recent_files.size());

    for (int i = 0; i < num_recent_files; ++i)
    {
        QString text = files[i];
        _action_recent_files[i]->setText(text);
        _action_recent_files[i]->setData(files[i]);
        _action_recent_files[i]->setVisible(true);
    }
    for (int j = num_recent_files; j < (int)_action_recent_files.size(); ++j) _action_recent_files[j]->setVisible(false);
}

void corboMainWindow::addToRecentFileAction(const QString& filename)
{
    QSettings settings(OrganizationName, ApplicationName);
    QStringList files = settings.value("recentFileList").toStringList();
    files.removeAll(filename);
    files.prepend(filename);
    while (files.size() > _action_recent_files.size()) files.removeLast();

    settings.setValue("recentFileList", files);

    updateRecentFileActions();
}

void corboMainWindow::aboutPopup()
{
    QMessageBox::about(this, tr("About"), tr("<p><b>Institute of Control Theory and Systems Engineering, TU Dortmund University.</b></p>"
                                             "<p>Written by Christoph Rösmann, <a "
                                             "href=\"mailto:christoph.roesmann@tu-dortmund.de\">christoph.roesmann@tu-dortmund.de</a></p>"));
}

void corboMainWindow::configureRpcWorker()
{
    displayStatus("trying to connect with rpc server at " + _rpc_server_address);

    RpcConnectionWorker* rpc_worker = new RpcConnectionWorker(_rpc_client);
    rpc_worker->moveToThread(&_rpc_worker_thread);
    connect(this, SIGNAL(requestRpcServerAddress(const QString&)), rpc_worker, SLOT(connectRpc(const QString&)));
    connect(&_rpc_worker_thread, SIGNAL(finished()), rpc_worker, SLOT(deleteLater()));
    connect(rpc_worker, SIGNAL(connectionResult(std::shared_ptr<MasterServiceClient>)), this,
            SLOT(updatedRpcClient(std::shared_ptr<MasterServiceClient>)));
    connect(rpc_worker, &RpcConnectionWorker::connectionResult, this,
            [this](std::shared_ptr<MasterServiceClient> /*client*/) { displayStatusReady(); });

    _rpc_worker_thread.start(QThread::LowPriority);

    emit requestRpcServerAddress(_rpc_server_address);
}

void corboMainWindow::updatedRpcClient(std::shared_ptr<MasterServiceClient> client)
{
    if (client)
    {
        // _rpc_client = client;
        // _toolbox_widget->setRpcClient(_rpc_client);
        _toolbox_widget->requestParametersFromService();
    }
}

}  // namespace gui
}  // namespace corbo
