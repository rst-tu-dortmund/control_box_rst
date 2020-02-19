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

#ifndef SRC_GUI_INCLUDE_CORBO_GUI_MAIN_WINDOW_H_
#define SRC_GUI_INCLUDE_CORBO_GUI_MAIN_WINDOW_H_

#include <corbo-communication/main_service_client.h>
#include <corbo-gui/rpc_connection_worker.h>
#include <corbo-gui/scope_collection_widget.h>
#include <corbo-gui/signal_collection_widget.h>
#include <corbo-gui/signal_helper.h>
#include <corbo-gui/toolbox_widget.h>

#include <QAction>
#include <QDockWidget>
#include <QMainWindow>
#include <QMenu>
#include <QThread>

#include <array>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace corbo {
namespace gui {

const QString OrganizationName = "rst";
const QString ApplicationName  = "control-box-rst";

/**
 * @brief GUI main window
 *
 * @ingroup gui
 *
 * The gui includes an instance of the MasterServiceClient
 * in order to communicate with a Master service.
 * Communication is performed via RPC and Protobuf messages
 *
 * @see Master MasterServiceClient
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 *
 * @todo GUI code documentation
 */
class corboMainWindow : public QMainWindow
{
    Q_OBJECT

 public:
    explicit corboMainWindow(QWidget* parent = 0);
    explicit corboMainWindow(std::shared_ptr<MasterServiceClient> rpc_client, QWidget* parent = 0);
    virtual ~corboMainWindow();

 signals:
    void requestRpcServerAddress(const QString& address);

 protected:
    void closeEvent(QCloseEvent* event) override;

    void createActions();
    void createMenus();
    void createStatusBar();
    void createDockWindows();
    void addMenuActions();
    void updateRecentFileActions();

    void configureRpcWorker();

    void readWindowSettings();
    void writeWindowSettings();

 protected slots:
    void setDefaultLayout();
    void addToRecentFileAction(const QString& filename);

 private slots:
    void displayStatus(const QString& status, double duration_ms = 0);
    void displayStatusReady() { displayStatus(tr("Ready")); }
    void aboutPopup();
    void updatedRpcClient(std::shared_ptr<MasterServiceClient> client);

 private:
    SignalHelper::Ptr _signal_helper;
    std::shared_ptr<MasterServiceClient> _rpc_client;

    QWidget* _central_widget;

    QAction* _action_save_parameters_file;
    QAction* _action_load_parameters_file;
    QAction* _action_remove_signals;
    QAction* _action_close_scopes;
    QAction* _action_rescale_scope_axes;
    QAction* _action_default_layout;
    QAction* _action_rpc_server_address;
    QAction* _action_exit_program;
    QAction* _action_about;
    std::array<QAction*, 6> _action_recent_files;  // number of recent files

    QMenu* _menu_file;
    QMenu* _menu_parameters;
    QMenu* _menu_recent_parameters;
    QMenu* _menu_scopes;
    QMenu* _menu_signals;
    QMenu* _menu_view;
    QMenu* _menu_help;

    QDockWidget* _dock_toolbox_widget;
    QDockWidget* _dock_signal_widget;
    QDockWidget* _dock_scope_widget;

    ToolboxWidget* _toolbox_widget;
    SignalCollectionWidget* _signal_collection_widget;
    ScopeCollectionWidget* _scope_collection_widget;

    QThread _rpc_worker_thread;
    QString _rpc_server_address = "localhost:50051";
};

}  // namespace gui
}  // namespace corbo

#endif  // SRC_GUI_INCLUDE_CORBO_GUI_MAIN_WINDOW_H_
