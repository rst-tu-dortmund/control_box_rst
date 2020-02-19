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

#ifndef SRC_GUI_INCLUDE_CORBO_GUI_SIGNAL_COLLECTION_WIDGET_H_
#define SRC_GUI_INCLUDE_CORBO_GUI_SIGNAL_COLLECTION_WIDGET_H_

#include <corbo-gui/collapsable_groupbox.h>
#include <corbo-gui/extended_tree_widget.h>
#include <corbo-gui/signal_helper.h>
#include <corbo-gui/signal_widget.h>
#include <QDateTime>
#include <QFileDialog>
#include <QHash>
#include <QList>
#include <QMenu>
#include <QString>
#include <QTreeWidget>
#include <QVBoxLayout>
#include <QWidget>

#include <tuple>

namespace corbo {
namespace gui {

class SignalCollectionWidget : public QWidget
{
    Q_OBJECT

 public:
    explicit SignalCollectionWidget(SignalHelper::Ptr signal_helper, QWidget* parent = 0);
    virtual ~SignalCollectionWidget();

    QSize sizeHint() const override;

 public slots:
    void addSignal(const QString& key, SignalHelper::SignalData& signal_data);
    void moveToRecent(const QDateTime& dtime);
    void removeTreeItem(const QString& key, int value_idx);

 signals:
    void requestTaskParameterBackup(int task_id);
    void taskRemoved(int task_id);

 protected:
    ExtendedTreeWidget* createSignalTree();
    void createHistoryArea();

    void getGroupInfo(const QTreeWidget* tree, const QTreeWidgetItem* item, QString* name_out, int* id_out);
    void removeSignal(QTreeWidget* tree, QTreeWidgetItem* item);

    void namespaceContextMenu(const QPoint& point, QTreeWidget* tree);
    void recentGroupContextMenu(const QPoint& point, CollapsableGroupBox* group, int task_id);

    void addContextActionsSingleValueSignal(SignalWidget& signal_widget, QMenu& menu);
    void addContextActionsMultiValueSignal(SignalInterface::Ptr signal, QMenu& menu, QMenu& menu_export);

    void addContextActionsTimeSeriesSignal(TimeSeriesSignal::Ptr ts_signal, QMenu& menu, int signal_value_idx, QMenu* menu_export = nullptr);
    void addContextActionsTimeSeriesSignal(TimeSeriesSignal::Ptr ts_signal, QMenu& menu, QMenu* menu_export = nullptr);
    void addContextActionsTimeSeriesSequenceSignal(TimeSeriesSequenceSignal::Ptr ts_seq_signal, QMenu& menu, int signal_value_idx,
                                                   QMenu* menu_export = nullptr);
    void addContextActionsTimeSeriesSequenceSignal(TimeSeriesSequenceSignal::Ptr ts_seq_signal, QMenu& menu, QMenu* menu_export = nullptr);
    void addContextActionsIndexedValuesSetSignal(IndexedValuesSetSignal::Ptr ivs_signal, QMenu& menu, QMenu* menu_export = nullptr);

    void addContextActionsMatrixSignal(MatrixSignal::Ptr mat_signal, QMenu& menu, QMenu* menu_export = nullptr);
    void addContextActionsMatrixSetSignal(MatrixSetSignal::Ptr mat_set_signal, QMenu& menu, QMenu* menu_export = nullptr);

 private:
    SignalHelper::Ptr _signal_helper;

    QVBoxLayout* _main_layout;

    QVector<std::tuple<CollapsableGroupBox*, QWidget*, int>> _recent_signals;  // groupbox + extended_tree_widget + task_id
    ExtendedTreeWidget* _signal_tree;
    int _current_task_id = 0;
};

}  // namespace gui
}  // namespace corbo

#endif  // SRC_GUI_INCLUDE_CORBO_GUI_SIGNAL_COLLECTION_WIDGET_H_
