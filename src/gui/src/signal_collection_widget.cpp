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

#include <corbo-core/corbo_core.h>
#include <corbo-gui/signal_collection_widget.h>
#include <corbo-gui/utilities.h>
#include <QStack>

#include <QLabel>
#include <QMessageBox>
#include <QMouseEvent>
#include <QTreeWidget>

#include <tuple>

namespace corbo {
namespace gui {

SignalCollectionWidget::SignalCollectionWidget(SignalHelper::Ptr signal_helper, QWidget* parent) : QWidget(parent)
{
    setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Expanding);

    _signal_helper = signal_helper;

    _main_layout = new QVBoxLayout(this);
    _main_layout->setAlignment(Qt::AlignTop | Qt::AlignLeft);
    _main_layout->setContentsMargins(0, 0, 0, 0);

    QHBoxLayout* title_layout = new QHBoxLayout;
    title_layout->setAlignment(Qt::AlignTop | Qt::AlignLeft);
    title_layout->setContentsMargins(5, 5, 5, 5);
    title_layout->addWidget(new QLabel(tr("Available Signals")));
    title_layout->addWidget(util::create_horizontal_line());
    _main_layout->addLayout(title_layout);

    _signal_tree = createSignalTree();
    _main_layout->addWidget(_signal_tree);

    _main_layout->addWidget(util::create_horizontal_line());
    _main_layout->addSpacing(10);
}

SignalCollectionWidget::~SignalCollectionWidget() {}

QSize SignalCollectionWidget::sizeHint() const { return QSize(200, 800); }

ExtendedTreeWidget* SignalCollectionWidget::createSignalTree()
{
    ExtendedTreeWidget* signal_tree = new ExtendedTreeWidget;
    signal_tree->viewport()->setBackgroundRole(QPalette::Background);
    signal_tree->setAttribute(Qt::WA_MacShowFocusRect, false);
    // signal_tree->setRootIsDecorated(false);
    // signal_tree->setColumnCount(1);
    // signal_tree->setIconSize(QSize(20, 20));
    signal_tree->setHeaderHidden(true);
    signal_tree->setIndentation(10);
    signal_tree->setStyleSheet("QTreeView { border: none; }");
    signal_tree->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);
    signal_tree->setContextMenuPolicy(Qt::CustomContextMenu);
    signal_tree->connect(signal_tree, &QTreeWidget::customContextMenuRequested,
                         [this, signal_tree](const QPoint& point) { this->namespaceContextMenu(point, signal_tree); });
    return signal_tree;
}

void SignalCollectionWidget::addSignal(const QString& key, SignalHelper::SignalData& signal_data)
{
    _current_task_id = _signal_helper->currentSeriesId();

    QStringList namespaces = signal_data.namespaces;
    namespaces.pop_back();  // remove own name

    QVector<QTreeWidgetItem*> namespace_items;

    for (const QString& group : namespaces)
    {
        // check if we have the current namespace
        QTreeWidgetItem* current_item = nullptr;
        if (namespace_items.empty())  // check top-level
        {
            QList<QTreeWidgetItem*> items    = _signal_tree->findItems(group, Qt::MatchExactly);
            if (!items.empty()) current_item = items.front();
        }
        else
        {
            // check children of last group
            for (int i = 0; i < namespace_items.back()->childCount(); ++i)
            {
                if (namespace_items.back()->child(i)->text(0) == group)
                {
                    current_item = namespace_items.back()->child(i);
                }
            }
        }

        if (!current_item)
        {
            // namespace not yet created
            QTreeWidgetItem* item = new QTreeWidgetItem;
            if (namespace_items.empty())
                _signal_tree->addTopLevelItem(item);  // new top-level item
            else
                namespace_items.back()->addChild(item);  // add as child to previous item
            item->setText(0, group);
            namespace_items.push_back(item);
        }
        else
        {
            // append item to the current item list
            namespace_items.push_back(current_item);
        }
    }

    // now add signal
    // first check if we have multiple dimensions
    int dimension = signal_data.signal->header.value_dimension;

    if (dimension == 0) return;

    QTreeWidgetItem* item = new QTreeWidgetItem;
    if (namespace_items.empty())  // new top-level item
        _signal_tree->addTopLevelItem(item);
    else
        namespace_items.back()->addChild(item);

    if (dimension == 1)
    {
        // add signal name directly
        SignalWidget* signal_widget = new SignalWidget(signal_data.short_name, key, 0);
        signal_data.active_indices.insert(0);
        _signal_tree->setItemWidget(item, 0, signal_widget);
    }
    else
    {
        // treat signal name as namespace and add 1D-subsignals
        item->setText(0, signal_data.namespaces.back());

        for (int i = 0; i < dimension; ++i)
        {
            // create 1D-subsignals
            QTreeWidgetItem* subitem = new QTreeWidgetItem;
            item->addChild(subitem);

            // get custom names if available
            QString sublabel = (i < signal_data.value_labels.size()) ? signal_data.value_labels[i] : signal_data.short_name + QString::number(i);
            SignalWidget* signal_widget = new SignalWidget(sublabel, key, i);
            signal_data.active_indices.insert(i);
            _signal_tree->setItemWidget(subitem, 0, signal_widget);
        }
    }

    _signal_tree->updateSizeHint();
}

void SignalCollectionWidget::moveToRecent(const QDateTime& dtime)
{
    ExtendedTreeWidget* new_tree = createSignalTree();
    QLayoutItem* current_tree    = _main_layout->replaceWidget(_signal_tree, new_tree, Qt::FindDirectChildrenOnly);
    _signal_tree                 = new_tree;

    // update name in common signal target
    if (_signal_helper)
    {
        QString groupname = "data_" + dtime.toString("yyyy_MM_dd_hh_mm_ss");  // prepend text, e.g., 'data', because otherwise some yaml importer fail
        _signal_helper->getSignalTreeRef()[_current_task_id].setTopLevelGroupName(groupname.toStdString());
    }

    if (current_tree)
    {
        QVBoxLayout* layout = new QVBoxLayout;
        layout->setContentsMargins(15, 0, 0, 0);
        layout->addItem(current_tree);

        QString title              = dtime.toString("dd.MM.yyyy hh:mm:ss");
        CollapsableGroupBox* group = new CollapsableGroupBox(title);

        // add context menu
        group->setContextMenuPolicy(Qt::CustomContextMenu);
        int task_id = _current_task_id;
        group->connect(group, &QTreeWidget::customContextMenuRequested,
                       [this, group, task_id](const QPoint& point) { this->recentGroupContextMenu(point, group, task_id); });

        group->groupBox()->setLayout(layout);

        _recent_signals.push_back(std::make_tuple(group, current_tree->widget(), _current_task_id));
        _main_layout->addWidget(group);

        ++_current_task_id;  // will also be update next time signals are added
    }
}

void SignalCollectionWidget::removeTreeItem(const QString& key, int value_idx)
{
    QString name;
    int id = 0;
    SignalHelper::key2Name(key, name, id);

    CollapsableGroupBox* group_box = nullptr;
    QTreeWidget* tree              = nullptr;
    auto group_it                  = _recent_signals.end();

    // search tree and group for current id;
    if (id == _current_task_id)
    {
        tree = _signal_tree;
    }
    else
    {
        group_it = std::find_if(_recent_signals.begin(), _recent_signals.end(),
                                [id](const std::tuple<CollapsableGroupBox*, QWidget*, int> item) { return id == std::get<2>(item); });
        if (group_it != _recent_signals.end())
        {
            group_box = std::get<0>(*group_it);
            tree      = static_cast<QTreeWidget*>(std::get<1>(*group_it));
        }
    }

    if (!tree) return;

    // search for current item:
    QStringList namespaces;
    SignalHelper::extractNamespace(name, namespaces);

    // find top-level child;
    QTreeWidgetItem* item      = nullptr;
    QTreeWidgetItem* candidate = nullptr;
    bool found                 = false;

    // n -> size+1 : also expend last namespace (signal name) in case of multi-value signals
    for (int n = 0; n < namespaces.size() + 1; ++n)
    {
        if (found) break;

        int no_items;
        if (n == 0)
            no_items = tree->topLevelItemCount();
        else if (item)
            no_items = item->childCount();
        else
            break;

        for (int i = 0; i < no_items; ++i)
        {
            QTreeWidgetItem* child = (n == 0) ? tree->topLevelItem(i) : item->child(i);

            // check leaf widget
            SignalWidget* signal_widget = dynamic_cast<SignalWidget*>(tree->itemWidget(child, 0));
            if (signal_widget)
            {
                if (signal_widget->key().compare(key) == 0)
                {
                    if (signal_widget->valueIdx() == value_idx || value_idx == SignalHelper::ALL_VALUES)
                    {
                        found     = true;
                        candidate = child;
                        break;
                    }
                }
            }

            // check namespace
            if (n < namespaces.size() && child->text(0).compare(namespaces[n]) == 0)
            {
                if (namespaces.size() == n + 1)
                {
                    if (value_idx == SignalHelper::ALL_VALUES)
                    {
                        found     = true;
                        candidate = child;
                        break;
                    }
                }
                item = child;  // expand
                break;
            }
        }
    }

    if (!candidate) return;

    // check if parent has no other children
    QTreeWidgetItem* parent           = candidate->parent();
    QTreeWidgetItem* parent_candidate = nullptr;
    while (parent != nullptr)
    {
        if (parent->childCount() > 1) break;
        parent_candidate = parent;
        parent           = parent->parent();
    }
    // delete highest parent that does not contain any further child
    if (parent_candidate)
        delete parent_candidate;
    else
        delete candidate;

    // check if our tree is empty, in that case erase group box as well
    if (group_box && group_it != _recent_signals.end() && tree->topLevelItemCount() == 0)
    {
        int task_id = std::get<2>(*group_it);
        group_box->deleteLater();
        _recent_signals.erase(group_it);
        emit taskRemoved(task_id);
    }
}

void SignalCollectionWidget::getGroupInfo(const QTreeWidget* tree, const QTreeWidgetItem* item, QString* name_out, int* id_out)
{
    // search if member of a "recent" group and extract task_id and group_box pointer
    int task_id   = _current_task_id;  // default construct with current_task_id
    auto group_it = std::find_if(_recent_signals.begin(), _recent_signals.end(),
                                 [tree](const std::tuple<CollapsableGroupBox*, QWidget*, int> item) { return tree == std::get<1>(item); });
    if (group_it != _recent_signals.end()) task_id = std::get<2>(*group_it);

    if (id_out) *id_out = task_id;

    // get namespace:
    if (name_out)
    {
        const QTreeWidgetItem* parent = item->parent();
        QString name                  = item->text(0);
        while (parent != nullptr)
        {
            name   = parent->text(0) + corbo::SIGNAL_NAMESPACE_DELIMITER + name;
            parent = parent->parent();
        }
        *name_out = name;
    }
}

void SignalCollectionWidget::removeSignal(QTreeWidget* tree, QTreeWidgetItem* item)
{
    // search if member of a "recent" group and extract task_id and group_box pointer
    int task_id   = _current_task_id;  // default construct with current_task_id
    auto group_it = std::find_if(_recent_signals.begin(), _recent_signals.end(),
                                 [tree](const std::tuple<CollapsableGroupBox*, QWidget*, int> item) { return tree == std::get<1>(item); });
    if (group_it != _recent_signals.end()) task_id = std::get<2>(*group_it);

    // first remove signals in the signal_helper
    SignalWidget* signal_widget = dynamic_cast<SignalWidget*>(tree->itemWidget(item, 0));
    if (signal_widget)
    {
        _signal_helper->removeSignal(signal_widget->key(), signal_widget->valueIdx());
    }
    else
    {
        // get namespace:
        QTreeWidgetItem* parent = item->parent();
        QString name            = item->text(0);
        while (parent != nullptr)
        {
            name   = parent->text(0) + corbo::SIGNAL_NAMESPACE_DELIMITER + name;
            parent = parent->parent();
        }
        _signal_helper->clearNamespace(name, task_id);
    }
}

void SignalCollectionWidget::recentGroupContextMenu(const QPoint& point, CollapsableGroupBox* group, int task_id)
{
    if (!group) return;

    QMenu menu(this);

    QAction* restore_param_action = new QAction(tr("Restore Parameters"), this);
    connect(restore_param_action, &QAction::triggered, [this, task_id]() {
        if (QMessageBox::Yes ==
            QMessageBox(QMessageBox::Information, "Restore Parameters", "Do you really want to restore parameters?",
                        QMessageBox::Yes | QMessageBox::No)
                .exec())
        {
            emit requestTaskParameterBackup(task_id);
        }
    });
    menu.addAction(restore_param_action);

    menu.addSeparator();

    // Export whole group
    QMenu* menu_export = menu.addMenu(tr("Export to ..."));

    auto it = _signal_helper->getSignalTreeRef().find(task_id);
    if (it != _signal_helper->getSignalTreeRef().end())
    {
        const CommonSignalTarget::SignalGroup* signal_group = &it->second.getSignals();

        // add all exporters that support 'SignalGroup export'
        for (const auto& item : DataExporterInterface::getFactory().getObjectMap())
        {
            DataExporterInterface::Ptr exporter_in_factory = item.second;
            if (exporter_in_factory->isSupportingSignalGroup())
            {
                QAction* export_group_action = new QAction(QString::fromStdString(exporter_in_factory->getFormatName()), this);
                connect(export_group_action, &QAction::triggered, [this, signal_group, exporter_in_factory]() {
                    DataExporterInterface::Ptr exporter =
                        exporter_in_factory->getInstance();  // create new exporter instance since we do not want to modify the one in the factory
                    QString format_name = QString::fromStdString(exporter->getFormatName());
                    QString file_suffix = QString::fromStdString(exporter->getFileSuffix());
                    std::string filename =
                        QFileDialog::getSaveFileName(this, "Save " + format_name, "data" + file_suffix, format_name + " File (* " + file_suffix + ")")
                            .toStdString();
                    exporter->exportSignalGroup(filename, (*signal_group));
                });
                menu_export->addAction(export_group_action);
            }
        }
    }

    menu.addSeparator();

    QAction* remove_action = new QAction(tr("Remove Signals"), this);
    connect(remove_action, &QAction::triggered, [this, task_id]() { _signal_helper->clearSeries(task_id); });
    menu.addAction(remove_action);

    menu.exec(group->mapToGlobal(point));
}

void SignalCollectionWidget::namespaceContextMenu(const QPoint& point, QTreeWidget* tree)
{
    if (!tree) return;

    // skip context menu for available signals, allow only for recent measurements!
    if (tree == _signal_tree) return;

    QTreeWidgetItem* item = tree->itemAt(point);

    if (!item) return;

    QMenu menu(this);

    SignalWidget* signal_widget = dynamic_cast<SignalWidget*>(tree->itemWidget(item, 0));
    if (signal_widget)
    {
        // this is a single signal (AND single value of a signal, e.g. a single component of the time-series object)
        addContextActionsSingleValueSignal(*signal_widget, menu);

        // add further actions
        QAction* remove_action = new QAction(tr("Remove Signal"), this);
        connect(remove_action, &QAction::triggered, [this, tree, item]() { removeSignal(tree, item); });
        menu.addAction(remove_action);
    }
    else
    {
        // get namespace:
        QString name;
        int id = 0;
        getGroupInfo(tree, item, &name, &id);

        // Export whole group
        QMenu* menu_export = menu.addMenu(tr("Export to ..."));

        auto it = _signal_helper->getSignalTreeRef().find(id);
        if (it != _signal_helper->getSignalTreeRef().end())
        {
            CommonSignalTarget::SignalGroup* signal_group = it->second.getGroup(name.toStdString());
            if (signal_group)
            {
                // add all exporters that support 'SignalGroup export'
                for (const auto& item : DataExporterInterface::getFactory().getObjectMap())
                {
                    DataExporterInterface::Ptr exporter_in_factory = item.second;
                    if (exporter_in_factory->isSupportingSignalGroup())
                    {
                        QAction* export_group_action = new QAction(QString::fromStdString(exporter_in_factory->getFormatName()), this);
                        connect(export_group_action, &QAction::triggered, [this, signal_group, exporter_in_factory]() {
                            DataExporterInterface::Ptr exporter =
                                exporter_in_factory
                                    ->getInstance();  // create new exporter instance since we do not want to modify the one in the factory
                            QString format_name  = QString::fromStdString(exporter->getFormatName());
                            QString file_suffix  = QString::fromStdString(exporter->getFileSuffix());
                            std::string filename = QFileDialog::getSaveFileName(this, "Save " + format_name, "data" + file_suffix,
                                                                                format_name + " File (* " + file_suffix + ")")
                                                       .toStdString();
                            exporter->exportSignalGroup(filename, (*signal_group));
                        });
                        menu_export->addAction(export_group_action);
                    }
                }
            }
            else
            {
                // we might have found a multi-valued signal which is not a group in the signal target
                SignalInterface::Ptr signal = it->second.getSignal(name.toStdString());
                if (signal)
                {
                    addContextActionsMultiValueSignal(signal, menu, *menu_export);
                }
                else
                {
                    PRINT_ERROR_NAMED("Cannot retrieve signal group " << name.toStdString());
                }
            }
        }

        menu.addSeparator();

        QAction* remove_action = new QAction(tr("Remove Signals"), this);
        connect(remove_action, &QAction::triggered, [this, tree, item]() { removeSignal(tree, item); });
        menu.addAction(remove_action);
    }

    menu.exec(tree->mapToGlobal(point));
}

void SignalCollectionWidget::addContextActionsSingleValueSignal(SignalWidget& signal_widget, QMenu& menu)
{
    // first check if we can access the signal data which we want to modify
    SignalHelper::SignalData* signal_data = _signal_helper->signalData(signal_widget.key());
    int signal_value_idx                  = signal_widget.valueIdx();
    if (!signal_data) return;

    // Test for TimeSeriesSignal
    TimeSeriesSignal::Ptr ts_signal = std::dynamic_pointer_cast<TimeSeriesSignal>(signal_data->signal);
    if (ts_signal)
    {
        addContextActionsTimeSeriesSignal(ts_signal, menu, signal_value_idx);
    }

    // Test for TimeSeriesSequenceSignal
    TimeSeriesSequenceSignal::Ptr ts_seq_signal = std::dynamic_pointer_cast<TimeSeriesSequenceSignal>(signal_data->signal);
    if (ts_seq_signal)
    {
        addContextActionsTimeSeriesSequenceSignal(ts_seq_signal, menu, signal_value_idx);
    }

    // Test for IndexedValuesSetSignal
    IndexedValuesSetSignal::Ptr ivs_signal = std::dynamic_pointer_cast<IndexedValuesSetSignal>(signal_data->signal);
    if (ivs_signal)
    {
        addContextActionsIndexedValuesSetSignal(ivs_signal, menu);
    }

    // Test for MatrixSignal
    MatrixSignal::Ptr matrix_signal = std::dynamic_pointer_cast<MatrixSignal>(signal_data->signal);
    if (matrix_signal)
    {
        addContextActionsMatrixSignal(matrix_signal, menu);
    }

    // Test for MatrixSetSignal
    MatrixSetSignal::Ptr matrix_set_signal = std::dynamic_pointer_cast<MatrixSetSignal>(signal_data->signal);
    if (matrix_set_signal)
    {
        addContextActionsMatrixSetSignal(matrix_set_signal, menu);
    }
}

void SignalCollectionWidget::addContextActionsMultiValueSignal(SignalInterface::Ptr signal, QMenu& menu, QMenu& menu_export)
{
    // Test for TimeSeriesSignal
    TimeSeriesSignal::Ptr ts_signal = std::dynamic_pointer_cast<TimeSeriesSignal>(signal);
    if (ts_signal)
    {
        addContextActionsTimeSeriesSignal(ts_signal, menu, &menu_export);
    }

    // Test for TimeSeriesSequenceSignal
    TimeSeriesSequenceSignal::Ptr ts_seq_signal = std::dynamic_pointer_cast<TimeSeriesSequenceSignal>(signal);
    if (ts_seq_signal)
    {
        addContextActionsTimeSeriesSequenceSignal(ts_seq_signal, menu, &menu_export);
    }

    // Test for IndexedValuesSetSignal
    IndexedValuesSetSignal::Ptr ivs_signal = std::dynamic_pointer_cast<IndexedValuesSetSignal>(signal);
    if (ivs_signal)
    {
        addContextActionsIndexedValuesSetSignal(ivs_signal, menu, &menu_export);
    }
}

void SignalCollectionWidget::addContextActionsTimeSeriesSignal(TimeSeriesSignal::Ptr ts_signal, QMenu& menu, int signal_value_idx, QMenu* menu_export)
{
    if (!ts_signal) return;

    // QMenu* menu_calc = menu.addMenu(tr("Calculations")); // add menu here. We might move that before the signal checks if we add other
    // calculations
    QMenu* menu_normalization = menu.addMenu(tr("Normalization"));

    // normalize by first value
    QAction* normalize_firstval_action = new QAction(tr("Normalize by First Value"), this);
    connect(normalize_firstval_action, &QAction::triggered,
            [ts_signal, signal_value_idx]() { ts_signal->getTimeSeriesRaw()->normalize(TimeSeries::Normalization::FirstValue, signal_value_idx); });
    menu_normalization->addAction(normalize_firstval_action);

    // normalize by absolute maximum
    QAction* normalize_absmax_action = new QAction(tr("Normalize by Absolute Maximum"), this);
    connect(normalize_absmax_action, &QAction::triggered, [ts_signal, signal_value_idx]() {
        ts_signal->getTimeSeriesRaw()->normalize(TimeSeries::Normalization::AbsoluteMaximum, signal_value_idx);
    });
    menu_normalization->addAction(normalize_absmax_action);

    // normalize by mean
    QAction* normalize_mean_action = new QAction(tr("Normalize by Mean"), this);
    connect(normalize_mean_action, &QAction::triggered,
            [ts_signal, signal_value_idx]() { ts_signal->getTimeSeriesRaw()->normalize(TimeSeries::Normalization::Mean, signal_value_idx); });
    menu_normalization->addAction(normalize_mean_action);

    // add separator after normalziation menu
    menu.addSeparator();

    QMenu* menu_export_aux = menu_export ? menu_export : menu.addMenu(tr("Export to ..."));

    // add all exporters that support 'TimeSeriesSignal export'
    for (const auto& item : DataExporterInterface::getFactory().getObjectMap())
    {
        DataExporterInterface::Ptr exporter_in_factory = item.second;
        if (exporter_in_factory->isSupportingTimeSeriesSignal())
        {
            QAction* export_group_action = new QAction(QString::fromStdString(exporter_in_factory->getFormatName()), this);
            connect(export_group_action, &QAction::triggered, [this, ts_signal, exporter_in_factory]() {
                DataExporterInterface::Ptr exporter =
                    exporter_in_factory->getInstance();  // create new exporter instance since we do not want to modify the one in the factory
                QString format_name = QString::fromStdString(exporter->getFormatName());
                QString file_suffix = QString::fromStdString(exporter->getFileSuffix());
                std::string filename =
                    QFileDialog::getSaveFileName(this, "Save " + format_name, "data" + file_suffix, format_name + " File (* " + file_suffix + ")")
                        .toStdString();
                exporter->exportTimeSeriesSignal(filename, (*ts_signal));
            });
            menu_export_aux->addAction(export_group_action);
        }
    }
}

void SignalCollectionWidget::addContextActionsTimeSeriesSignal(TimeSeriesSignal::Ptr ts_signal, QMenu& menu, QMenu* menu_export)
{
    if (!ts_signal) return;

    QMenu* menu_export_aux = menu_export ? menu_export : menu.addMenu(tr("Export to ..."));

    // add all exporters that support 'TimeSeriesSignal export'
    for (const auto& item : DataExporterInterface::getFactory().getObjectMap())
    {
        DataExporterInterface::Ptr exporter_in_factory = item.second;
        if (exporter_in_factory->isSupportingTimeSeriesSignal())
        {
            QAction* export_group_action = new QAction(QString::fromStdString(exporter_in_factory->getFormatName()), this);
            connect(export_group_action, &QAction::triggered, [this, ts_signal, exporter_in_factory]() {
                DataExporterInterface::Ptr exporter =
                    exporter_in_factory->getInstance();  // create new exporter instance since we do not want to modify the one in the factory
                QString format_name = QString::fromStdString(exporter->getFormatName());
                QString file_suffix = QString::fromStdString(exporter->getFileSuffix());
                std::string filename =
                    QFileDialog::getSaveFileName(this, "Save " + format_name, "data" + file_suffix, format_name + " File (* " + file_suffix + ")")
                        .toStdString();
                exporter->exportTimeSeriesSignal(filename, (*ts_signal));
            });
            menu_export_aux->addAction(export_group_action);
        }
    }
}

void SignalCollectionWidget::addContextActionsTimeSeriesSequenceSignal(TimeSeriesSequenceSignal::Ptr ts_seq_signal, QMenu& menu, int signal_value_idx,
                                                                       QMenu* menu_export)
{
    if (!ts_seq_signal) return;

    QMenu* menu_export_aux = menu_export ? menu_export : menu.addMenu(tr("Export to ..."));

    // add all exporters that support 'TimeSeriesSignal export'
    for (const auto& item : DataExporterInterface::getFactory().getObjectMap())
    {
        DataExporterInterface::Ptr exporter_in_factory = item.second;
        if (exporter_in_factory->isSupportingTimeSeriesSequenceSignal())
        {
            QAction* export_group_action = new QAction(QString::fromStdString(exporter_in_factory->getFormatName()), this);
            connect(export_group_action, &QAction::triggered, [this, ts_seq_signal, exporter_in_factory]() {
                DataExporterInterface::Ptr exporter =
                    exporter_in_factory->getInstance();  // create new exporter instance since we do not want to modify the one in the factory
                QString format_name = QString::fromStdString(exporter->getFormatName());
                QString file_suffix = QString::fromStdString(exporter->getFileSuffix());
                std::string filename =
                    QFileDialog::getSaveFileName(this, "Save " + format_name, "data" + file_suffix, format_name + " File (* " + file_suffix + ")")
                        .toStdString();
                exporter->exportTimeSeriesSequenceSignal(filename, (*ts_seq_signal));
            });
            menu_export_aux->addAction(export_group_action);
        }
    }
}

void SignalCollectionWidget::addContextActionsTimeSeriesSequenceSignal(TimeSeriesSequenceSignal::Ptr ts_seq_signal, QMenu& menu, QMenu* menu_export)
{
    if (!ts_seq_signal) return;

    QMenu* menu_export_aux = menu_export ? menu_export : menu.addMenu(tr("Export to ..."));

    // add all exporters that support 'TimeSeriesSignal export'
    for (const auto& item : DataExporterInterface::getFactory().getObjectMap())
    {
        DataExporterInterface::Ptr exporter_in_factory = item.second;
        if (exporter_in_factory->isSupportingTimeSeriesSequenceSignal())
        {
            QAction* export_group_action = new QAction(QString::fromStdString(exporter_in_factory->getFormatName()), this);
            connect(export_group_action, &QAction::triggered, [this, ts_seq_signal, exporter_in_factory]() {
                DataExporterInterface::Ptr exporter =
                    exporter_in_factory->getInstance();  // create new exporter instance since we do not want to modify the one in the factory
                QString format_name = QString::fromStdString(exporter->getFormatName());
                QString file_suffix = QString::fromStdString(exporter->getFileSuffix());
                std::string filename =
                    QFileDialog::getSaveFileName(this, "Save " + format_name, "data" + file_suffix, format_name + " File (* " + file_suffix + ")")
                        .toStdString();
                exporter->exportTimeSeriesSequenceSignal(filename, (*ts_seq_signal));
            });
            menu_export_aux->addAction(export_group_action);
        }
    }
}

void SignalCollectionWidget::addContextActionsIndexedValuesSetSignal(IndexedValuesSetSignal::Ptr ivs_signal, QMenu& menu, QMenu* menu_export)
{
    if (!ivs_signal) return;

    QMenu* menu_export_aux = menu_export ? menu_export : menu.addMenu(tr("Export to ..."));

    // add all exporters that support 'IndexedValuesSetSignal export'
    for (const auto& item : DataExporterInterface::getFactory().getObjectMap())
    {
        DataExporterInterface::Ptr exporter_in_factory = item.second;
        if (exporter_in_factory->isSupportingIndexedValuesSetSignal())
        {
            QAction* export_group_action = new QAction(QString::fromStdString(exporter_in_factory->getFormatName()), this);
            connect(export_group_action, &QAction::triggered, [this, ivs_signal, exporter_in_factory]() {
                DataExporterInterface::Ptr exporter =
                    exporter_in_factory->getInstance();  // create new exporter instance since we do not want to modify the one in the factory
                QString format_name = QString::fromStdString(exporter->getFormatName());
                QString file_suffix = QString::fromStdString(exporter->getFileSuffix());
                std::string filename =
                    QFileDialog::getSaveFileName(this, "Save " + format_name, "data" + file_suffix, format_name + " File (* " + file_suffix + ")")
                        .toStdString();
                exporter->exportIndexedValuesSetSignal(filename, (*ivs_signal));
            });
            menu_export_aux->addAction(export_group_action);
        }
    }
}

void SignalCollectionWidget::addContextActionsMatrixSignal(MatrixSignal::Ptr mat_signal, QMenu& menu, QMenu* menu_export)
{
    if (!mat_signal) return;

    QMenu* menu_export_aux = menu_export ? menu_export : menu.addMenu(tr("Export to ..."));

    // add all exporters that support 'TimeSeriesSignal export'
    for (const auto& item : DataExporterInterface::getFactory().getObjectMap())
    {
        DataExporterInterface::Ptr exporter_in_factory = item.second;
        if (exporter_in_factory->isSupportingMatrixSignal())
        {
            QAction* export_group_action = new QAction(QString::fromStdString(exporter_in_factory->getFormatName()), this);
            connect(export_group_action, &QAction::triggered, [this, mat_signal, exporter_in_factory]() {
                DataExporterInterface::Ptr exporter =
                    exporter_in_factory->getInstance();  // create new exporter instance since we do not want to modify the one in the factory
                QString format_name = QString::fromStdString(exporter->getFormatName());
                QString file_suffix = QString::fromStdString(exporter->getFileSuffix());
                std::string filename =
                    QFileDialog::getSaveFileName(this, "Save " + format_name, "data" + file_suffix, format_name + " File (* " + file_suffix + ")")
                        .toStdString();
                exporter->exportMatrixSignal(filename, (*mat_signal));
            });
            menu_export_aux->addAction(export_group_action);
        }
    }
}

void SignalCollectionWidget::addContextActionsMatrixSetSignal(MatrixSetSignal::Ptr mat_set_signal, QMenu& menu, QMenu* menu_export)
{
    if (!mat_set_signal) return;

    QMenu* menu_export_aux = menu_export ? menu_export : menu.addMenu(tr("Export to ..."));

    // add all exporters that support 'TimeSeriesSignal export'
    for (const auto& item : DataExporterInterface::getFactory().getObjectMap())
    {
        DataExporterInterface::Ptr exporter_in_factory = item.second;
        if (exporter_in_factory->isSupportingMatrixSetSignal())
        {
            QAction* export_group_action = new QAction(QString::fromStdString(exporter_in_factory->getFormatName()), this);
            connect(export_group_action, &QAction::triggered, [this, mat_set_signal, exporter_in_factory]() {
                DataExporterInterface::Ptr exporter =
                    exporter_in_factory->getInstance();  // create new exporter instance since we do not want to modify the one in the factory
                QString format_name = QString::fromStdString(exporter->getFormatName());
                QString file_suffix = QString::fromStdString(exporter->getFileSuffix());
                std::string filename =
                    QFileDialog::getSaveFileName(this, "Save " + format_name, "data" + file_suffix, format_name + " File (* " + file_suffix + ")")
                        .toStdString();
                exporter->exportMatrixSetSignal(filename, (*mat_set_signal));
            });
            menu_export_aux->addAction(export_group_action);
        }
    }
}

}  // namespace gui
}  // namespace corbo
