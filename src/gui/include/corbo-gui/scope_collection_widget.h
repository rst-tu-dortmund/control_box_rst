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

#ifndef SRC_GUI_INCLUDE_CORBO_GUI_SCOPE_COLLECTION_WIDGET_H_
#define SRC_GUI_INCLUDE_CORBO_GUI_SCOPE_COLLECTION_WIDGET_H_

#include <corbo-gui/signal_collection_widget.h>
#include <corbo-gui/signal_helper.h>
#include <corbo-gui/utilities.h>
#include <QTimer>
#include <QVBoxLayout>
#include <QWidget>

namespace corbo {
namespace gui {

class ScopeCollectionWidget : public QWidget
{
    Q_OBJECT

 public:
    explicit ScopeCollectionWidget(SignalHelper::ConstPtr signal_helper, QWidget* parent = 0);
    virtual ~ScopeCollectionWidget();

    QSize sizeHint() const override;

 public slots:
    // void addSignal(const QString& key, const SignalInterface& signal, SignalHelper::SignalData& signal_data, bool updated);
    void addMeasurement(const QString& key, Measurement::ConstPtr measurement, SignalHelper::SignalData& signal_data, bool first);
    void removeSignal(const QString& key, int value_idx);
    void initializeTask(int task_id);
    void resizeScopeAxes();
    void addScope();
    void closeAllScopes();

 signals:
    void scopeMeasurementUpdate(const QString& key, Measurement::ConstPtr measurement, SignalHelper::SignalData& signal_data, bool first);
    void scopeSignalRemoval(const QString& key, int value_idx);
    void scopeTaskInitialization(int task_id, bool inherit_signals);
    void requestScopeAxesResize();
    void previewTimeUpdate(double preview_time);

 protected:
    void createMenu();
    void createScopeArea();

 private:
    QVBoxLayout* _main_layout;
    QVBoxLayout* _scope_layout;

    SignalHelper::ConstPtr _signal_helper;

    bool _inherit_signals = true;

    // can be changed (e.g. by the dial) to change the current time for TimeSeriesSequence plots:
    double _current_preview_time     = 0.0;
    QTimer* _dial_preview_time_timer = nullptr;
};

}  // namespace gui
}  // namespace corbo

#endif  // SRC_GUI_INCLUDE_CORBO_GUI_SCOPE_COLLECTION_WIDGET_H_
