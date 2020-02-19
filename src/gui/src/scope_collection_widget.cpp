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

#include <corbo-gui/scope_collection_widget.h>

#include <corbo-core/console.h>
#include <corbo-gui/label_edit_widget.h>
#include <corbo-gui/scope_widget.h>
#include <corbo-gui/slider_center_stuck_widget.h>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QPushButton>
#include <QScrollArea>

namespace corbo {
namespace gui {

ScopeCollectionWidget::ScopeCollectionWidget(SignalHelper::ConstPtr signal_helper, QWidget* parent) : QWidget(parent), _signal_helper(signal_helper)
{
    setSizePolicy(QSizePolicy::Expanding, QSizePolicy::MinimumExpanding);

    _main_layout = new QVBoxLayout(this);
    // _layout->setContentsMargins(0, 0, 0, 0);
    _main_layout->setAlignment(Qt::AlignTop);

    createMenu();
    createScopeArea();
}

ScopeCollectionWidget::~ScopeCollectionWidget() {}

QSize ScopeCollectionWidget::sizeHint() const { return QSize(1500, 800); }

void ScopeCollectionWidget::createMenu()
{
    QHBoxLayout* menu_layout = new QHBoxLayout;

    // create button for adding new scopes
    QPushButton* new_scope_btn = new QPushButton(tr("Add Scope"));
    new_scope_btn->setMaximumWidth(100);

    QFont font = new_scope_btn->font();
    font.setPointSize(10);
    new_scope_btn->setFont(font);
    connect(new_scope_btn, &QPushButton::clicked, [this](bool) { addScope(); });

    menu_layout->addWidget(new_scope_btn);

    menu_layout->addStretch();

    // add text field to select the current preview time
    LabelEditWidget* tfs_edit = new LabelEditWidget(tr("Preview time"), "0");
    tfs_edit->setMaximumWidth(130);
    auto tfs_edit_changed = [tfs_edit, this]() {
        // get values from line edit
        bool ok;
        double value = tfs_edit->widgetLineEdit()->text().toDouble(&ok);
        if (!ok || value < 0)
        {
            value = 0;
            tfs_edit->setText("0");
        }
        _current_preview_time = value;
        emit previewTimeUpdate(_current_preview_time);
    };
    connect(tfs_edit->widgetLineEdit(), &QLineEdit::editingFinished, tfs_edit_changed);
    menu_layout->addWidget(tfs_edit);

    // add slider widget to provide the user with an oppertunity to scroll through time (e.g. for plotting TimeSeriesSequences)
    SliderCenterStuck* slider = new SliderCenterStuck(Qt::Horizontal);
    // slider->setMaximumHeight(10);
    slider->setMaximumWidth(80);
    slider->setToCenter();
    slider->setStuckToCenter(10);
    connect(slider, &SliderCenterStuck::valueChanged, this, [slider, tfs_edit, this](int value) {
        int speed = value - slider->getCenter();

        if (_dial_preview_time_timer)
        {
            _dial_preview_time_timer->stop();
            delete _dial_preview_time_timer;
            _dial_preview_time_timer = nullptr;
        }
        if (speed != 0)
        {
            _dial_preview_time_timer = new QTimer(this);

            double dt = (double)speed * (double)std::abs(speed) * 0.0001;  // create quadratic scaling behavior
            dt        = std::trunc(dt * 100) / 100;                        // only two decimal places
            connect(_dial_preview_time_timer, &QTimer::timeout, tfs_edit, [tfs_edit, dt, this]() {
                bool ok;
                double time   = tfs_edit->widgetLineEdit()->text().toDouble(&ok);
                if (!ok) time = 0;
                time += dt;
                if (time < 0)
                {
                    time = 0;
                    _dial_preview_time_timer->stop();
                }
                tfs_edit->setText(QString::number(time));
                _current_preview_time = time;
                emit previewTimeUpdate(_current_preview_time);
            });
            _dial_preview_time_timer->start(250);
        }
    });
    menu_layout->addWidget(slider);

    menu_layout->addSpacing(20);

    // checkbox for remembering signal association
    QCheckBox* inherit_signals = new QCheckBox(tr("Inherit Signals"));
    inherit_signals->setChecked(true);
    inherit_signals->setToolTip(tr("If selected, the signal-to-scope association is remembered for future tasks."));
    connect(inherit_signals, &QCheckBox::toggled, [this](bool checked) { _inherit_signals = checked; });

    menu_layout->addWidget(inherit_signals);

    _main_layout->addLayout(menu_layout);

    QFrame* hline = new QFrame;
    hline->setFrameShape(QFrame::HLine);
    hline->setFrameShadow(QFrame::Sunken);
    hline->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Maximum);
    _main_layout->addWidget(hline);
}

void ScopeCollectionWidget::createScopeArea()
{
    _scope_layout = new QVBoxLayout;
    _scope_layout->setAlignment(Qt::AlignTop | Qt::AlignLeft);
    _scope_layout->setContentsMargins(0, 0, 0, 0);

    // add scroll area for the content
    QScrollArea* scroll_area = new QScrollArea;
    // we need to create a group since the scrollarea only works with a singla widget rather than a whole layout
    QWidget* group = new QWidget;
    group->setLayout(_scope_layout);
    // scroll_area->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::MinimumExpanding);
    scroll_area->setWidgetResizable(true);
    scroll_area->setWidget(group);
    _main_layout->addWidget(scroll_area);
}

void ScopeCollectionWidget::addScope()
{
    if (!_signal_helper) return;

    ScopeWidget* scope_widget = new ScopeWidget(_signal_helper);
    scope_widget->setPreviewTime(_current_preview_time);

    // connect scope signal update
    connect(this, &ScopeCollectionWidget::scopeMeasurementUpdate, scope_widget,
            [scope_widget](const QString& key, Measurement::ConstPtr measurement, SignalHelper::SignalData& data, bool /*first*/) {
                scope_widget->addMeasurement(key, measurement, data);
                // else
                //    scope_widget->addSignal(key, signal);
                // TODO(roesmann) we have no value_idx here, since we currently
                //    only support drag and drop insertion
            });
    // conect scope signal removal
    connect(this, &ScopeCollectionWidget::scopeSignalRemoval, scope_widget, &ScopeWidget::removeSignal);
    // connect scope task update
    connect(this, &ScopeCollectionWidget::scopeTaskInitialization, scope_widget, &ScopeWidget::initializeTask);
    // connect scope axes resize
    connect(this, &ScopeCollectionWidget::requestScopeAxesResize, scope_widget, &ScopeWidget::rescaleAxes);
    // connect preview time / time from start
    connect(this, &ScopeCollectionWidget::previewTimeUpdate, scope_widget, &ScopeWidget::setPreviewTime);

    _scope_layout->addWidget(scope_widget);
}

void ScopeCollectionWidget::closeAllScopes()
{
    //    for (auto widget : _scope_layout->findChildren<QWidget*>(QString(), Qt::FindDirectChildrenOnly))
    //    {
    //        // This does not seem to work right now
    //        widget->deleteLater();
    //    }
    QLayoutItem* item;
    while ((item = _scope_layout->takeAt(0)) != nullptr)
    {
        delete item->widget();
        delete item;
    }
}

void ScopeCollectionWidget::addMeasurement(const QString& key, Measurement::ConstPtr measurement, SignalHelper::SignalData& signal_data, bool first)
{
    emit scopeMeasurementUpdate(key, measurement, signal_data, first);
}

void ScopeCollectionWidget::removeSignal(const QString& key, int value_idx) { emit scopeSignalRemoval(key, value_idx); }

void ScopeCollectionWidget::initializeTask(int task_id) { emit scopeTaskInitialization(task_id, _inherit_signals); }

void ScopeCollectionWidget::resizeScopeAxes() { emit requestScopeAxesResize(); }

}  // namespace gui
}  // namespace corbo
