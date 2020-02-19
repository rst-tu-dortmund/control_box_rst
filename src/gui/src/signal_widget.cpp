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

#include <corbo-gui/signal_widget.h>

#include <QApplication>
#include <QDrag>
#include <QLabel>
#include <QMimeData>
#include <QMouseEvent>
#include <QStyleHints>

#include <corbo-core/console.h>

namespace corbo {
namespace gui {

SignalWidget::SignalWidget(const QString& title, const QString& key, int value_idx, QWidget* parent)
    : QWidget(parent), _key(key), _value_idx(value_idx)
{
    setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);

    _layout = new QHBoxLayout(this);
    _layout->setContentsMargins(0, 0, 0, 0);
    _layout->setAlignment(Qt::AlignTop | Qt::AlignLeft);

    QLabel* label = new QLabel(title);
    _layout->addWidget(label);
}

SignalWidget::~SignalWidget() {}

QSize SignalWidget::sizeHint() const { return QSize(200, 15); }

void SignalWidget::mousePressEvent(QMouseEvent* event)
{
    if (event->button() == Qt::LeftButton) _drag_start_position = event->pos();

    QWidget::mousePressEvent(event);
}

void SignalWidget::mouseMoveEvent(QMouseEvent* event)
{
    if (!(event->buttons() & Qt::LeftButton)) return;
    if ((event->pos() - _drag_start_position).manhattanLength() < QApplication::startDragDistance()) return;

    // drag and drop action detected:
    QDrag* drag          = new QDrag(this);
    QMimeData* mime_data = new QMimeData;
    mime_data->setText(_key + "::" + QString::number(_value_idx));
    drag->setMimeData(mime_data);

    Qt::DropAction drop_action = drag->exec();
}

}  // namespace gui
}  // namespace corbo
