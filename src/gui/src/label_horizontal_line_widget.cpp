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

#include <corbo-gui/label_horizontal_line_widget.h>
#include <QHBoxLayout>

namespace corbo {
namespace gui {

LabelHLineWidget::LabelHLineWidget(const QString& label, bool hline, QWidget* parent) : QWidget(parent)
{
    QHBoxLayout* box = new QHBoxLayout(this);
    box->setContentsMargins(0, 0, 0, 0);
    box->setAlignment(Qt::AlignLeft);
    _label = new QLabel(label);
    box->addWidget(_label);
    if (hline)
    {
        QFrame* frame = new QFrame;
        // frame->setMinimumWidth(50);
        frame->setFrameShape(QFrame::HLine);
        frame->setFrameShadow(QFrame::Sunken);
        box->addWidget(frame, 0, Qt::AlignLeft | Qt::AlignVCenter);
    }
}

QSize LabelHLineWidget::sizeHint() const { return QSize(200, 20); }

}  // namespace gui
}  // namespace corbo
