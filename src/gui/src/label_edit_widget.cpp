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

#include <corbo-gui/label_edit_widget.h>
#include <QHBoxLayout>

namespace corbo {
namespace gui {

LabelEditWidget::LabelEditWidget(const QString& label, const QString& default_text, QWidget* parent) : QWidget(parent)
{
    setSizePolicy(QSizePolicy::Policy::Preferred, QSizePolicy::Policy::Preferred);

    QHBoxLayout* box = new QHBoxLayout(this);
    box->setContentsMargins(5, 0, 0, 0);
    box->setAlignment(Qt::AlignLeft);
    _label = new QLabel(label);
    box->addWidget(_label);
    _edit = new QLineEdit;
    _edit->setText(default_text);

    _label->setBuddy(_edit);
    box->addWidget(_edit);
}

QSize LabelEditWidget::sizeHint() const { return QSize(200, 20); }

}  // namespace gui
}  // namespace corbo
