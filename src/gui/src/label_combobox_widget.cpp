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

#include <corbo-gui/label_combobox_widget.h>
#include <QHBoxLayout>

namespace corbo {
namespace gui {

LabelComboBoxWidget::LabelComboBoxWidget(const QString& label, QWidget* parent) : QWidget(parent)
{
    setSizePolicy(QSizePolicy::Policy::Preferred, QSizePolicy::Policy::Preferred);

    QHBoxLayout* box = new QHBoxLayout(this);
    box->setContentsMargins(0, 0, 0, 0);
    box->setAlignment(Qt::AlignLeft);
    _label = new QLabel(label);
    box->addWidget(_label);
    _combobox = new QComboBox;

    _label->setBuddy(_combobox);
    box->addWidget(_combobox);
}

QSize LabelComboBoxWidget::sizeHint() const { return QSize(200, 20); }

}  // namespace gui
}  // namespace corbo
