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

#ifndef SRC_GUI_INCLUDE_CORBO_GUI_LABEL_COMBOBOX_WIDGET_H_
#define SRC_GUI_INCLUDE_CORBO_GUI_LABEL_COMBOBOX_WIDGET_H_

#include <QComboBox>
#include <QLabel>

namespace corbo {
namespace gui {

class LabelComboBoxWidget : public QWidget
{
    Q_OBJECT

 public:
    explicit LabelComboBoxWidget(QWidget* parent = 0) : LabelComboBoxWidget("", parent) {}
    explicit LabelComboBoxWidget(const QString& label, QWidget* parent = 0);
    ~LabelComboBoxWidget() {}

    QSize sizeHint() const override;

    void setLabel(const QString& label) { _label->setText(label); }

    QComboBox* widgetComboBox() { return _combobox; }
    const QComboBox* widgetComboBox() const { return _combobox; }

    QLabel* widgetLabel() { return _label; }
    const QLabel* widgetLabel() const { return _label; }

 protected:
 private:
    QComboBox* _combobox;
    QLabel* _label;
};

}  // namespace gui
}  // namespace corbo

#endif  // SRC_GUI_INCLUDE_CORBO_GUI_LABEL_COMBOBOX_WIDGET_H_
