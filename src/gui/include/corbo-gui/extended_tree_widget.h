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

#ifndef SRC_GUI_INCLUDE_CORBO_GUI_EXTENDED_TREE_WIDGET_H_
#define SRC_GUI_INCLUDE_CORBO_GUI_EXTENDED_TREE_WIDGET_H_

#include <QTreeWidget>

namespace corbo {
namespace gui {

// the common tree widget does not support deselection of items when clicking somewhere else in the widget
// Hence, we add this functionality:
class ExtendedTreeWidget : public QTreeWidget
{
    Q_OBJECT

 public:
    explicit ExtendedTreeWidget(QWidget* parent = nullptr);
    virtual ~ExtendedTreeWidget() {}

    QSize sizeHint() const override;

 public slots:
    void updateSizeHint() { updateGeometry(); }

 protected:
    int recursiveHeightHint(QTreeWidgetItem* item) const;

 private:
    void mousePressEvent(QMouseEvent* event) override;
};

}  // namespace gui
}  // namespace corbo

#endif  // SRC_GUI_INCLUDE_CORBO_GUI_EXTENDED_TREE_WIDGET_H_
