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

#include <corbo-gui/extended_tree_widget.h>
#include <QHeaderView>

#include <algorithm>

namespace corbo {
namespace gui {

ExtendedTreeWidget::ExtendedTreeWidget(QWidget* parent) : QTreeWidget(parent)
{
    connect(this, SIGNAL(itemExpanded(QTreeWidgetItem*)), this, SLOT(updateSizeHint()));
    connect(this, SIGNAL(itemCollapsed(QTreeWidgetItem*)), this, SLOT(updateSizeHint()));

    // it is convenient for the user to also expand the child tree if only a single child is present:
    connect(this, &ExtendedTreeWidget::itemExpanded, [](QTreeWidgetItem* item) {
        if (item->childCount() == 1) item->child(0)->setExpanded(true);
    });
}

QSize ExtendedTreeWidget::sizeHint() const
{
    QSize size = QTreeWidget::sizeHint();
    int height = 2 * frameWidth();
    if (!isHeaderHidden() && header()) height += header()->sizeHint().height();

    for (int i = 0; i < topLevelItemCount(); ++i) height += recursiveHeightHint(topLevelItem(i));

    size.setHeight(height);
    return size;
}

int ExtendedTreeWidget::recursiveHeightHint(QTreeWidgetItem* item) const
{
    int item_height        = 0;
    int item_widget_height = 0;
    // get highest column
    for (int col_idx = 0; col_idx < item->columnCount(); ++col_idx)
    {
        // take font size into account for plain text items
        // hack: we add 5 for the rest of the frame
        // FIXME(roesmann)
        // int item_height_aux = std::max(item->sizeHint(col_idx).height(), item->font(col_idx).pointSize() + 5);
        int item_height_aux = visualItemRect(item).height();

        item_height = std::max(item_height, item_height_aux);
    }
    QWidget* item_widget                = itemWidget(item, 0);  // note, itemWidget does not imply a positive item->columnCount() !
    if (item_widget) item_widget_height = std::max(item_widget_height, item_widget->geometry().height());
    // TODO(roesmann) or item_widget->sizeHint().height()

    int height = std::max(item_height, item_widget_height);

    if (item->isExpanded())
    {
        int num_children = item->childCount();
        for (int child_idx = 0; child_idx < num_children; ++child_idx) height += recursiveHeightHint(item->child(child_idx));
    }
    return height;
}

void ExtendedTreeWidget::mousePressEvent(QMouseEvent* event)
{
    clearSelection();
    // QModelIndex item = indexAt(event->pos());
    // bool selected    = selectionModel()->isSelected(item);
    QTreeWidget::mousePressEvent(event);
    // if (selected) selectionModel()->select(item, QItemSelectionModel::Deselect);
}

}  // namespace gui
}  // namespace corbo
