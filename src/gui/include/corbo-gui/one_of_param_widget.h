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

#ifndef SRC_GUI_INCLUDE_CORBO_GUI_ONE_OF_PARAM_WIDGET_H_
#define SRC_GUI_INCLUDE_CORBO_GUI_ONE_OF_PARAM_WIDGET_H_

#include <corbo-gui/label_combobox_widget.h>
#include <corbo-gui/parameter_cache.h>

#include <corbo-communication/message_parser.h>

#include <QHash>
#include <QPair>
#include <QSet>
#include <QString>
#include <QVBoxLayout>

#include <list>
#include <string>

namespace corbo {
namespace gui {

class OneOfParamWidget : public QWidget
{
    Q_OBJECT

 public:
    explicit OneOfParamWidget(ParameterCache* cache = nullptr, QWidget* parent = nullptr) : OneOfParamWidget("", cache, parent) {}
    explicit OneOfParamWidget(const QString& label, ParameterCache* cache = nullptr, QWidget* parent = nullptr);
    virtual ~OneOfParamWidget();

    QSize sizeHint() const override;

    const std::list<std::string>& nestedFieldName() const { return _nested_field_name; }
    std::list<std::string>& nestedFieldName() { return _nested_field_name; }
    std::string nestedFieldNameUnrolled(const std::string& delimiter);
    
    void addDescription(const QString& description);

 public slots:

    void registerParameter(const QString& label, const MessageParser::FieldInformation& info);
    void selectParameter(const QString& label);

 signals:
    void currentParameterChanged(const QString& label);
    void signalUpdateRequested();

 protected:
    void addParametersToCache();
    bool restoreParametersFromCache();

 private:
    QVBoxLayout* _layout;
    QWidget* _param_widget = nullptr;
    LabelComboBoxWidget* _combobox;

    std::list<std::string> _nested_field_name;

    google::protobuf::Message* _message = nullptr;
    QHash<QString, const google::protobuf::FieldDescriptor*> _fields;

    QString _selected_item;

    bool _selected_item_parsed = false;  //!< Keep track if we have already parsed the selected item to avoid expanding a wrong protobuf message

    ParameterCache* _param_cache;
};

}  // namespace gui
}  // namespace corbo

#endif  // SRC_GUI_INCLUDE_CORBO_GUI_ONE_OF_PARAM_WIDGET_H_
