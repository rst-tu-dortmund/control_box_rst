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

#ifndef SRC_GUI_INCLUDE_CORBO_GUI_PARAMETER_WIDGET_H_
#define SRC_GUI_INCLUDE_CORBO_GUI_PARAMETER_WIDGET_H_

#include <QHash>
#include <QPair>
#include <QStack>
#include <QVBoxLayout>
#include <QWidget>

#include <corbo-communication/message_parser.h>
#include <corbo-gui/collapsable_groupbox.h>
#include <corbo-gui/one_of_param_widget.h>
#include <corbo-gui/parameter_cache.h>
#include <QGroupBox>

#include <memory>
#include <string>
#include <tuple>
#include <vector>

namespace corbo {
namespace gui {

class ParameterWidget : public QWidget
{
    Q_OBJECT

 public:
    explicit ParameterWidget(ParameterCache* cache = nullptr, QWidget* parent = nullptr);
    virtual ~ParameterWidget();

    QSize sizeHint() const override;

 public slots:

    void generateFromMessage(std::shared_ptr<google::protobuf::Message> message);

    void generateFromAllocatedField(google::protobuf::Message* message, const google::protobuf::FieldDescriptor* field);

    const google::protobuf::Message* getMessage() const { return _param_message.get(); }
    google::protobuf::Message* getMessage() { return _param_message.get(); }

    void addParameterInt32(int value, const MessageParser::FieldInformation& info);
    void addParameterInt32Array(const std::vector<int>& values, const MessageParser::FieldInformation& info);
    void addParameterDouble(double value, const MessageParser::FieldInformation& info);
    void addParameterDoubleArray(const std::vector<double>& values, const MessageParser::FieldInformation& info);
    void addParameterBool(bool value, const MessageParser::FieldInformation& info);
    void addParameterBoolArray(const std::vector<bool>& values, const MessageParser::FieldInformation& info);
    void addParameterEnum(const std::string& value, const MessageParser::FieldInformation& info);
    void addParameterString(std::string value, const MessageParser::FieldInformation& info);
    void addInfoText(const std::string& text, const MessageParser::FieldInformation& info);

    void startGroup(const MessageParser::FieldInformation& info);
    void endGroup();

    bool hasParameters() const;

    void clearElements();

    const std::list<std::string>& nestedParentFieldNames() const { return _nested_parent_fields; }
    std::list<std::string>& nestedParentFieldNames() { return _nested_parent_fields; }
    QString parentFieldNames() const { return QString::fromStdString(MessageParser::nestedNamesToString(_nested_parent_fields)); }

 signals:
    void updatedOneOfField(const QString& text);
    void parameterInt32Updated(const QString& parameter, int value);
    void signalUpdateRequested();

 protected:
    void generateElements();
    void parse(google::protobuf::Message* message, const google::protobuf::FieldDescriptor* field = nullptr);

 private:
    void addSubWidget(QWidget* widget);
    void addOneOfField(const MessageParser::FieldInformation& info);

    std::shared_ptr<google::protobuf::Message> _param_message;

    QVBoxLayout* _layout;
    QStack<std::tuple<CollapsableGroupBox*, QLayout*, bool>> _groups;  //!< groupbox ptr, layout ptr, new msg flag (e.g. for loading default params)
    QHash<QString, OneOfParamWidget*> _oneof_widgets;

    ParameterCache* _param_cache;

    std::list<std::string> _nested_parent_fields;

    bool _has_parameters = false;
};

}  // namespace gui
}  // namespace corbo

#endif  // SRC_GUI_INCLUDE_CORBO_GUI_PARAMETER_WIDGET_H_
