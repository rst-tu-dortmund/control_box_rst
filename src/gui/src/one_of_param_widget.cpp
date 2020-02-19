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

#include <corbo-communication/message_parser.h>
#include <corbo-gui/one_of_param_widget.h>
#include <corbo-gui/parameter_widget.h>

#include <corbo-core/console.h>
#include <QVBoxLayout>

#include <memory>
#include <string>

namespace corbo {
namespace gui {

OneOfParamWidget::OneOfParamWidget(const QString& label, ParameterCache* cache, QWidget* parent) : QWidget(parent), _param_cache(cache)
{
    setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

    _layout = new QVBoxLayout(this);
    _layout->setContentsMargins(5, 5, 0, 0);
    _layout->setAlignment(Qt::AlignTop | Qt::AlignLeft);

    _combobox = new LabelComboBoxWidget(label);

    _layout->addWidget(_combobox);

    connect(_combobox->widgetComboBox(), SIGNAL(currentTextChanged(const QString&)), this, SLOT(selectParameter(const QString&)));
}

OneOfParamWidget::~OneOfParamWidget() {}

QSize OneOfParamWidget::sizeHint() const { return QSize(200, 50); }

void OneOfParamWidget::registerParameter(const QString& label, const MessageParser::FieldInformation& info)
{
    _message = info.message_raw;
    _fields.insert(label, info.field_raw);

    // check if the current parameter is active
    // we also check if no one-of field is selected, then this field becomes active
    bool selected_field = (info.field_raw == info.oneof_selected_field) || !info.oneof_selected_field;

    if (selected_field) _selected_item_parsed = true;

    _combobox->widgetComboBox()->addItem(label);
    // this triggers selectParameter() if it is the first Item
    // (-> we might expand and allocate a wrong protobuf message, hence we
    // are keeping track of the selected value with member _selected_item_parsed.

    // if this parameter is active switch combobox item
    if (selected_field)
    {
        _combobox->widgetComboBox()->setCurrentText(label);  // should trigger selectParameter() callback
    }
}

void OneOfParamWidget::selectParameter(const QString& label)
{
    // see descirption in registerParameter()
    if (!_selected_item_parsed) return;

    // store old paramters (_selected_item)
    addParametersToCache();
    _selected_item = label;

    if (_param_widget)
    {
        _layout->removeWidget(_param_widget);
        delete _param_widget;
        _param_widget = nullptr;
    }

    restoreParametersFromCache();

    auto field = _fields.find(label);
    if (field != _fields.end())
    {
        ParameterWidget* params = new ParameterWidget(_param_cache);
        // forward nested names in case they are required for reconstructing the full parameter path/name
        params->nestedParentFieldNames() = _nested_field_name;
        params->generateFromAllocatedField(_message, field.value());
        if (params->hasParameters())
        {
            _layout->addWidget(params);
            _param_widget = params;
            params->connect(params, &ParameterWidget::signalUpdateRequested, [this]() { emit signalUpdateRequested(); });
            params->connect(params, &ParameterWidget::updatedOneOfField, [this](const QString& name) { emit currentParameterChanged(name); });
        }
        else
        {
            delete params;
            _param_widget = nullptr;
        }

        emit currentParameterChanged(label);
    }
}

void OneOfParamWidget::addParametersToCache()
{
    if (!_param_cache || _selected_item.isEmpty()) return;

    // store child oneof-fields as well
    QList<OneOfParamWidget*> widgets = _param_widget->findChildren<OneOfParamWidget*>("", Qt::FindDirectChildrenOnly);
    for (auto widget : widgets) widget->addParametersToCache();

    // store parameter of this one-if
    std::string param_path = ParameterCache::nestedNameListToKey(_nested_field_name);
    // add selected parameter
    param_path += "/" + _selected_item.toStdString();
    // store complete message, but we later only use the corresonding one-of-field
    // TODO(roesmann) here is place for improving efficency and required memory size
    _param_cache->toCache(param_path, *_message);
}

bool OneOfParamWidget::restoreParametersFromCache()
{
    std::string param_path = ParameterCache::nestedNameListToKey(_nested_field_name);
    // add selected parameter
    param_path += "/" + _selected_item.toStdString();

    std::unique_ptr<google::protobuf::Message> cache = _param_cache->fromCache(param_path);
    if (!cache) return false;

    const google::protobuf::FieldDescriptor* current_field = _fields[_selected_item];

    // only copy the required one-of-field
    _message->GetReflection()->SwapFields(_message, cache.get(), {current_field});
    // _message->CopyFrom(*cache);
    return true;
}

std::string OneOfParamWidget::nestedFieldNameUnrolled(const std::string& delimiter)
{
    auto it            = _nested_field_name.begin();
    std::string nnames = *it;
    for (std::advance(it, 1); it != _nested_field_name.end(); ++it) nnames += delimiter + *it;
    return nnames;
}
    
void OneOfParamWidget::addDescription(const QString& description)
{
    _layout->addWidget(new QLabel(description));
}

}  // namespace gui
}  // namespace corbo
