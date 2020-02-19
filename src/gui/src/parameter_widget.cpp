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

#include <corbo-gui/parameter_widget.h>

#include <corbo-communication/message_parser.h>
#include <corbo-core/console.h>
#include <corbo-core/utilities.h>
#include <corbo-gui/collapsable_groupbox.h>
#include <corbo-gui/horizontal_button_group.h>
#include <corbo-gui/label_combobox_widget.h>
#include <corbo-gui/label_edit_widget.h>
#include <corbo-gui/label_slider_widget.h>
#include <corbo-gui/one_of_param_widget.h>
#include <corbo-gui/utilities.h>

#include <QCheckBox>
#include <QLabel>  // For testing only
#include <QMessageBox>
#include <QStringList>
#include <QVBoxLayout>

#include <google/protobuf/reflection.h>

#include <functional>
#include <memory>

namespace corbo {
namespace gui {

ParameterWidget::ParameterWidget(ParameterCache* cache, QWidget* parent) : QWidget(parent), _param_cache(cache)
{
    setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Maximum);

    _layout = new QVBoxLayout(this);
    _layout->setContentsMargins(5, 5, 5, 5);
    _layout->setAlignment(Qt::AlignTop | Qt::AlignLeft);
}

ParameterWidget::~ParameterWidget() {}

QSize ParameterWidget::sizeHint() const { return QSize(300, 50); }

void ParameterWidget::addOneOfField(const MessageParser::FieldInformation& info)
{
    if (!info.oneof_descriptor) return;  // no one of element;

    _has_parameters = true;  // even if we do not have real parameters, one-ofs should be displayed (they can also be selected)

    // create oneof widget if not yet created for a previous parameter
    OneOfParamWidget* oneof_widget;

    // QString oneof_widget_key = QString::fromStdString(info.oneof_descriptor->full_name());
    QString oneof_widget_key;
    auto it          = info.nested_names.begin();
    oneof_widget_key = QString::fromStdString(*it);
    std::advance(it, 1);
    for (; it != std::prev(info.nested_names.end()); ++it) oneof_widget_key += "." + QString::fromStdString(*it);

    // get one-of field name in parent message if available in order to allow multiple one-of of the
    // same type in a single message
    // QString one_of_field_name;
    // auto nested_name_it = info.nested_names.end();
    // take the entry before the last, since the last is already a field of the given one-of type.
    // if (info.nested_names.size() > 1)
    //{
    //    std::advance(nested_name_it, -2);
    //    one_of_field_name = QString::fromStdString(*nested_name_it);
    //}
    // oneof_widget_key += "." + one_of_field_name;  // for top-level-one-of, one_of_field_name is empty.

    // check if key is already present
    auto it_oneof = _oneof_widgets.find(oneof_widget_key);
    if (it_oneof != _oneof_widgets.end())
    {
        // widget found
        oneof_widget = it_oneof.value();
    }
    else
    {
        // create new oneof widget
        QString widget_label            = QString::fromStdString(info.oneof_descriptor->name());  // TODO(roesmann): separate label in protobuf
        oneof_widget                    = new OneOfParamWidget(widget_label, _param_cache);
        oneof_widget->nestedFieldName() = info.nested_names;

        oneof_widget->connect(oneof_widget, &OneOfParamWidget::currentParameterChanged,
                              [this](const QString& label) { emit updatedOneOfField(label); });
        oneof_widget->connect(oneof_widget, &OneOfParamWidget::signalUpdateRequested, [this]() { emit signalUpdateRequested(); });

        _oneof_widgets.insert(oneof_widget_key, oneof_widget);

        addSubWidget(oneof_widget);
    }
    // now register current field
    QString field_label = QString::fromStdString(info.field_raw->name());  // TODO(roesmann) actual label
    oneof_widget->registerParameter(field_label, info);
    return;
}

void ParameterWidget::addParameterInt32(int value, const MessageParser::FieldInformation& info)
{
    _has_parameters = true;

    messages::GuiType type = info.gui_type.first ? info.gui_type.second : messages::GuiType::TEXTEDIT;
    QString label          = info.label.first ? QString::fromStdString(info.label.second) : QString::fromStdString(info.field_name);
    QString description    = info.description.first ? QString::fromStdString(info.description.second) : "";
    int min                = info.min_value.first ? static_cast<int>(info.min_value.second) : CORBO_MIN_INT;
    int max                = info.max_value.first ? static_cast<int>(info.max_value.second) : CORBO_MAX_INT;

    // add field description label if specified
    if (info.print_description) addSubWidget(new QLabel(description));

    // check if we need to load a default value
    if (info.default_value.first && !_groups.empty() && std::get<2>(_groups.back()))
    {
        // overwrite value
        bool ok = false;
        value   = util::qstring_to_int(QString::fromStdString(info.default_value.second), &ok);
        if (!ok)
        {
            PRINT_ERROR_NAMED("Cannot convert default parameter: " << info.default_value.second << " to int.");
            value = 0;
        }
        info.message_raw->GetReflection()->SetInt32(info.message_raw, info.field_raw, value);
    }

    switch (type)
    {
        case messages::GuiType::SLIDER:
        {
            LabelSliderWidget* param_widget = new LabelSliderWidget(label, value, 1, 0);
            param_widget->setToolTip(description);
            // if min and max are not specified, set some reasonable default limits:
            if (!info.min_value.first) min = 0;
            if (!info.max_value.first) max = 100;
            param_widget->setMinMax(min, max);
            addSubWidget(param_widget);

            // set up modifer signal
            // message and field pointers must be valid for the lifetime of this object!
            auto fun_write_param = [info, this](double value) {
                int value_int = static_cast<int>(value);
                info.message_raw->GetReflection()->SetInt32(info.message_raw, info.field_raw, value_int);
                emit parameterInt32Updated(QString::fromStdString(MessageParser::nestedNamesToString(info.nested_names)), value_int);
            };
            connect(param_widget, &LabelSliderWidget::valueChanged, fun_write_param);
            break;
        }
        case messages::GuiType::TEXTEDIT:
        default:
        {
            LabelEditWidget* param_widget = new LabelEditWidget(label, util::int_to_qstring(value));
            param_widget->setToolTip(description);
            QRegExp regexp_no_decimal("^((?!\\.).)*$");
            param_widget->widgetLineEdit()->setValidator(new QRegExpValidator(regexp_no_decimal, param_widget->widgetLineEdit()));
            addSubWidget(param_widget);

            // set up modifer signal
            // message and field pointers must be valid for the lifetime of this object!
            auto fun_write_param = [info, param_widget, min, max, this]() {
                // get values from line edit
                bool ok;
                int value         = util::qstring_to_int(param_widget->widgetLineEdit()->text(), &ok);
                bool out_of_range = false;

                // check limits
                if (ok) out_of_range = !corbo::util::is_in_bounds(value, min, max);

                // set values if ok or revert if not
                if (ok && !out_of_range)
                {
                    info.message_raw->GetReflection()->SetInt32(info.message_raw, info.field_raw, value);
                    emit parameterInt32Updated(
                        parentFieldNames() + "." + QString::fromStdString(MessageParser::nestedNamesToString(info.nested_names)), value);
                }
                else
                {
                    // revert to original value
                    param_widget->widgetLineEdit()->setText(
                        util::double_to_qstring(info.message_raw->GetReflection()->GetInt32(*info.message_raw, info.field_raw)));
                }

                // We need to check this after resetting to the previous value in order to prevent displaying the warning twice.
                if (out_of_range)
                    QMessageBox::warning(this, tr("Out of Range"),
                                         "Parameter must be in the interval [" + util::int_to_qstring(min) + ", " + util::int_to_qstring(max) + "]");
            };
            connect(param_widget->widgetLineEdit(), &QLineEdit::editingFinished, fun_write_param);
            break;
        }
    }
}

void ParameterWidget::addParameterInt32Array(const std::vector<int>& values, const MessageParser::FieldInformation& info)
{
    _has_parameters = true;

    messages::GuiType type = info.gui_type.first ? info.gui_type.second : messages::GuiType::TEXTEDIT;
    QString label          = info.label.first ? QString::fromStdString(info.label.second) : QString::fromStdString(info.field_name);
    QString description    = info.description.first ? QString::fromStdString(info.description.second) : "";
    int min                = info.min_value.first ? info.min_value.second : CORBO_MIN_INT;
    int max                = info.max_value.first ? info.max_value.second : CORBO_MAX_INT;

    std::vector<int> values_mod = values;

    // add field description label if specified
    if (info.print_description) addSubWidget(new QLabel(description));

    // check if we need to load default values
    if (info.default_value.first && !_groups.empty() && std::get<2>(_groups.back()))
    {
        // overwrite value
        bool ok = false;
        util::qstring_to_container(QString::fromStdString(info.default_value.second), values_mod, &ok);
        if (ok)
        {
            info.message_raw->GetReflection()->GetMutableRepeatedFieldRef<int>(info.message_raw, info.field_raw).CopyFrom(values_mod);
        }
        else
        {
            PRINT_ERROR_NAMED("Cannot convert default parameter: " << info.default_value.second << " to int array.");
            values_mod = values;
        }
    }

    switch (type)
    {
        case messages::GuiType::TEXTEDIT:
        default:
        {
            LabelEditWidget* param_widget = new LabelEditWidget(label, util::int_container_to_qstring(values_mod.begin(), values_mod.end()));
            param_widget->setToolTip(description);
            QRegExp regexp_no_decimal("^((?!\\.).)*$");
            param_widget->widgetLineEdit()->setValidator(new QRegExpValidator(regexp_no_decimal, param_widget->widgetLineEdit()));
            addSubWidget(param_widget);

            // set up modifer signal
            // message and field pointers must be valid for the lifetime of this object!
            auto fun_write_param = [info, param_widget, min, max, this]() {
                // get values from line edit
                bool ok;
                std::vector<int> values;
                util::qstring_to_container(param_widget->widgetLineEdit()->text(), values, &ok);

                // check dimension
                bool dimension_range_exceeded = false;
                int field_size                = -1;
                if (info.dynamic_size)
                {
                    field_size = (int)values.size();
                    if (info.min_size.first && field_size < info.min_size.second)
                    {
                        dimension_range_exceeded = true;
                        ok                       = false;
                    }
                    else if (info.max_size.first && field_size < info.max_size.second)
                    {
                        dimension_range_exceeded = true;
                        ok                       = false;
                    }
                }
                else
                {
                    int field_size = info.message_raw->GetReflection()->FieldSize(*info.message_raw, info.field_raw);
                    if (ok && (int)values.size() != field_size) ok = false;
                }

                // check limits
                bool out_of_range = false;
                if (ok) out_of_range = !corbo::util::is_in_bounds_all(values.begin(), values.end(), min, max);

                // set values if ok or revert if not
                if (ok && !out_of_range)
                {
                    if (info.dynamic_size)
                    {
                        info.message_raw->GetReflection()
                            ->GetMutableRepeatedFieldRef<int>(info.message_raw, info.field_raw)
                            .CopyFrom(values);  // this method resizes the field to the correct size
                    }
                    else
                    {
                        for (int idx = 0; idx < (int)values.size(); ++idx)
                        {
                            info.message_raw->GetReflection()->SetRepeatedInt32(info.message_raw, info.field_raw, idx, values[idx]);
                        }
                    }
                }
                else
                {
                    // revert to original value
                    auto orig_values = info.message_raw->GetReflection()->GetRepeatedFieldRef<int>(*info.message_raw, info.field_raw);
                    QString test     = util::int_container_to_qstring(orig_values.begin(), orig_values.end());
                    param_widget->widgetLineEdit()->setText(test);
                }

                // We need to check this after resetting to the previous value in order to prevent displaying the warning twice.
                if ((int)values.size() != field_size)
                {
                    QMessageBox::warning(this, tr("Invalid Number of Elements"),
                                         "Number of elements must be " + QString::number(field_size) + ", but " + QString::number(values.size()) +
                                             (values.size() == 1 ? " is" : " are") + " provided.");
                }
                else if (dimension_range_exceeded)
                {
                    QMessageBox::warning(this, tr("Invalid Number of Elements"),
                                         "Dimension must be between " + util::int_to_qstring(info.min_size.second) + " and " +
                                             util::int_to_qstring(info.max_size.second) + ".");
                }
                else if (out_of_range)
                {
                    QMessageBox::warning(this, tr("Out of Range"),
                                         "Parameters must be in the interval [" + util::int_to_qstring(min) + ", " + util::int_to_qstring(max) + "]");
                }
            };
            connect(param_widget->widgetLineEdit(), &QLineEdit::editingFinished, fun_write_param);
            break;
        }
    }
}

void ParameterWidget::addParameterDouble(double value, const MessageParser::FieldInformation& info)
{
    _has_parameters = true;

    messages::GuiType type = info.gui_type.first ? info.gui_type.second : messages::GuiType::TEXTEDIT;
    QString label          = info.label.first ? QString::fromStdString(info.label.second) : QString::fromStdString(info.field_name);
    QString description    = info.description.first ? QString::fromStdString(info.description.second) : "";
    double min             = info.min_value.first ? info.min_value.second : CORBO_MIN_DBL;
    double max             = info.max_value.first ? info.max_value.second : CORBO_MAX_DBL;

    // add field description label if specified
    if (info.print_description) addSubWidget(new QLabel(description));

    // check if we need to load a default value
    if (info.default_value.first && !_groups.empty() && std::get<2>(_groups.back()))
    {
        // overwrite value
        bool ok = false;
        value   = util::qstring_to_double(QString::fromStdString(info.default_value.second), &ok);
        if (!ok)
        {
            PRINT_ERROR_NAMED("Cannot convert default parameter: " << info.default_value.second << " to double.");
            value = 0;
        }
        info.message_raw->GetReflection()->SetDouble(info.message_raw, info.field_raw, value);
    }

    switch (type)
    {
        case messages::GuiType::SLIDER:
        {
            LabelSliderWidget* param_widget = new LabelSliderWidget(label, value, 1);
            param_widget->setToolTip(description);
            // if min and max are not specified, set some reasonable default limits:
            if (!info.min_value.first) min = 0;
            if (!info.max_value.first) max = 100;
            param_widget->setMinMax(min, max);
            addSubWidget(param_widget);

            // set up modifer signal
            // message and field pointers must be valid for the lifetime of this object!
            auto fun_write_param = [this, info](double value) {
                info.message_raw->GetReflection()->SetDouble(info.message_raw, info.field_raw, value);
                if (info.update_signals) emit signalUpdateRequested();
            };
            connect(param_widget, &LabelSliderWidget::valueChanged, fun_write_param);
            break;
        }
        case messages::GuiType::TEXTEDIT:
        default:
        {
            LabelEditWidget* param_widget = new LabelEditWidget(label, util::double_to_qstring(value));
            param_widget->setToolTip(description);
            addSubWidget(param_widget);

            // set up modifer signal
            // message and field pointers must be valid for the lifetime of this object!
            auto fun_write_param = [info, param_widget, min, max, this]() {
                // get values from line edit
                bool ok;
                double value      = util::qstring_to_double(param_widget->widgetLineEdit()->text(), &ok);
                bool out_of_range = false;

                // check limits
                if (ok) out_of_range = !corbo::util::is_in_bounds(value, min, max);

                // set values if ok or revert if not
                if (ok && !out_of_range)
                    info.message_raw->GetReflection()->SetDouble(info.message_raw, info.field_raw, value);
                else
                {
                    // revert to original value
                    param_widget->widgetLineEdit()->setText(
                        util::double_to_qstring(info.message_raw->GetReflection()->GetDouble(*info.message_raw, info.field_raw)));
                }

                // We need to check this after resetting to the previous value in order to prevent displaying the warning twice.
                if (out_of_range)
                    QMessageBox::warning(
                        this, tr("Out of Range"),
                        "Parameter must be in the interval [" + util::double_to_qstring(min) + ", " + util::double_to_qstring(max) + "]");

                if (info.update_signals) emit signalUpdateRequested();
            };
            connect(param_widget->widgetLineEdit(), &QLineEdit::editingFinished, fun_write_param);
            break;
        }
    }
}

void ParameterWidget::addParameterDoubleArray(const std::vector<double>& values, const MessageParser::FieldInformation& info)
{
    _has_parameters = true;

    messages::GuiType type = info.gui_type.first ? info.gui_type.second : messages::GuiType::TEXTEDIT;
    QString label          = info.label.first ? QString::fromStdString(info.label.second) : QString::fromStdString(info.field_name);
    QString description    = info.description.first ? QString::fromStdString(info.description.second) : "";
    double min             = info.min_value.first ? info.min_value.second : CORBO_MIN_DBL;
    double max             = info.max_value.first ? info.max_value.second : CORBO_MAX_DBL;

    std::vector<double> values_mod = values;

    // add field description label if specified
    if (info.print_description) addSubWidget(new QLabel(description));

    // check if we need to load default values
    if (info.default_value.first && !_groups.empty() && std::get<2>(_groups.back()))
    {
        // overwrite value
        bool ok = false;
        util::qstring_to_container(QString::fromStdString(info.default_value.second), values_mod, &ok);
        if (ok)
        {
            info.message_raw->GetReflection()->GetMutableRepeatedFieldRef<double>(info.message_raw, info.field_raw).CopyFrom(values_mod);
        }
        else
        {
            PRINT_ERROR_NAMED("Cannot convert default parameter: " << info.default_value.second << " to double array.");
            values_mod = values;
        }
    }

    switch (type)
    {
        case messages::GuiType::TEXTEDIT:
        default:
        {
            LabelEditWidget* param_widget = new LabelEditWidget(label, util::double_container_to_qstring(values_mod.begin(), values_mod.end()));
            param_widget->setToolTip(description);
            addSubWidget(param_widget);

            // set up modifer signal
            // message and field pointers must be valid for the lifetime of this object!
            auto fun_write_param = [info, param_widget, min, max, this]() {
                // get values from line edit
                bool ok;
                std::vector<double> values;
                util::qstring_to_container(param_widget->widgetLineEdit()->text(), values, &ok);

                // check dimension
                int field_size                = -1;
                bool dimension_range_exceeded = false;
                if (info.dynamic_size)
                {
                    field_size = (int)values.size();
                    if (info.min_size.first && field_size < info.min_size.second)
                    {
                        dimension_range_exceeded = true;
                        ok                       = false;
                    }
                    else if (info.max_size.first && field_size < info.max_size.second)
                    {
                        dimension_range_exceeded = true;
                        ok                       = false;
                    }
                }
                else
                {
                    field_size = info.message_raw->GetReflection()->FieldSize(*info.message_raw, info.field_raw);
                    if (ok && (int)values.size() != field_size) ok = false;
                }

                // check limits
                bool out_of_range = false;
                if (ok) out_of_range = !corbo::util::is_in_bounds_all(values.begin(), values.end(), min, max);

                // set values if ok or revert if not
                if (ok && !out_of_range)
                {
                    if (info.dynamic_size)
                    {
                        info.message_raw->GetReflection()
                            ->GetMutableRepeatedFieldRef<double>(info.message_raw, info.field_raw)
                            .CopyFrom(values);  // this method resizes the field to the correct size
                    }
                    else
                    {
                        for (int idx = 0; idx < (int)values.size(); ++idx)
                        {
                            info.message_raw->GetReflection()->SetRepeatedDouble(info.message_raw, info.field_raw, idx, values[idx]);
                        }
                    }
                }
                else
                {
                    // revert to original value
                    auto orig_values = info.message_raw->GetReflection()->GetRepeatedFieldRef<double>(*info.message_raw, info.field_raw);
                    QString text     = util::double_container_to_qstring(orig_values.begin(), orig_values.end());
                    param_widget->widgetLineEdit()->setText(text);
                }

                // We need to check this after resetting to the previous value in order to prevent displaying the warning twice.
                if ((int)values.size() != field_size)
                {
                    QMessageBox::warning(this, tr("Invalid Number of Elements"),
                                         "Number of elements must be " + QString::number(field_size) + ", but " + QString::number(values.size()) +
                                             (values.size() == 1 ? " is" : " are") + " provided.");
                }
                else if (dimension_range_exceeded)
                {
                    QMessageBox::warning(this, tr("Invalid Number of Elements"),
                                         "Dimension must be between " + util::int_to_qstring(info.min_size.second) + " and " +
                                             util::int_to_qstring(info.max_size.second) + ".");
                }
                else if (out_of_range)
                {
                    QMessageBox::warning(
                        this, tr("Out of Range"),
                        "Parameters must be in the interval [" + util::double_to_qstring(min) + ", " + util::double_to_qstring(max) + "]");
                }
                if (info.update_signals) emit signalUpdateRequested();
            };
            connect(param_widget->widgetLineEdit(), &QLineEdit::editingFinished, fun_write_param);
            break;
        }
    }
}

void ParameterWidget::addParameterBool(bool value, const MessageParser::FieldInformation& info)
{
    _has_parameters = true;

    messages::GuiType type = info.gui_type.first ? info.gui_type.second : messages::GuiType::CHECKBOX;
    QString label          = info.label.first ? QString::fromStdString(info.label.second) : QString::fromStdString(info.field_name);
    QString description    = info.description.first ? QString::fromStdString(info.description.second) : "";

    // add field description label if specified
    if (info.print_description) addSubWidget(new QLabel(description));

    // check if we need to load a default value
    if (info.default_value.first && !_groups.empty() && std::get<2>(_groups.back()))
    {
        // overwrite value
        bool ok = false;
        value   = util::qstring_to_bool(QString::fromStdString(info.default_value.second), &ok);
        if (!ok)
        {
            PRINT_ERROR_NAMED("Cannot convert default parameter: " << info.default_value.second << " to bool.");
            value = 0;
        }
        info.message_raw->GetReflection()->SetBool(info.message_raw, info.field_raw, value);
    }

    switch (type)
    {
        case messages::GuiType::CHECKBOX:
        default:
        {
            QCheckBox* param_widget = new QCheckBox(label);
            param_widget->setToolTip(description);
            param_widget->setChecked(value);
            addSubWidget(param_widget);

            // set up modifer signal
            // message and field pointers must be valid for the lifetime of this object!
            auto fun_write_param = [info, this](bool toggled) {
                info.message_raw->GetReflection()->SetBool(info.message_raw, info.field_raw, toggled);
                if (info.update_signals) emit signalUpdateRequested();
            };
            connect(param_widget, &QCheckBox::toggled, fun_write_param);

            break;
        }
    }
}

void ParameterWidget::addParameterBoolArray(const std::vector<bool>& values, const MessageParser::FieldInformation& info)
{
    _has_parameters = true;

    messages::GuiType type = info.gui_type.first ? info.gui_type.second : messages::GuiType::CHECKBOX;
    QString label          = info.label.first ? QString::fromStdString(info.label.second) : QString::fromStdString(info.field_name);
    QString description    = info.description.first ? QString::fromStdString(info.description.second) : "";

    std::vector<bool> values_mod = values;

    // add field description label if specified
    if (info.print_description) addSubWidget(new QLabel(description));

    // check if we need to load default values
    if (info.default_value.first && !_groups.empty() && std::get<2>(_groups.back()))
    {
        // overwrite value
        bool ok = false;
        util::qstring_to_container(QString::fromStdString(info.default_value.second), values_mod, &ok);
        if (ok)
        {
            info.message_raw->GetReflection()->GetMutableRepeatedFieldRef<bool>(info.message_raw, info.field_raw).CopyFrom(values_mod);
        }
        else
        {
            PRINT_ERROR_NAMED("Cannot convert default parameter: " << info.default_value.second << " to bool array.");
            values_mod = values;
        }
    }

    // check if we have multiple labels separated by comma
    QStringList label_list = label.split(',');

    QString group_name = label_list.empty() ? label : label_list.first();

    if (type == messages::GuiType::TEXTEDIT)
    {
        // text edit version
        LabelEditWidget* param_widget = new LabelEditWidget(label, util::bool_container_to_qstring(values_mod.begin(), values_mod.end()));
        param_widget->setToolTip(description);
        addSubWidget(param_widget);

        // set up modifer signal
        // message and field pointers must be valid for the lifetime of this object!
        auto fun_write_param = [info, param_widget, this]() {
            // get values from line edit
            bool ok;
            std::vector<bool> values;
            util::qstring_to_container(param_widget->widgetLineEdit()->text(), values, &ok);

            // check dimension
            int field_size                = -1;
            bool dimension_range_exceeded = false;
            if (info.dynamic_size)
            {
                field_size = (int)values.size();
                if (info.min_size.first && field_size < info.min_size.second)
                {
                    dimension_range_exceeded = true;
                    ok                       = false;
                }
                else if (info.max_size.first && field_size < info.max_size.second)
                {
                    dimension_range_exceeded = true;
                    ok                       = false;
                }
            }
            else
            {
                field_size = info.message_raw->GetReflection()->FieldSize(*info.message_raw, info.field_raw);
                if (ok && (int)values.size() != field_size) ok = false;
            }

            // set values if ok or revert if not
            if (ok)
            {
                if (info.dynamic_size)
                {
                    info.message_raw->GetReflection()
                        ->GetMutableRepeatedFieldRef<bool>(info.message_raw, info.field_raw)
                        .CopyFrom(values);  // this method resizes the field to the correct size
                }
                else
                {
                    for (int idx = 0; idx < (int)values.size(); ++idx)
                    {
                        info.message_raw->GetReflection()->SetRepeatedBool(info.message_raw, info.field_raw, idx, values[idx]);
                    }
                }
            }
            else
            {
                // revert to original value
                auto orig_values = info.message_raw->GetReflection()->GetRepeatedFieldRef<bool>(*info.message_raw, info.field_raw);
                QString text     = util::bool_container_to_qstring(orig_values.begin(), orig_values.end());
                param_widget->widgetLineEdit()->setText(text);
            }

            // We need to check this after resetting to the previous value in order to prevent displaying the warning twice.
            if ((int)values.size() != field_size)
            {
                QMessageBox::warning(this, tr("Invalid Number of Elements"),
                                     "Number of elements must be " + QString::number(field_size) + ", but " + QString::number(values.size()) +
                                         (values.size() == 1 ? " is" : " are") + " provided.");
            }
            else if (dimension_range_exceeded)
            {
                QMessageBox::warning(this, tr("Invalid Number of Elements"),
                                     "Dimension must be between " + util::int_to_qstring(info.min_size.second) + " and " +
                                         util::int_to_qstring(info.max_size.second) + ".");
            }
            if (info.update_signals) emit signalUpdateRequested();
        };
        connect(param_widget->widgetLineEdit(), &QLineEdit::editingFinished, fun_write_param);
    }
    else
    {
        // horizontal button group
        bool exclusive = type != messages::GuiType::CHECKBOX;  // only if we have a checkbox, multiple checked booleans are allowed

        HoriontalButtonGroup* param_widget = new HoriontalButtonGroup(group_name, exclusive);
        param_widget->setToolTip(description);

        for (int idx = 0; idx < (int)values_mod.size(); ++idx)
        {
            QString field_label  = (idx + 1 < label_list.size()) ? label_list[idx + 1] : "";  // first element is the group name!
            QAbstractButton* btn = param_widget->addButton(values_mod[idx], field_label);

            // set up modifer signal
            // message and field pointers must be valid for the lifetime of this object!
            auto fun_write_param = [info, idx, this](bool toggled) {
                info.message_raw->GetReflection()->SetRepeatedBool(info.message_raw, info.field_raw, idx, toggled);
                if (info.update_signals) emit signalUpdateRequested();
            };
            connect(btn, &QAbstractButton::toggled, fun_write_param);
        }
        addSubWidget(param_widget);
    }
}

void ParameterWidget::addParameterEnum(const std::string& value, const MessageParser::FieldInformation& info)
{
    _has_parameters = true;

    QString label       = info.label.first ? QString::fromStdString(info.label.second) : QString::fromStdString(info.field_name);
    QString description = info.description.first ? QString::fromStdString(info.description.second) : "";

    // add field description label if specified
    if (info.print_description) addSubWidget(new QLabel(description));

    LabelComboBoxWidget* param_widget = new LabelComboBoxWidget(label);
    param_widget->setToolTip(description);
    addSubWidget(param_widget);

    // search for all enum types
    const google::protobuf::EnumDescriptor* enum_descr = info.field_raw->enum_type();
    if (!enum_descr)
    {
        PRINT_WARNING("ParameterWidget::addParameterEnum(): provided parameter '" << info.field_name << "' is no enumerator!");
        return;
    }
    for (int idx = 0; idx < enum_descr->value_count(); ++idx)
    {
        param_widget->widgetComboBox()->addItem(QString::fromStdString(enum_descr->value(idx)->name()));
    }
    // select current value
    param_widget->widgetComboBox()->setCurrentText(QString::fromStdString(value));

    // set up modifer signal
    // message and field pointers must be valid for the lifetime of this object!
    auto fun_write_param = [info, this](int idx) {
        info.message_raw->GetReflection()->SetEnumValue(info.message_raw, info.field_raw, idx);
        if (info.update_signals) emit signalUpdateRequested();
    };
    connect(param_widget->widgetComboBox(), static_cast<void (QComboBox::*)(int)>(&QComboBox::currentIndexChanged), fun_write_param);
}

void ParameterWidget::addParameterString(std::string value, const MessageParser::FieldInformation& info)
{
    _has_parameters = true;

    QString label       = info.label.first ? QString::fromStdString(info.label.second) : QString::fromStdString(info.field_name);
    QString description = info.description.first ? QString::fromStdString(info.description.second) : "";

    // add field description label if specified
    if (info.print_description) addSubWidget(new QLabel(description));

    // check if we need to load a default value
    if (info.default_value.first && !_groups.empty() && std::get<2>(_groups.back()))
    {
        // overwrite value
        value = info.default_value.second;

        info.message_raw->GetReflection()->SetString(info.message_raw, info.field_raw, value);
    }

    LabelEditWidget* param_widget = new LabelEditWidget(label, QString::fromStdString(value));
    param_widget->setToolTip(description);
    addSubWidget(param_widget);

    // set up modifer signal
    // message and field pointers must be valid for the lifetime of this object!
    auto fun_write_param = [param_widget, info, this]() {
        std::string value = param_widget->widgetLineEdit()->text().toUtf8().constData();
        info.message_raw->GetReflection()->SetString(info.message_raw, info.field_raw, value);
        if (info.update_signals) emit signalUpdateRequested();
    };
    connect(param_widget->widgetLineEdit(), &QLineEdit::editingFinished, fun_write_param);
}

void ParameterWidget::addInfoText(const std::string& text, const MessageParser::FieldInformation& info)
{
    _has_parameters = true;

    QString text_q = QString::fromStdString(text);

    if (info.default_value.first && text.empty()) text_q = QString::fromStdString(info.default_value.second);

    addSubWidget(new QLabel(text_q));
}

void ParameterWidget::startGroup(const MessageParser::FieldInformation& info)
{
    // start group
    QString name               = info.label.first ? QString::fromStdString(info.label.second) : QString::fromStdString(info.field_name);
    CollapsableGroupBox* group = new CollapsableGroupBox(name);
    group->setCollapsed(info.collapsed);

    QVBoxLayout* layout = new QVBoxLayout;
    layout->setContentsMargins(5, 15, 5, 0);
    layout->setAlignment(Qt::AlignTop | Qt::AlignLeft);
    _groups.push(std::make_tuple(group, layout, info.new_message));
}
void ParameterWidget::endGroup()
{
    if (_groups.empty()) return;

    std::tuple<CollapsableGroupBox*, QLayout*, bool> top = _groups.pop();

    // apply layout to group
    std::get<0>(top)->groupBox()->setLayout(std::get<1>(top));
    std::get<0>(top)->setCollapsed(false);

    // append to parent group or widget
    if (_groups.empty())
        _layout->addWidget(std::get<0>(top));
    else
        std::get<1>(_groups.top())->addWidget(std::get<0>(top));
}

void ParameterWidget::addSubWidget(QWidget* widget)
{
    // set font size for parameters
    // QFont font = widget->font();
    // font.setPointSize(11);
    // widget->setFont(font);
    widget->setStyleSheet("font-size: 11px");  // hereby, we also set the font size of child widgets!

    if (_groups.empty())
    {
        _layout->addWidget(widget);
        // TODO(roesmann) this command is called quite often with the same widget!!!
        // _layout->addWidget(new QLabel("test", this));
    }
    else if (std::get<1>(_groups.top()))
    {
        std::get<1>(_groups.top())->addWidget(widget);
        // _layout->addWidget(new QLabel("test2", this));
    }
}

void ParameterWidget::generateFromMessage(std::shared_ptr<google::protobuf::Message> message)
{
    clearElements();
    _param_message = message;
    generateElements();
}

void ParameterWidget::generateFromAllocatedField(google::protobuf::Message* message, const google::protobuf::FieldDescriptor* field)
{
    clearElements();
    parse(message, field);
}

void ParameterWidget::generateElements()
{
    if (!_param_message) return;
    parse(_param_message.get());
}

void ParameterWidget::parse(google::protobuf::Message* message, const google::protobuf::FieldDescriptor* field)
{
    MessageParser parser;

    using std::placeholders::_1;
    using std::placeholders::_2;

    parser.setCallbackValueInt32(std::bind(&ParameterWidget::addParameterInt32, this, _1, _2));
    parser.setCallbackValueInt32Array(std::bind(&ParameterWidget::addParameterInt32Array, this, _1, _2));
    parser.setCallbackValueDouble(std::bind(&ParameterWidget::addParameterDouble, this, _1, _2));
    parser.setCallbackValueDoubleArray(std::bind(&ParameterWidget::addParameterDoubleArray, this, _1, _2));
    parser.setCallbackValueBool(std::bind(&ParameterWidget::addParameterBool, this, _1, _2));
    parser.setCallbackValueBoolArray(std::bind(&ParameterWidget::addParameterBoolArray, this, _1, _2));
    parser.setCallbackValueEnum(std::bind(&ParameterWidget::addParameterEnum, this, _1, _2));
    parser.setCallbackValueString(std::bind(&ParameterWidget::addParameterString, this, _1, _2));
    parser.setCallbackValueStringInfo(std::bind(&ParameterWidget::addInfoText, this, _1, _2));

    parser.setCallbackMessageEvent([this, field](corbo::MessageParser::MessageEvent ev, const corbo::MessageParser::FieldInformation& info) {
        if (ev == corbo::MessageParser::MessageEvent::MessageStart)
            startGroup(info);
        else if (ev == corbo::MessageParser::MessageEvent::MessageEnd)
            endGroup();
        else if (ev == corbo::MessageParser::MessageEvent::OneOf)
        {
            if (info.field_raw != field) addOneOfField(info);
        }
    });

    std::list<std::string> nested_field_name;  // required for storing nested field names while parsing

    // prepend any previous parent field names
    nested_field_name = _nested_parent_fields;
    if (!nested_field_name.empty()) nested_field_name.pop_back();  // erase last to avoid duplicate

    if (field)
    {
        parser.parseField(message, field, false, true, false,
                          &nested_field_name);  // options: parse first oneof even if not active since we want to force parse the provided field
    }
    else
        parser.parse(message, false, true, &nested_field_name);
}

bool ParameterWidget::hasParameters() const
{
    if (_has_parameters) return true;

    // also check child widgets (subgroups)
    for (auto widget : this->findChildren<ParameterWidget*>(QString(), Qt::FindChildrenRecursively))
    {
        if (static_cast<ParameterWidget*>(widget)->hasParameters()) return true;
    }
    return false;
}

void ParameterWidget::clearElements()
{
    for (auto widget : this->findChildren<QWidget*>(QString(), Qt::FindDirectChildrenOnly)) delete widget;

    _groups.clear();
    _oneof_widgets.clear();
}

}  // namespace gui
}  // namespace corbo
