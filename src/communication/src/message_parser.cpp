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

#ifdef MESSAGE_SUPPORT

#include <corbo-communication/message_parser.h>
#include <corbo-core/console.h>

#include <limits>
#include <list>
#include <string>
#include <utility>
#include <vector>

namespace p = google::protobuf;

namespace corbo {

void MessageParser::consumeFieldMessage(p::Message* message, const p::Reflection* reflection, const p::FieldDescriptor* field, bool expand_oneof,
                                        bool skip_inactive_oneof, std::list<std::string>* nested_names)
{
    FieldInformation info;
    getFieldInformation(message, field, info, nested_names, true);

    if (field->is_repeated())
    {
        for (int i = 0; i < reflection->FieldSize(*message, field); ++i)
        {
            p::Message* new_message = reflection->MutableRepeatedMessage(message, field, i);

            if (_fun_msg_event) _fun_msg_event(MessageEvent::MessageStart, info);

            parse(new_message, expand_oneof, skip_inactive_oneof, nested_names);

            if (_fun_msg_event) _fun_msg_event(MessageEvent::MessageEnd, info);
        }
    }
    else
    {
        p::Message* new_message = reflection->MutableMessage(message, field);

        if (_fun_msg_event) _fun_msg_event(MessageEvent::MessageStart, info);

        parse(new_message, expand_oneof, skip_inactive_oneof, nested_names);

        if (_fun_msg_event) _fun_msg_event(MessageEvent::MessageEnd, info);
    }
}

void MessageParser::consumeFieldValue(p::Message* message, const p::Reflection* reflection, const p::FieldDescriptor* field,
                                      std::list<std::string>* nested_names)
{
    switch (field->cpp_type())
    {
        case p::FieldDescriptor::CPPTYPE_INT32:
        {
            if (field->is_repeated())
            {
                std::vector<std::int32_t> values;
                for (int i = 0; i < reflection->FieldSize(*message, field); ++i) values.push_back(reflection->GetRepeatedInt32(*message, field, i));

                if (_fun_value_int32_array)
                {
                    FieldInformation info;
                    getFieldInformation(message, field, info, nested_names, false);
                    _fun_value_int32_array(values, info);
                }
                else
                    PRINT_WARNING_COND(_verbose, "MessageParser: found value '" << field->name() << "' [int32 array] but no callback is registered.");
            }
            else
            {
                std::int32_t value = reflection->GetInt32(*message, field);
                if (_fun_value_int32)
                {
                    FieldInformation info;
                    getFieldInformation(message, field, info, nested_names, false);
                    _fun_value_int32(value, info);
                }
                else
                    PRINT_WARNING_COND(_verbose, "MessageParser: found value '" << field->name() << "' [int32] but no callback is registered.");
            }
            break;
        }
        case p::FieldDescriptor::CPPTYPE_UINT32:
        {
            if (field->is_repeated())
            {
                std::vector<std::uint32_t> values;
                for (int i = 0; i < reflection->FieldSize(*message, field); ++i) values.push_back(reflection->GetRepeatedUInt32(*message, field, i));

                if (_fun_value_uint32_array)
                {
                    FieldInformation info;
                    getFieldInformation(message, field, info, nested_names, false);
                    _fun_value_uint32_array(values, info);
                }
                else
                    PRINT_WARNING_COND(_verbose,
                                       "MessageParser: found value '" << field->name() << "' [uint32 array] but no callback is registered.");
            }
            else
            {
                std::uint32_t value = reflection->GetUInt32(*message, field);
                if (_fun_value_uint32)
                {
                    FieldInformation info;
                    getFieldInformation(message, field, info, nested_names, false);
                    _fun_value_uint32(value, info);
                }
                else
                    PRINT_WARNING_COND(_verbose, "MessageParser: found value '" << field->name() << "' [uint32] but no callback is registered.");
            }
            break;
        }
        case p::FieldDescriptor::CPPTYPE_INT64:
        {
            if (field->is_repeated())
            {
                std::vector<std::int64_t> values;
                for (int i = 0; i < reflection->FieldSize(*message, field); ++i) values.push_back(reflection->GetRepeatedInt64(*message, field, i));

                if (_fun_value_int64_array)
                {
                    FieldInformation info;
                    getFieldInformation(message, field, info, nested_names, false);
                    _fun_value_int64_array(values, info);
                }
                else
                    PRINT_WARNING_COND(_verbose, "MessageParser: found value '" << field->name() << "' [int64 array] but no callback is registered.");
            }
            else
            {
                std::int64_t value = reflection->GetInt64(*message, field);
                if (_fun_value_int64)
                {
                    FieldInformation info;
                    getFieldInformation(message, field, info, nested_names, false);
                    _fun_value_int64(value, info);
                }
                else
                    PRINT_WARNING_COND(_verbose, "MessageParser: found value '" << field->name() << "' [int64] but no callback is registered.");
            }
            break;
        }
        case p::FieldDescriptor::CPPTYPE_UINT64:
        {
            if (field->is_repeated())
            {
                std::vector<std::uint64_t> values;
                for (int i = 0; i < reflection->FieldSize(*message, field); ++i) values.push_back(reflection->GetRepeatedUInt64(*message, field, i));

                if (_fun_value_uint64_array)
                {
                    FieldInformation info;
                    getFieldInformation(message, field, info, nested_names, false);
                    _fun_value_uint64_array(values, info);
                }
                else
                    PRINT_WARNING_COND(_verbose,
                                       "MessageParser: found value '" << field->name() << "' [uint64 array] but no callback is registered.");
            }
            else
            {
                std::uint64_t value = reflection->GetUInt64(*message, field);
                if (_fun_value_uint64)
                {
                    FieldInformation info;
                    getFieldInformation(message, field, info, nested_names, false);
                    _fun_value_uint64(value, info);
                }
                else
                    PRINT_WARNING_COND(_verbose, "MessageParser: found value '" << field->name() << "' [uint64] but no callback is registered.");
            }
            break;
        }
        case p::FieldDescriptor::CPPTYPE_FLOAT:
        {
            if (field->is_repeated())
            {
                std::vector<float> values;
                for (int i = 0; i < reflection->FieldSize(*message, field); ++i) values.push_back(reflection->GetRepeatedFloat(*message, field, i));

                if (_fun_value_float_array)
                {
                    FieldInformation info;
                    getFieldInformation(message, field, info, nested_names, false);
                    _fun_value_float_array(values, info);
                }
                else
                    PRINT_WARNING_COND(_verbose, "MessageParser: found value '" << field->name() << "' [float array] but no callback is registered.");
            }
            else
            {
                float value = reflection->GetFloat(*message, field);
                if (_fun_value_float)
                {
                    FieldInformation info;
                    getFieldInformation(message, field, info, nested_names, false);
                    _fun_value_float(value, info);
                }
                else
                    PRINT_WARNING_COND(_verbose, "MessageParser: found value '" << field->name() << "' [float] but no callback is registered.");
            }
            break;
        }
        case p::FieldDescriptor::CPPTYPE_DOUBLE:
        {
            if (field->is_repeated())
            {
                std::vector<double> values;
                for (int i = 0; i < reflection->FieldSize(*message, field); ++i) values.push_back(reflection->GetRepeatedDouble(*message, field, i));

                if (_fun_value_double_array)
                {
                    FieldInformation info;
                    getFieldInformation(message, field, info, nested_names, false);
                    _fun_value_double_array(values, info);
                }
                else
                    PRINT_WARNING_COND(_verbose,
                                       "MessageParser: found value '" << field->name() << "' [double array] but no callback is registered.");
            }
            else
            {
                double value = reflection->GetDouble(*message, field);
                if (_fun_value_double)
                {
                    FieldInformation info;
                    getFieldInformation(message, field, info, nested_names, false);
                    _fun_value_double(value, info);
                }
                else
                    PRINT_WARNING_COND(_verbose, "MessageParser: found value '" << field->name() << "' [double] but no callback is registered.");
            }
            break;
        }
        case p::FieldDescriptor::CPPTYPE_STRING:
        {
            if (field->is_repeated())
            {
                std::vector<std::string> values;
                for (int i = 0; i < reflection->FieldSize(*message, field); ++i) values.push_back(reflection->GetRepeatedString(*message, field, i));

                if (_fun_value_string_array)
                {
                    FieldInformation info;
                    getFieldInformation(message, field, info, nested_names, false);
                    _fun_value_string_array(values, info);
                }
                else
                    PRINT_WARNING_COND(_verbose,
                                       "MessageParser: found value '" << field->name() << "' [string array] but no callback is registered.");
            }
            else
            {
                std::string text = reflection->GetString(*message, field);
                FieldInformation info;
                getFieldInformation(message, field, info, nested_names, false);

                if (info.info_text.first && info.info_text.second)
                {
                    if (_fun_value_string_info)
                        _fun_value_string_info(text, info);
                    else
                        PRINT_WARNING_COND(_verbose,
                                           "MessageParser: found value '" << field->name() << "' [string info] but no callback is registered.");
                }
                else
                {
                    if (_fun_value_string)
                        _fun_value_string(text, info);
                    else
                        PRINT_WARNING_COND(_verbose, "MessageParser: found value '" << field->name() << "' [string] but no callback is registered.");
                }
            }
            break;
        }
        case p::FieldDescriptor::CPPTYPE_BOOL:
        {
            if (field->is_repeated())
            {
                std::vector<bool> values;
                for (int i = 0; i < reflection->FieldSize(*message, field); ++i) values.push_back(reflection->GetRepeatedBool(*message, field, i));

                if (_fun_value_bool_array)
                {
                    FieldInformation info;
                    getFieldInformation(message, field, info, nested_names, false);
                    _fun_value_bool_array(values, info);
                }
                else
                    PRINT_WARNING_COND(_verbose, "MessageParser: found value '" << field->name() << "' [bool array] but no callback is registered.");
            }
            else
            {
                bool flag = reflection->GetBool(*message, field);
                if (_fun_value_bool)
                {
                    FieldInformation info;
                    getFieldInformation(message, field, info, nested_names, false);
                    _fun_value_bool(flag, info);
                }
                else
                    PRINT_WARNING_COND(_verbose, "MessageParser: found value '" << field->name() << "' [bool] but no callback is registered.");
            }
            break;
        }
        case p::FieldDescriptor::CPPTYPE_ENUM:
        {
            if (field->is_repeated())
            {
                std::vector<std::string> values;
                for (int i = 0; i < reflection->FieldSize(*message, field); ++i)
                {
                    const p::EnumValueDescriptor* enum_desc = reflection->GetRepeatedEnum(*message, field, i);
                    values.push_back(enum_desc->name());
                }

                if (_fun_value_enum_array)
                {
                    FieldInformation info;
                    getFieldInformation(message, field, info, nested_names, false);
                    _fun_value_enum_array(values, info);
                }
                else
                    PRINT_WARNING_COND(_verbose, "MessageParser: found value '" << field->name() << "' [enum array] but no callback is registered.");
            }
            else
            {
                const p::EnumValueDescriptor* enum_desc = reflection->GetEnum(*message, field);

                std::string name = enum_desc->name();
                if (_fun_value_enum)
                {
                    FieldInformation info;
                    getFieldInformation(message, field, info, nested_names, false);
                    _fun_value_enum(name, info);
                }
                else
                    PRINT_WARNING_COND(_verbose, "MessageParser: found value '" << field->name() << "' [enum] but no callback is registered.");
            }
            break;
        }
        case p::FieldDescriptor::CPPTYPE_MESSAGE:
        {
            PRINT_WARNING_COND(_verbose,
                               "MessageParser: CPPTYPE_MESSAGE '" << field->name() << "' should be handeld in a different parser method...");
            break;
        }
        default:
        {
            PRINT_WARNING_COND(_verbose, "MessageParser: Field '" << field->name() << "' unknown / not implemented.");
            break;
        }
    }
}

void MessageParser::parseField(p::Message* message, const p::FieldDescriptor* field, bool expand_oneof, bool expand_first_oneof,
                               bool skip_inactive_oneof, std::list<std::string>* nested_names)
{
    if (nested_names) nested_names->push_back(field->name());

    bool expand_field = true;

    // check for one-of field
    const p::OneofDescriptor* oneof = field->containing_oneof();
    if (oneof)
    {
        // inform about oneof
        if (_fun_msg_event)
        {
            FieldInformation info;
            getFieldInformation(message, field, info, nested_names);
            _fun_msg_event(MessageEvent::OneOf, info);
        }
        if (expand_oneof || expand_first_oneof)
        {
            if (message->GetReflection()->HasOneof(*message, oneof))
            {
                const p::FieldDescriptor* oneof_field = message->GetReflection()->GetOneofFieldDescriptor(*message, oneof);
                // skip if this is not the active oneof-field
                if (skip_inactive_oneof && field != oneof_field) expand_field = false;
            }
            else
            {
                if (skip_inactive_oneof) expand_field = false;  // skip also if no oneof is active
            }
        }
        else
            expand_field = false;
    }

    if (expand_field)
    {
        // check if we have a separate message or a value
        if (field->cpp_type() == p::FieldDescriptor::CPPTYPE_MESSAGE)
        {
            consumeFieldMessage(message, message->GetReflection(), field, expand_oneof, skip_inactive_oneof, nested_names);
        }
        else
        {
            consumeFieldValue(message, message->GetReflection(), field, nested_names);
        }
    }

    if (nested_names) nested_names->pop_back();
}

void MessageParser::parse(p::Message* message, bool expand_oneof, bool skip_inactive_oneof, std::list<std::string>* nested_names)
{
    // iterate
    const p::Descriptor* desc = message->GetDescriptor();
    // const p::Reflection* refl = message.GetReflection();

    int num_fields = desc->field_count();
    for (int i = 0; i < num_fields; ++i)
    {
        const p::FieldDescriptor* field = desc->field(i);
        parseField(message, field, expand_oneof, expand_oneof, skip_inactive_oneof, nested_names);
    }
}

void MessageParser::getFieldInformation(p::Message* message, const p::FieldDescriptor* field, FieldInformation& info,
                                        std::list<std::string>* nested_names, bool is_message)
{
    // extract name
    info.field_name = field->name();
    info.msg_type   = message->GetTypeName();

    // check if group is new
    if (is_message)
    {
        if (field->is_repeated())
        {
            info.new_message = message->GetReflection()->FieldSize(*message, field) == 0;
        }
        else
        {
            info.new_message = !message->GetReflection()->HasField(*message, field);
        }
    }
    else
        info.new_message = false;

    // store nested names if provided
    if (nested_names) info.nested_names = *nested_names;

    // extract one-of descriptor (nullptr if no oneof)
    info.oneof_descriptor = field->containing_oneof();
    if (info.oneof_descriptor)
        info.oneof_selected_field = message->GetReflection()->GetOneofFieldDescriptor(*message, info.oneof_descriptor);  // returns 0 if not set

    // extract enum fields if field defines ane enumeration
    const p::EnumDescriptor* enum_desc = field->enum_type();
    if (enum_desc)
    {
        for (int i = 0; i < enum_desc->value_count(); ++i)
        {
            info.enum_list.push_back(enum_desc->value(i)->name());
        }
    }

    // support reflection api from outside and hence store relevant information:
    info.field_raw   = field;
    info.message_raw = message;

    // check field options
    const p::FieldOptions& options = field->options();

    // default value:
    if (options.HasExtension(messages::default_value))
        info.default_value = std::make_pair(true, options.GetExtension(messages::default_value));
    else
        info.default_value = std::make_pair(false, 0.0);

    // min value:
    if (options.HasExtension(messages::min_value))
        info.min_value = std::make_pair(true, options.GetExtension(messages::min_value));
    else
        info.min_value = std::make_pair(false, 0.0);

    // max value:
    if (options.HasExtension(messages::max_value))
        info.max_value = std::make_pair(true, options.GetExtension(messages::max_value));
    else
        info.max_value = std::make_pair(false, 0.0);

    // dynamic size:
    if (options.HasExtension(messages::dynamic_size))
        info.dynamic_size = options.GetExtension(messages::dynamic_size);
    else
        info.dynamic_size = false;

    // min size:
    if (options.HasExtension(messages::min_size))
        info.min_size = std::make_pair(true, options.GetExtension(messages::min_size));
    else
        info.min_size = std::make_pair(false, 0);

    // max size:
    if (options.HasExtension(messages::max_size))
        info.max_size = std::make_pair(true, options.GetExtension(messages::max_size));
    else
        info.max_size = std::make_pair(false, std::numeric_limits<int>::max());

    // info flag:
    if (options.HasExtension(messages::info))
        info.info_text = std::make_pair(true, options.GetExtension(messages::info));
    else
        info.info_text = std::make_pair(false, false);

    // description
    if (options.HasExtension(messages::description))
        info.description = std::make_pair(true, options.GetExtension(messages::description));
    else
        info.description = std::make_pair(false, "");

    // additional label:
    if (options.HasExtension(messages::label))
        info.label = std::make_pair(true, options.GetExtension(messages::label));
    else
        info.label = std::make_pair(false, "");

    // preferred gui type
    if (options.HasExtension(messages::gui_type))
        info.gui_type = std::make_pair(true, options.GetExtension(messages::gui_type));
    else
        info.gui_type = std::make_pair(false, (messages::GuiType)0);

    // collapsed
    if (options.HasExtension(messages::collapsed)) info.collapsed = options.GetExtension(messages::collapsed);

    // update signals
    if (options.HasExtension(messages::update_signals)) info.update_signals = options.GetExtension(messages::update_signals);

    // print_description
    if (options.HasExtension(messages::print_description)) info.print_description = options.GetExtension(messages::print_description);

    // message only options
    if (info.new_message)
    {
        const p::MessageOptions& msg_options = message->GetDescriptor()->options();

        // msg_description
        if (msg_options.HasExtension(messages::msg_description))
            info.msg_description = std::make_pair(true, msg_options.GetExtension(messages::msg_description));
        else
            info.msg_description = std::make_pair(false, "");
    }
}

std::string MessageParser::nestedNamesToString(const std::list<std::string>& namespaces)
{
    std::string text;
    if (namespaces.empty()) return text;
    auto it = namespaces.begin();
    text += *it;
    for (it = std::next(it); it != namespaces.end(); ++it) text += "." + *it;
    return text;
}

}  // namespace corbo

// MESSAGE_SUPPORT
#endif
