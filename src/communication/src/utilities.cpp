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

#include <corbo-communication/utilities.h>
#include <corbo-core/console.h>

#include <vector>

#ifdef WIN32
// WinUser.h defines GetMessage which is also used by protobuf
#undef GetMessage
#endif

namespace corbo {

#ifdef MESSAGE_SUPPORT
bool util::get_oneof_field_type(const google::protobuf::Message& message, const std::string& oneof_name, std::string& item_name,
                                bool include_namespace)
{
    item_name.clear();
    const google::protobuf::OneofDescriptor* oneof = nullptr;
    if (oneof_name.empty())
    {
        // check if we have only a single oneof
        if (message.GetDescriptor()->oneof_decl_count() == 1)
        {
            oneof = message.GetDescriptor()->oneof_decl(0);
        }
        else
        {
            PRINT_ERROR_NAMED("Multiple oneofs found in the message but no oneof_name provided, I don't know which one to choose.");
            return false;
        }
    }
    else
    {
        oneof = message.GetDescriptor()->FindOneofByName(oneof_name);
    }
    if (oneof)
    {
        const google::protobuf::FieldDescriptor* oneof_field = message.GetReflection()->GetOneofFieldDescriptor(message, oneof);
        if (oneof_field && oneof_field->cpp_type() == google::protobuf::FieldDescriptor::CPPTYPE_MESSAGE)
        {
            // ok, we have found the field
            item_name = message.GetReflection()->GetMessage(message, oneof_field).GetTypeName();
            if (!include_namespace)
            {
                std::size_t found = item_name.find_last_of(".");
                item_name         = item_name.substr(found + 1);
            }
            return true;
        }
    }
    if (item_name.empty()) PRINT_ERROR_NAMED("no item found for oneof-name " + oneof_name + ".");
    return false;
}

bool util::get_oneof_field_type_expand_isolated(const google::protobuf::Message& top_message, const google::protobuf::OneofDescriptor* top_oneof,
                                                std::string& item_name, bool include_namespace, int max_depth)
{
    const google::protobuf::FieldDescriptor* top_oneof_field = top_message.GetReflection()->GetOneofFieldDescriptor(top_message, top_oneof);
    if (top_oneof_field && top_oneof_field->cpp_type() != google::protobuf::FieldDescriptor::CPPTYPE_MESSAGE) return false;

    const google::protobuf::Message& message = top_message.GetReflection()->GetMessage(top_message, top_oneof_field);

    int num_oneofs = message.GetDescriptor()->oneof_decl_count();

    if (num_oneofs == 1 && max_depth > 0)
    {
        get_oneof_field_type_expand_isolated(message, message.GetDescriptor()->oneof_decl(0), item_name, include_namespace, max_depth - 1);
    }
    else
    {
        // ok, we have found the field
        item_name = message.GetTypeName();
        if (!include_namespace)
        {
            std::size_t found = item_name.find_last_of(".");
            item_name         = item_name.substr(found + 1);
        }
    }
    return true;
}

bool util::get_oneof_field_type_expand_isolated(const google::protobuf::Message& message, const std::string& top_oneof_name, std::string& item_name,
                                                bool include_namespace, int max_depth)
{
    const google::protobuf::OneofDescriptor* top_oneof = nullptr;
    if (top_oneof_name.empty())
    {
        // check if we have only a single oneof
        if (message.GetDescriptor()->oneof_decl_count() == 1)
        {
            top_oneof = message.GetDescriptor()->oneof_decl(0);
        }
        else
        {
            PRINT_ERROR_NAMED("Multiple oneofs found in the message but no top_oneof_name provided, I don't know which one to choose.");
            return false;
        }
    }
    else
    {
        top_oneof = message.GetDescriptor()->FindOneofByName(top_oneof_name);
    }
    if (top_oneof)
    {
        return get_oneof_field_type_expand_isolated(message, top_oneof, item_name, include_namespace, max_depth);
    }
    return false;
}

#endif

}  // namespace corbo
