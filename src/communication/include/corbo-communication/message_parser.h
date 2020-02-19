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

#ifndef SRC_COMMUNICATION_INCLUDE_CORBO_COMMUNICATION_MESSAGE_PARSER_H_
#define SRC_COMMUNICATION_INCLUDE_CORBO_COMMUNICATION_MESSAGE_PARSER_H_

#ifdef MESSAGE_SUPPORT

#include <corbo-communication/messages/descriptor_extensions.pb.h>
#include <google/protobuf/message.h>
#include <functional>
#include <list>
#include <string>
#include <utility>
#include <vector>

namespace corbo {

/**
 * @brief Parser for protobuf messages
 *
 * @ingroup communication
 *
 * This class provides methods for recursively iterating
 * a protobuf message and triggering user-defined callback
 * functions each time a field or submessage is expanded.
 * A special MessageParser::FieldInformation struct is filled
 * with some extra information that can be relevant
 * for post-processing purposes.
 *
 * Check the available member functions for supported callbacks.
 *
 * @author Christoph Rösmann (christoph.roesmann@tu-dortmund.de)
 */
class MessageParser
{
 public:
    //! Field information that is passed to each callback
    struct FieldInformation
    {
        std::string msg_type;                 //!< The message type including message namespaces
        std::string field_name;               //!< Name of the current field
        std::vector<std::string> enum_list;   //!< Contains available enum items if current field is part of an enum
        std::list<std::string> nested_names;  //!< Store field names of parent messages

        // support manipulation of the field from outside
        // (using reflection interface):
        const google::protobuf::FieldDescriptor* field_raw            = nullptr;  //!< Raw pointer to the current field descriptor
        google::protobuf::Message* message_raw                        = nullptr;  //!< Pointer to the current message (can be a submessage!)
        const google::protobuf::OneofDescriptor* oneof_descriptor     = nullptr;  //!< if part of oneof store pointer to oneof-descriptor
        const google::protobuf::FieldDescriptor* oneof_selected_field = nullptr;  //!< if part of oneof store here the selected field

        bool new_message = true;  //!< true, if the field is newly constructed (should work only for message types rather than for primitive types)

        // The following list contains extra options that are defined in 'messages/descriptor_extensions.proto'
        // better would be std::optional (C++17)
        std::pair<bool, std::string> label;            //!< A user defined label for the current field
        std::pair<bool, std::string> description;      //!< Description of the current field
        std::pair<bool, std::string> msg_description;  //!< Description of the current message (only in use with a message rather than a field)
        std::pair<bool, std::string> default_value;    //!< Desired default value
        std::pair<bool, double> min_value;             //!< Expected minimum value
        std::pair<bool, double> max_value;             //!< Expected maximum value
        bool dynamic_size = false;        //!< If type is an array / repeated field: is the size dynamic or fixed according to initialization
        std::pair<bool, int> min_size;    //!< If dynamic_size==1: expected minimum size
        std::pair<bool, int> max_size;    //!< If dynamic_size==1: expected maximum size
        std::pair<bool, bool> info_text;  //!< Flags the value of a string param as info text label

        std::pair<bool, messages::GuiType> gui_type;  //!< Desired element in a possible gui according to messages::GuiType (TEXTEDIT, SLIDER, ...)
        bool collapsed         = false;               //!< True if this field (submessage) is intended to be initiliazed collapsed in a gui
        bool update_signals    = false;               //!< True if a modification of this field might change signals for SignalTargetInterface
        bool print_description = false;               //!< If true, the desciption is meant to be printed in addition to the tooltip
    };

    //! Possible events that can occure in message callbacks (see setCallbackMessageEvent)
    enum class MessageEvent { MessageStart, MessageEnd, OneOf };

    template <typename T>
    using CbValue = std::function<void(const T&, const FieldInformation&)>;

    using CbMessageEvent = std::function<void(const MessageEvent&, const FieldInformation&)>;

    //! Callback for int32 fields with prototype void(int,FieldInformation)
    void setCallbackValueInt32(CbValue<std::int32_t> fun_value_int32) { _fun_value_int32 = fun_value_int32; }
    //! Callback for repeated int32 fields with prototype void(std::vector<int>,FieldInformation)
    void setCallbackValueInt32Array(CbValue<std::vector<std::int32_t>> fun_value_int32_array) { _fun_value_int32_array = fun_value_int32_array; }

    //! Callback for int32 fields with prototype void(unsigned int,FieldInformation)
    void setCallbackValueUInt32(CbValue<std::uint32_t> fun_value_uint32) { _fun_value_uint32 = fun_value_uint32; }
    //! Callback for repeated uint32 fields with prototype void(std::vector<unsigned int>,FieldInformation)
    void setCallbackValueUInt32Array(CbValue<std::vector<std::uint32_t>> fun_value_uint32_array) { _fun_value_uint32_array = fun_value_uint32_array; }

    //! Callback for int64 fields with prototype void(long int,FieldInformation)
    void setCallbackValueInt64(CbValue<std::int64_t> fun_value_int64) { _fun_value_int64 = fun_value_int64; }
    //! Callback for repeated int64 fields with prototype void(std::vector<long int>,FieldInformation)
    void setCallbackValueInt64Array(CbValue<std::vector<std::int64_t>> fun_value_int64_array) { _fun_value_int64_array = fun_value_int64_array; }

    //! Callback for uint64 fields with prototype void(long unsigned int,FieldInformation)
    void setCallbackValueUInt64(CbValue<std::uint64_t> fun_value_uint64) { _fun_value_uint64 = fun_value_uint64; }
    //! Callback for repeated uint64 fields with prototype void(std::vector<long unsigned int>,FieldInformation)
    void setCallbackValueUInt64Array(CbValue<std::vector<std::uint64_t>> fun_value_uint64_array) { _fun_value_uint64_array = fun_value_uint64_array; }

    //! Callback for float fields with prototype void(float,FieldInformation)
    void setCallbackValueFloat(CbValue<float> fun_value_float) { _fun_value_float = fun_value_float; }
    //! Callback for repeated double fields with prototype void(std::vector<double>,FieldInformation)
    void setCallbackValueFloatArray(CbValue<std::vector<float>> fun_value_float_array) { _fun_value_float_array = fun_value_float_array; }

    //! Callback for double fields with prototype void(double,FieldInformation)
    void setCallbackValueDouble(CbValue<double> fun_value_double) { _fun_value_double = fun_value_double; }
    //! Callback for repeated double fields with prototype void(std::vector<double>,FieldInformation)
    void setCallbackValueDoubleArray(CbValue<std::vector<double>> fun_value_double_array) { _fun_value_double_array = fun_value_double_array; }

    //! Callback for string fields with prototype void(std::string,FieldInformation)
    void setCallbackValueString(CbValue<std::string> fun_value_string) { _fun_value_string = fun_value_string; }
    //! Callback for repeated string fields with prototype void(std::vector<std::string>,FieldInformation)
    void setCallbackValueStringArray(CbValue<std::vector<std::string>> fun_value_string_array) { _fun_value_string_array = fun_value_string_array; }

    //! Callback for string info fields with prototype void(std::string,FieldInformation)
    void setCallbackValueStringInfo(CbValue<std::string> fun_value_string_info) { _fun_value_string_info = fun_value_string_info; }

    //! Callback for bool fields with prototype void(bool,FieldInformation)
    void setCallbackValueBool(CbValue<bool> fun_value_bool) { _fun_value_bool = fun_value_bool; }
    //! Callback for repeated bool fields with prototype void(std::vector<bool>,FieldInformation)
    void setCallbackValueBoolArray(CbValue<std::vector<bool>> fun_value_bool_array) { _fun_value_bool_array = fun_value_bool_array; }

    //! Callback for enumeration fields with prototype void(std::string,FieldInformation), the string contains the currently selected field name
    void setCallbackValueEnum(CbValue<std::string> fun_value_enum) { _fun_value_enum = fun_value_enum; }
    //! Callback for repeated enumeration fields with prototype void(std::vector<std::string>,FieldInformation)
    void setCallbackValueEnumArray(CbValue<std::vector<std::string>> fun_value_enum_array) { _fun_value_enum_array = fun_value_enum_array; }

    //! Callback for a new (sub-)message event according to MessageEvent
    void setCallbackMessageEvent(CbMessageEvent fun_msg_event) { _fun_msg_event = fun_msg_event; }

    /**
     * @brief Parse a protobuf (sub-) message and invoke relevant callbacks
     * @param[in]  message               Message to be parsed (writable in order to allow setting parameters from outside)
     * @param[in]  expand_oneof          Expand one-of fields while parsing
     * @param[in]  skip_inactive_oneof   Skip uninitialized one-of fields (no one-of selected)
     * @param[out] nested_names          Temporary cache for parent messages names
     *                                   since they are not available from submessages (visitor list for messages) [optional]
     */
    void parse(google::protobuf::Message* message, bool expand_oneof = true, bool skip_inactive_oneof = true,
               std::list<std::string>* nested_names = nullptr);

    /**
     * @brief Parse a protobuf field and invoke relevant callbacks
     * @param[in]  message               Current (parent-)message (writable in order to allow setting parameters from outside)
     * @param[in]  field                 Field (descriptor) to be parsed
     * @param[in]  expand_oneof          Expand one-of fields while parsing
     * @param[in]  expand_first_oneof    In case \c expand_oneof == false: expand only the first oneof:
     *                                   this is necessary if \c field is already part of a oneof
     * @param[in]  skip_inactive_oneof   Skip uninitialized one-of fields (no one-of selected)
     * @param[out] nested_names          Temporary cache for parent messages names
     *                                   since they are not available from submessages (visitor list for messages) [optional]
     */
    void parseField(google::protobuf::Message* message, const google::protobuf::FieldDescriptor* field, bool expand_oneof = true,
                    bool expand_first_oneof = true, bool skip_inactive_oneof = true, std::list<std::string>* nested_names = nullptr);

    //! Activate verbosity: add console information in case a message field is found for which no callback is registered
    void setVerbose(bool active) { _verbose = active; }

    //! Convert nested_names / namespaces to a single string, e.g. my/name/space/field
    static std::string nestedNamesToString(const std::list<std::string>& namespaces);

 protected:
    //! Expand field in case it is a new submessage
    void consumeFieldMessage(google::protobuf::Message* message, const google::protobuf::Reflection* reflection,
                             const google::protobuf::FieldDescriptor* field, bool expand_oneof = true, bool skip_inactive_oneof = true,
                             std::list<std::string>* nested_names = nullptr);
    //! Expand field in case it is a single or repeated value
    void consumeFieldValue(google::protobuf::Message* message, const google::protobuf::Reflection* reflection,
                           const google::protobuf::FieldDescriptor* field, std::list<std::string>* nested_names = nullptr);

    //! Collect extra field information which is passed to the callbacks
    void getFieldInformation(google::protobuf::Message* message, const google::protobuf::FieldDescriptor* field, FieldInformation& info,
                             std::list<std::string>* nested_names = nullptr, bool is_message = false);

    CbValue<std::int32_t> _fun_value_int32;
    CbValue<std::vector<std::int32_t>> _fun_value_int32_array;
    CbValue<std::uint32_t> _fun_value_uint32;
    CbValue<std::vector<std::uint32_t>> _fun_value_uint32_array;
    CbValue<std::int64_t> _fun_value_int64;
    CbValue<std::vector<std::int64_t>> _fun_value_int64_array;
    CbValue<std::uint64_t> _fun_value_uint64;
    CbValue<std::vector<std::uint64_t>> _fun_value_uint64_array;
    CbValue<float> _fun_value_float;
    CbValue<std::vector<float>> _fun_value_float_array;
    CbValue<double> _fun_value_double;
    CbValue<std::vector<double>> _fun_value_double_array;
    CbValue<std::string> _fun_value_string;
    CbValue<std::vector<std::string>> _fun_value_string_array;
    CbValue<std::string> _fun_value_string_info;
    CbValue<bool> _fun_value_bool;
    CbValue<std::vector<bool>> _fun_value_bool_array;
    CbValue<std::string> _fun_value_enum;
    CbValue<std::vector<std::string>> _fun_value_enum_array;

    CbMessageEvent _fun_msg_event;

    bool _verbose = true;
};

}  // namespace corbo

#endif  // MESSAGE_SUPPORT

#endif  // SRC_COMMUNICATION_INCLUDE_CORBO_COMMUNICATION_MESSAGE_PARSER_H_
