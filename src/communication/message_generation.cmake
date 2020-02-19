## control-box-rst messages and services generation
macro(generate_corbo_proto_msgs arg_proto_dir arg_generated_dir arg_msg_target arg_srv_target arg_rpc_target arg_rpc_support arg_custom_proto_include)

set(proto_gen_cpp_dir ${arg_generated_dir})
set(proto_gen_python_dir ${arg_generated_dir})

# Collect all proto files
set(proto_dir ${arg_proto_dir})
file(GLOB_RECURSE proto_msg_files "${proto_dir}/messages/*.proto")
file(GLOB_RECURSE proto_srv_files "${proto_dir}/services/*.proto" )

# Create lists of message files to be generated.
set(proto_msg_gen_cpp_headers "")
set(proto_msg_gen_cpp_src "")
set(proto_msg_gen_python_files "")
foreach(proto_file ${proto_msg_files})
    file(RELATIVE_PATH proto_name "${proto_dir}/messages" ${proto_file})
    # remove proto file ending
    string(REGEX REPLACE "\\.[^.]*$" "" proto_name ${proto_name})
    list(APPEND proto_msg_gen_cpp_headers ${proto_gen_cpp_dir}/messages/${proto_name}.pb.h)
    list(APPEND proto_msg_gen_cpp_src ${proto_gen_cpp_dir}/messages/${proto_name}.pb.cc)
    list(APPEND proto_msg_gen_python_files ${proto_gen_python_dir}/messages/${proto_name}_pb2.py)
endforeach(proto_file ${proto_msg_files})


# Create lists of service files to be generated.
set(proto_srv_gen_cpp_headers "")
set(proto_srv_gen_cpp_src "")
set(proto_grpc_gen_cpp_headers "")
set(proto_grpc_gen_cpp_src "")
set(proto_srv_gen_python_files "")
foreach(proto_file ${proto_srv_files})
    file(RELATIVE_PATH proto_name "${proto_dir}/services" ${proto_file})
    # remove proto file ending
    string(REGEX REPLACE "\\.[^.]*$" "" proto_name ${proto_name})
    list(APPEND proto_srv_gen_cpp_headers ${proto_gen_cpp_dir}/services/${proto_name}.pb.h)
    list(APPEND proto_srv_gen_cpp_src ${proto_gen_cpp_dir}/services/${proto_name}.pb.cc)
    list(APPEND proto_grpc_gen_cpp_headers ${proto_gen_cpp_dir}/services/${proto_name}.grpc.pb.h)
    list(APPEND proto_grpc_gen_cpp_src ${proto_gen_cpp_dir}/services/${proto_name}.grpc.pb.cc)
    list(APPEND proto_srv_gen_python_files ${proto_gen_python_dir}/services/${proto_name}_pb2.py)
endforeach(proto_file ${proto_srv_files})

# get executable from protobuf package
get_target_property(PROTOBUF_PROTOC_EXECUTABLE libprotobuf PROTOBUF_PROTOC_EXECUTABLE)

# import official protobuf include directory (for extending those files)
get_target_property(PROTOBUF_INCLUDE_DIR libprotobuf INTERFACE_INCLUDE_DIRECTORIES)

# Check for custom protobuf message include dirs
# TODO allow lists with multiple include dirs
if (EXISTS ${arg_custom_proto_include})
     message(STATUS "Generating custom protobuf messages in ${arg_custom_proto_include}")
     file(GLOB_RECURSE custom_proto_msg_files "${arg_custom_proto_include}/*.proto")

     foreach(proto_file ${custom_proto_msg_files})
         file(RELATIVE_PATH proto_name "${arg_custom_proto_include}/messages/custom" ${proto_file})
         # remove proto file ending
         string(REGEX REPLACE "\\.[^.]*$" "" proto_name ${proto_name})
         list(APPEND proto_msg_gen_cpp_headers ${proto_gen_cpp_dir}/messages/custom/${proto_name}.pb.h)
         list(APPEND proto_msg_gen_cpp_src ${proto_gen_cpp_dir}/messages/custom/${proto_name}.pb.cc)
         list(APPEND proto_msg_gen_python_files ${proto_gen_python_dir}/messages/custom/${proto_name}_pb2.py)
         # append  ${custom_proto_msg_files} with files, since we need the list in the custom_command
         list(APPEND proto_msg_files ${proto_file})
     endforeach(proto_file ${custom_proto_msg_files})
     # also set custom include dir for protoc executable
     set(CUSTOM_PROTOC_INCLUDE_DIR -I=${arg_custom_proto_include})
endif (EXISTS ${arg_custom_proto_include})

# Generate message files
# Create custom command to invoke protoc and generate proto headers and sources
# The forloop is a workaround, since protoc seems to get problems if the number of proto file inputs is large.
# Refer to the services custom_command below to see how it was working before...
list(LENGTH proto_msg_files count)
math(EXPR count "${count}-1")
#foreach(protoitem ${proto_msg_files})
foreach(i RANGE ${count})
	list(GET proto_msg_files ${i} item_proto)
	list(GET proto_msg_gen_cpp_headers ${i} item_h)
	list(GET proto_msg_gen_cpp_src ${i} item_cpp)
	list(GET proto_msg_gen_python_files ${i} item_py)
	add_custom_command(
		COMMAND ${PROTOBUF_PROTOC_EXECUTABLE} -I=${proto_dir} -I=${PROTOBUF_INCLUDE_DIR} ${CUSTOM_PROTOC_INCLUDE_DIR} --cpp_out=${proto_gen_cpp_dir} --python_out=${proto_gen_python_dir} ${item_proto}
		WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
		DEPENDS  ${item_proto} libprotobuf
		OUTPUT ${item_h} ${item_cpp} ${item_py}
	)
endforeach(i)
add_custom_target(${arg_msg_target} ALL
    DEPENDS ${proto_msg_gen_cpp_headers} ${proto_msg_gen_cpp_src} ${proto_msg_gen_python_files} libprotobuf
)

if (${arg_rpc_support})
  # Generate service files

  # create custom command to invoke protoc and generate proto headers and sources
  add_custom_command(
      COMMAND ${PROTOBUF_PROTOC_EXECUTABLE} -I=${proto_dir} -I=${PROTOBUF_INCLUDE_DIR} ${CUSTOM_PROTOC_INCLUDE_DIR} --cpp_out=${proto_gen_cpp_dir} --python_out=${proto_gen_python_dir} ${proto_srv_files}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      DEPENDS  ${proto_srv_files} libprotobuf
      OUTPUT ${proto_srv_gen_cpp_headers} ${proto_srv_gen_cpp_src} ${proto_srv_gen_python_files}
  )
  add_custom_target(${arg_srv_target} ALL
      DEPENDS ${proto_srv_gen_cpp_headers} ${proto_srv_gen_cpp_src} ${proto_srv_gen_python_files} libprotobuf
  )

  # Generate GRPC files
  add_custom_command(
      COMMAND ${PROTOBUF_PROTOC_EXECUTABLE} -I=${proto_dir} -I=${PROTOBUF_INCLUDE_DIR} ${CUSTOM_PROTOC_INCLUDE_DIR} --grpc_out=${proto_gen_cpp_dir} --plugin=protoc-gen-grpc=${GRPC_CPP_PLUGIN_EXECUTABLE} ${proto_srv_files}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      DEPENDS  ${proto_srv_files} libprotobuf libgrpc
      OUTPUT ${proto_grpc_gen_cpp_headers} ${proto_grpc_gen_cpp_src}
  )
  add_custom_target(${arg_rpc_target} ALL
      DEPENDS ${proto_grpc_gen_cpp_headers} ${proto_grpc_gen_cpp_src} libprotobuf libgrpc
  )
else (${arg_rpc_support})
# clear variables
set(proto_srv_gen_cpp_src "")
set(proto_grpc_gen_cpp_src "")
endif (${arg_rpc_support})

endmacro()
