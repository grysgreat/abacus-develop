remove_definitions(-D__DEEPKS)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
remove_definitions(-D__EXX)
remove_definitions(-DUSE_PAW)

AddTest(
  TARGET io_json_output_json
  SOURCES ../json.cpp
)