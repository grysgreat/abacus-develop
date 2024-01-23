list(APPEND objects_json
    abacusjson.cpp
    json_general_info.cpp

)

add_library(
    json_output
    OBJECT
    ${objects_json}
)

if(ENABLE_COVERAGE)
  add_coverage(json_output)
endif()

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
