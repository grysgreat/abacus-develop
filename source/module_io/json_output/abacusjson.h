#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#ifdef __RAPIDJSON
#include "rapidjson/document.h"
#include "rapidjson/prettywriter.h"
#include "rapidjson/stringbuffer.h"

namespace Json
{

using jsonValue = rapidjson::Value;
// This class is used to construct a json value, and output some key values to a json file.
class AbacusJson
{
  public:
    // Add a key value pair to the json
    template <typename T>
    static void add_json(std::vector<std::string> keys, const T& value,bool IsArray)
    {
        if (!doc.IsObject())
        {
            doc.SetObject();
        }
        rapidjson::Value val(value);
        add_nested_member(keys.begin(), keys.end(), val, doc, doc.GetAllocator(),IsArray);
    }

    // Output the json to a file
    static void write_to_json(std::string filename);

    /**
    *   Recursive solution, using template derivation mechanism, 
    *   take out the first one at a time, shorten the size of the parameter package.
    */
    
    jsonValue create_jsonObj(){
      jsonValue object(rapidjson::kObjectType);
      return object;
    }



  private:
    static rapidjson::Document doc;

    static void add_nested_member(std::vector<std::string>::iterator begin,
                                  std::vector<std::string>::iterator end,
                                  rapidjson::Value& val,
                                  rapidjson::Value& parent,
                                  rapidjson::Document::AllocatorType& allocator,
                                  bool IsArray = false
                                  );
};
template <>
void AbacusJson::add_json(std::vector<std::string> keys, const std::string& value,bool IsArray);

template <>
void AbacusJson::add_json(std::vector<std::string> keys, const rapidjson::Value& value,bool IsArray);


} // namespace Json
#endif