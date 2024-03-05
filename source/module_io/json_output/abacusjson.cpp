#include "abacusjson.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
namespace Json
{

#ifdef __RAPIDJSON
rapidjson::Document AbacusJson::doc;

bool isNum(std::string str)  
{  
	std::stringstream sin;  
    sin<<str;
	double d;  
	char c;  
	if(!(sin >> d))  
		return false;
	
	if (sin >> c) 
		return false;
	return true;  
}


void AbacusJson::add_nested_member(std::vector<std::string>::iterator begin,
                                   std::vector<std::string>::iterator end,
                                   rapidjson::Value& val,
                                   rapidjson::Value& parent,
                                   rapidjson::Document::AllocatorType& allocator,
                                   bool IsArray
                                   )
{
    if (begin != end)
    {
        rapidjson::Value key((*begin).c_str(), allocator);
        std::string key_str = *begin;



        if (begin + 1 == end)
        {

            if(isNum(key_str)){
                std::string::size_type sz;
                int index = std::stoi(key_str,&sz);
                
                if(index!=0){
                    parent[index-1] = val;
                }
                else {
                    int arr_size = parent.Size();
                    parent[arr_size-1] = val;
                }
            }
            // if key exists, then overwrite it
            else if (parent.HasMember(key))
            {
                if(parent[key].IsArray()){
                    parent[key].PushBack(val, allocator);
                }else{
                    // if key is an object, then warn the user
                    if (parent[key].IsObject())
                    {
                        std::cout << "Warning: write to json, key " << *begin
                                << " exist and is an object, and abacus will overwrite it with a value." << std::endl;
                    }
                    parent[key] = val;
                }
            }
            else{
                if(IsArray==true){
                    rapidjson::Value arr(rapidjson::kArrayType);
                    arr.PushBack(val, allocator);
                    parent.AddMember(key, arr, allocator);
                } else{
                    parent.AddMember(key, val, allocator);
                }
                
            }
        }
        else
        {
            if(isNum(key_str)){
                std::string::size_type sz;
                int index = std::stoi(key_str,&sz);
                
                if(index!=0){
                    add_nested_member(begin + 1, end, val, parent[index-1], allocator,IsArray);
                }
                else {
                    int arr_size = parent.Size();
                    add_nested_member(begin + 1, end, val, parent[arr_size-1], allocator,IsArray);
                }
            }
            // need to check if the key exists
            else if (parent.HasMember(key))
            {
                // this key should be an object
                if (!parent[key].IsObject())
                {
                    std::cout << "Warning: write to json, key " << *begin
                              << " exist and is not an object, and abacus will add it as a middle node." << std::endl;
                }
                add_nested_member(begin + 1, end, val, parent[key], allocator,IsArray);
            }
            else
            {
                rapidjson::Value paraent_val(rapidjson::kObjectType);
                add_nested_member(begin + 1, end, val, paraent_val, allocator,IsArray);
                parent.AddMember(key, paraent_val, allocator);
            }
        }
    }
}
// Output the json to a file
void AbacusJson::write_to_json(std::string filename)
{
    rapidjson::StringBuffer buffer;
    rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(buffer);
    doc.Accept(writer);

    std::ofstream ofs(filename);
    ofs << buffer.GetString();
    ofs.close();
};
  template <>
  void AbacusJson::add_json(std::vector<std::string> keys, const std::string& value,bool IsArray)
  {
      if (!doc.IsObject())
      {
          doc.SetObject();
      }
      rapidjson::Value val(value.c_str(), doc.GetAllocator());
      add_nested_member(keys.begin(), keys.end(), val, doc, doc.GetAllocator(),IsArray);
  }


// Overloaded template functions for json class objects
  template <>
  void AbacusJson::add_json(std::vector<std::string> keys, const rapidjson::Value& value,bool IsArray)
  {

        if (!doc.IsObject())
        {
            doc.SetObject();
        }

        rapidjson::Value val(value,doc.GetAllocator());
        add_nested_member(keys.begin(), keys.end(), val, doc, doc.GetAllocator(),IsArray);
  }
#endif
} // namespace Json
