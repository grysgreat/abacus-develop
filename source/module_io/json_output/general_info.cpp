#include "general_info.h"
#include "../para_json.h"
#include "abacusjson.h"
#include "module_base/parallel_global.h"
#include "module_io/input.h"

namespace Json
{
void gen_general_info(Input& input)
{

#ifdef VERSION
    const std::string version = VERSION;
#else
    const std::string version = "unknown";
#endif
#ifdef COMMIT
    const std::string commit = COMMIT;
#else
    const std::string commit = "unknown";
#endif

    // start_time
    std::time_t start_time = input.get_start_time();
    std::string start_time_str;
    convert_time(start_time, start_time_str);

    // end_time
    std::time_t time_now = std::time(NULL);
    std::string end_time_str;
    convert_time(time_now, end_time_str);

    int mpi_num = Parallel_Global::mpi_number;
    int omp_num = Parallel_Global::omp_number;

    AbacusJson::add_json({"general_info", "version"}, version);
    AbacusJson::add_json({"general_info", "commit"}, commit);
    AbacusJson::add_json({"general_info", "device"}, input.device);
    AbacusJson::add_json({"general_info", "mpi_num"}, mpi_num);
    AbacusJson::add_json({"general_info", "omp_num"}, omp_num);
    AbacusJson::add_json({"general_info", "pseudo_dir"}, input.pseudo_dir);
    AbacusJson::add_json({"general_info", "orbital_dir"}, input.orbital_dir);
    AbacusJson::add_json({"general_info", "stru_file"}, input.stru_file);
    AbacusJson::add_json({"general_info", "kpt_file"}, input.kpoint_file);
    AbacusJson::add_json({"general_info", "start_time"}, start_time_str);
    AbacusJson::add_json({"general_info", "end_time"}, end_time_str);
}
} // namespace Json