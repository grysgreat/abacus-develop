#include "readin_info.h"
#include "module_io/input.h"
#include "../para_json.h"
#include "abacusjson.h"


//Add json objects to init
namespace Json
{

#ifdef __RAPIDJSON


void gen_stru(UnitCell *ucell){

    int ntype = ucell->ntype;

    //atom mag
    for(int i=0;i<ntype;i++){
        int na = ucell->atoms[i].na;

        Json::jsonValue magArray(JarrayType);
        for(int j=0;j<na;j++){
            magArray.JPushBack(ucell->atoms[i].mag[j]);
        }
        Json::AbacusJson::add_json({"readin","stru","mag"}, magArray,true);
    }
    
    
    //file name of pseudopotential file
    std::string* pseudo_fn = ucell->pseudo_fn;
    for(int i=0;i<ntype;i++){
        std::string str = pseudo_fn[i];
        Json::AbacusJson::add_json({"readin","stru","pp"}, str,true);
    }


    //file name of orbital file
    std::string* orbital_fn = ucell->orbital_fn;
    for(int i=0;i<ntype;i++){
        std::string str = GlobalV::global_orbital_dir + orbital_fn[i];
        if(!str.compare("")){
            str = "null";
        }
        Json::AbacusJson::add_json({"readin","stru","orb"}, str,true);
    }

    //atom label
    std::string* label = ucell->atom_label;
    for(int i=0;i<ntype;i++){
        std::string str = label[i];
        Json::AbacusJson::add_json({"readin","stru","label"}, str,true);
    }

    //atom element
    for(int i=0;i<ntype;i++){
        std::string str = ucell->atoms[i].label;
        Json::AbacusJson::add_json({"readin","stru","element"}, str,true);
    }
    

    //coordinate
    double lat0 = ucell->lat0;
    for(int i=0;i<ntype;i++){
        ModuleBase::Vector3<double>* tau = ucell->atoms[i].tau;
        int na = ucell->atoms[i].na;
        for(int j=0;j<na;j++){
            Json::jsonValue coordinateArray(JarrayType);
            coordinateArray.JPushBack(tau[j][0]*lat0);
            coordinateArray.JPushBack(tau[j][1]*lat0);
            coordinateArray.JPushBack(tau[j][2]*lat0);
            Json::AbacusJson::add_json({"readin","stru","coordinate"}, coordinateArray,true);
        }
    }


    //cell 
    {
        Json::jsonValue cellArray1(JarrayType);
        Json::jsonValue cellArray2(JarrayType);
        Json::jsonValue cellArray3(JarrayType);
        cellArray1.JPushBack(ucell->latvec.e11);
        cellArray1.JPushBack(ucell->latvec.e12);
        cellArray1.JPushBack(ucell->latvec.e13);
        cellArray2.JPushBack(ucell->latvec.e21);
        cellArray2.JPushBack(ucell->latvec.e22);
        cellArray2.JPushBack(ucell->latvec.e23);
        cellArray3.JPushBack(ucell->latvec.e31);
        cellArray3.JPushBack(ucell->latvec.e32);
        cellArray3.JPushBack(ucell->latvec.e33);
        Json::AbacusJson::add_json({"readin","stru","cell"}, cellArray1,true);
        Json::AbacusJson::add_json({"readin","stru","cell"}, cellArray2,true);
        Json::AbacusJson::add_json({"readin","stru","cell"}, cellArray3,true);
    }
    return;
}


#endif
} // namespace Json