#include "veff_lcao.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

namespace hamilt
{

template class Veff<OperatorLCAO<double, double>>;

template class Veff<OperatorLCAO<std::complex<double>, double>>;

template class Veff<OperatorLCAO<std::complex<double>, std::complex<double>>>;

template<typename TK, typename TR>
Veff<OperatorLCAO<TK, TR>>::~Veff()
{
    //do nothing
}

template<typename TK, typename TR>
void Veff<OperatorLCAO<TK, TR>>::contributeHR()
{
    ModuleBase::TITLE("Veff", "contributeHR");
    return;
}

template<typename TK, typename TR>
void Veff<OperatorLCAO<TK, TR>>::contributeHk(int ik)
{
    ModuleBase::TITLE("Veff", "contributeHk");
    ModuleBase::timer::tick("Veff", "contributeHk");
    if(std::is_same<TK, double>::value)
    {
        //-----------------------------------------
        //(1) prepare data for this k point.
        // copy the local potential from array.
        //-----------------------------------------
        const double* vr_eff1 = this->pot->get_effective_v(GlobalV::CURRENT_SPIN);
        const double* vofk_eff1 = this->pot->get_effective_vofk(GlobalV::CURRENT_SPIN);

        //--------------------------------------------
        // (3) folding matrix,
        // and diagonalize the H matrix (T+Vl+Vnl).
        //--------------------------------------------

        if(XC_Functional::get_func_type()==3 || XC_Functional::get_func_type()==5)
        {
            Gint_inout inout(vr_eff1, vofk_eff1, Gint_Tools::job_type::vlocal_meta);
            this->GG->cal_vlocal(&inout, this->LM, this->new_e_iteration);
        }
        else
        {
            Gint_inout inout(vr_eff1, Gint_Tools::job_type::vlocal);
            this->GG->cal_vlocal(&inout, this->LM, this->new_e_iteration);
        }

        this->new_e_iteration = false;
    }
    else if(std::is_same<TK, std::complex<double>>::value)
    {
        //-----------------------------------------
        //(1) prepare data for this k point.
        // copy the local potential from array.
        //-----------------------------------------
        double* vr_eff1 = this->pot->get_effective_v(GlobalV::CURRENT_SPIN);
        double* vofk_eff1 = this->pot->get_effective_vofk(GlobalV::CURRENT_SPIN);

        //--------------------------------------------
        //(2) check if we need to calculate
        // pvpR = < phi0 | v(spin) | phiR> for a new spin.
        //--------------------------------------------
        if (GlobalV::CURRENT_SPIN == this->GK->get_spin())
        {
            // GlobalV::ofs_running << " Same spin, same vlocal integration." << std::endl;
        }
        else
        {
            // GlobalV::ofs_running << " (spin change)" << std::endl;
            this->GK->reset_spin(GlobalV::CURRENT_SPIN);

            // if you change the place of the following code,
            // rememeber to delete the #include
            if(XC_Functional::get_func_type()==3 || XC_Functional::get_func_type()==5)
            {
                Gint_inout inout(vr_eff1, vofk_eff1, 0, Gint_Tools::job_type::vlocal_meta);
                this->GK->cal_gint(&inout);
            }
            else
            {
                // vlocal = Vh[rho] + Vxc[rho] + Vl(pseudo)
                Gint_inout inout(vr_eff1, 0, Gint_Tools::job_type::vlocal);
                this->GK->cal_gint(&inout);
            }

            // added by zhengdy-soc, for non-collinear case
            // integral 4 times, is there any method to simplify?
            if (GlobalV::NSPIN == 4)
            {
                for (int is = 1; is < 4; is++)
                {
                    vr_eff1 = this->pot->get_effective_v(is);
                    if(XC_Functional::get_func_type()==3 || XC_Functional::get_func_type()==5)
                    {
                        vofk_eff1 = this->pot->get_effective_vofk(is);
                    }
                    
                    if(XC_Functional::get_func_type()==3 || XC_Functional::get_func_type()==5)
                    {
                        Gint_inout inout(vr_eff1, vofk_eff1, is, Gint_Tools::job_type::vlocal_meta);
                        this->GK->cal_gint(&inout);
                    }
                    else
                    {
                        Gint_inout inout(vr_eff1, is, Gint_Tools::job_type::vlocal);
                        this->GK->cal_gint(&inout);
                    }
                }
            }
        }

        this->GK->folding_vl_k(ik, this->LM, this->kvec_d);
    }
    else
    {
        throw std::runtime_error("Veff::contributeHk: type not supported.");
    }
    

    ModuleBase::timer::tick("Veff", "contributeHk");
}

}