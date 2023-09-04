#ifndef VEFFLCAO_H
#define VEFFLCAO_H
#include "module_base/timer.h"
#include "module_elecstate/potentials/potential_new.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.h"
#include "module_hamilt_lcao/module_gint/gint_gamma.h"
#include "module_hamilt_lcao/module_gint/gint_k.h"
#include "operator_lcao.h"

namespace hamilt
{

#ifndef __VEFFTEMPLATE
#define __VEFFTEMPLATE

template <class T>
class Veff : public T
{
};

#endif

template <typename TK, typename TR>
class Veff<OperatorLCAO<TK, TR>> : public OperatorLCAO<TK, TR>
{
  public:
    Veff<OperatorLCAO<TK, TR>>(Gint_k* GK_in,
                          Local_Orbital_Charge* loc_in,
                          LCAO_Matrix* LM_in,
                          const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                          elecstate::Potential* pot_in,
                          hamilt::HContainer<TR>* hR_in,
                          std::vector<TK>* hK_in)
        : GK(GK_in),
          loc(loc_in),
          pot(pot_in),
          OperatorLCAO<TK, TR>(LM_in, kvec_d_in, hR_in, hK_in)
    {
        this->cal_type = lcao_gint;
    }
    Veff<OperatorLCAO<TK, TR>>(Gint_Gamma* GG_in,
                          Local_Orbital_Charge* loc_in,
                          LCAO_Matrix* LM_in,
                          const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                          elecstate::Potential* pot_in,
                          hamilt::HContainer<TR>* hR_in,
                          std::vector<TK>* hK_in
                          )
        : GG(GG_in), loc(loc_in), pot(pot_in),
        OperatorLCAO<TK, TR>(LM_in, kvec_d_in, hR_in, hK_in)
    {
        this->cal_type = lcao_gint;
    }

    ~Veff<OperatorLCAO<TK, TR>>();

    virtual void contributeHR() override;

    virtual void contributeHk(int ik) override;

  private:
    // used for k-dependent grid integration.
    Gint_k* GK = nullptr;

    // used for gamma only algorithms.
    Gint_Gamma* GG = nullptr;

    // Charge calculating method in LCAO base and contained grid base calculation: DM_R, DM, pvpR_reduced
    Local_Orbital_Charge* loc = nullptr;

    elecstate::Potential* pot = nullptr;
};

} // namespace hamilt
#endif