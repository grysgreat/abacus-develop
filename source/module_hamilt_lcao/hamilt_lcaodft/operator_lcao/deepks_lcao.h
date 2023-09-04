#ifndef DEEPKSLCAO_H
#define DEEPKSLCAO_H
#include "module_base/timer.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.h"
#include "operator_lcao.h"

namespace hamilt
{

#ifndef __DEEPKSTEMPLATE
#define __DEEPKSTEMPLATE

template <class T>
class DeePKS : public T
{
};

#endif

template <typename TK, typename TR>
class DeePKS<OperatorLCAO<TK, TR>> : public OperatorLCAO<TK, TR>
{
  public:
    DeePKS<OperatorLCAO<TK, TR>>(Local_Orbital_Charge* loc_in,
                            LCAO_Matrix* LM_in,
                            const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                            HContainer<TR>* hR_in,
                            std::vector<TK>* hK_in,
                            const int& nks_in)
        : loc(loc_in),
          nks(nks_in),
          OperatorLCAO<TK, TR>(LM_in, kvec_d_in, hR_in, hK_in)
    {
        this->cal_type = lcao_deepks;
    }

    virtual void contributeHR() override;

    virtual void contributeHk(int ik) override;

  private:
    Local_Orbital_Charge* loc;

    bool HR_fixed_done = false;

    const int& nks;
};

} // namespace hamilt
#endif