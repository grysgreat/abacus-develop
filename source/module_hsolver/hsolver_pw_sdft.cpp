#include "hsolver_pw_sdft.h"

#include "module_base/global_function.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_elecstate/module_charge/symmetry_rho.h"

#include <algorithm>

namespace hsolver {
void HSolverPW_SDFT::solve(hamilt::Hamilt<std::complex<double>>* pHamilt,
                           psi::Psi<std::complex<double>>& psi,
                           elecstate::ElecState* pes,
                           ModulePW::PW_Basis_K* wfc_basis,
                           Stochastic_WF& stowf,
                           const int istep,
                           const int iter,
                           const std::string method_in,

                           const int scf_iter_in,
                           const bool need_subspace_in,
                           const int diag_iter_max_in,
                           const double iter_diag_thr_in,

                           const bool skip_charge)
{
    ModuleBase::TITLE("HSolverPW_SDFT", "solve");
    ModuleBase::timer::tick("HSolverPW_SDFT", "solve");

    const int npwx = psi.get_nbasis();
    const int nbands = psi.get_nbands();
    const int nks = psi.get_nk();

    this->scf_iter = scf_iter_in;
    this->need_subspace = need_subspace_in;
    this->diag_iter_max = diag_iter_max_in;
    this->iter_diag_thr = iter_diag_thr_in;

    // prepare for the precondition of diagonalization
    std::vector<double> precondition(psi.get_nbasis(), 0.0);

    // select the method of diagonalization
    this->method = method_in;
    // report if the specified diagonalization method is not supported
    const std::initializer_list<std::string> _methods
        = {"cg", "dav", "dav_subspace", "bpcg"};
    if (std::find(std::begin(_methods), std::end(_methods), this->method)
        == std::end(_methods)) {
        ModuleBase::WARNING_QUIT("HSolverPW::solve",
                                 "This method of DiagH is not supported!");
    }

    // part of KSDFT to get KS orbitals
    for (int ik = 0; ik < nks; ++ik) {
        pHamilt->updateHk(ik);
        if (nbands > 0 && GlobalV::MY_STOGROUP == 0) {
            this->updatePsiK(pHamilt, psi, ik);
            // template add precondition calculating here
            update_precondition(precondition, ik, this->wfc_basis->npwk[ik]);
            /// solve eigenvector and eigenvalue for H(k)
            double* p_eigenvalues = &(pes->ekb(ik, 0));
            this->hamiltSolvePsiK(pHamilt, psi, precondition, p_eigenvalues);
        }

        stoiter.stohchi.current_ik = ik;

#ifdef __MPI
        if (nbands > 0) {
            MPI_Bcast(&psi(ik, 0, 0),
                      npwx * nbands,
                      MPI_DOUBLE_COMPLEX,
                      0,
                      PARAPW_WORLD);
            MPI_Bcast(&(pes->ekb(ik, 0)), nbands, MPI_DOUBLE, 0, PARAPW_WORLD);
        }
#endif
        stoiter.orthog(ik, psi, stowf);
        stoiter.checkemm(ik, istep, iter, stowf); // check and reset emax & emin
    }

    this->output_iterInfo();

    // psi only should be initialed once for PW
    if (!this->initialed_psi)
    {
        this->initialed_psi = true;
    }

    for (int ik = 0; ik < nks; ik++) {
        // init k
        if (nks > 1) {
            pHamilt->updateHk(ik);
}
        stoiter.stohchi.current_ik = ik;
        stoiter.calPn(ik, stowf);
    }

    stoiter.itermu(iter, pes);
    stoiter.calHsqrtchi(stowf);
    if (skip_charge) {
        ModuleBase::timer::tick("HSolverPW_SDFT", "solve");
        return;
    }
    //(5) calculate new charge density
    // calculate KS rho.
    if (nbands > 0) {
        pes->psiToRho(psi);
#ifdef __MPI
        MPI_Bcast(&pes->f_en.eband, 1, MPI_DOUBLE, 0, PARAPW_WORLD);
#endif
    } else {
        for (int is = 0; is < GlobalV::NSPIN; is++) {
            ModuleBase::GlobalFunc::ZEROS(pes->charge->rho[is],
                                          pes->charge->nrxx);
        }
    }
    // calculate stochastic rho
    stoiter.sum_stoband(stowf, pes, pHamilt, wfc_basis);

    // will do rho symmetry and energy calculation in esolver
    ModuleBase::timer::tick("HSolverPW_SDFT", "solve");
    return;
}

double HSolverPW_SDFT::set_diagethr(double diag_ethr_in,
                                    const int istep,
                                    const int iter,
                                    const double drho) {
    if (iter == 1) {
        if (istep == 0) {
            if (GlobalV::init_chg == "file") {
                diag_ethr_in = 1.0e-5;
            }
            diag_ethr_in = std::max(diag_ethr_in, GlobalV::PW_DIAG_THR);
        } else {
            diag_ethr_in = std::max(diag_ethr_in, 1.0e-5);
}
    } else {
        if (GlobalV::NBANDS > 0 && this->stoiter.KS_ne > 1e-6) {
            diag_ethr_in
                = std::min(diag_ethr_in,
                           0.1 * drho / std::max(1.0, this->stoiter.KS_ne));
        } else {
            diag_ethr_in = 0.0;
}
    }

    return diag_ethr_in;
}
} // namespace hsolver