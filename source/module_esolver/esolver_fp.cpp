#include "esolver_fp.h"

#include "module_base/global_variable.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/cube_io.h"
#include "module_io/output_log.h"
#include "module_io/write_elecstat_pot.h"
#include "module_parameter/parameter.h"
namespace ModuleESolver
{

ESolver_FP::ESolver_FP()
{
    // pw_rho = new ModuleBase::PW_Basis();
    pw_rho = new ModulePW::PW_Basis_Big(GlobalV::device_flag, GlobalV::precision_flag);

    if (GlobalV::double_grid)
    {
        pw_rhod = new ModulePW::PW_Basis_Big(GlobalV::device_flag, GlobalV::precision_flag);
    }
    else
    {
        pw_rhod = pw_rho;
    }

    // temporary, it will be removed
    pw_big = static_cast<ModulePW::PW_Basis_Big*>(pw_rhod);
    pw_big->setbxyz(PARAM.inp.bx, PARAM.inp.by, PARAM.inp.bz);
    sf.set(pw_rhod, PARAM.inp.nbspline);

    GlobalC::ucell.symm.epsilon = GlobalC::ucell.symm.epsilon_input = PARAM.inp.symmetry_prec;
}

ESolver_FP::~ESolver_FP()
{
    delete pw_rho;
    if (GlobalV::double_grid)
    {
        delete pw_rhod;
    }
    delete this->pelec;
}

void ESolver_FP::before_all_runners(const Input_para& inp, UnitCell& cell)
{
    ModuleBase::TITLE("ESolver_FP", "before_all_runners");

    //! 1) read pseudopotentials
    if (!GlobalV::use_paw)
    {
        cell.read_pseudo(GlobalV::ofs_running);
    }

    //! 2) initialie the plane wave basis for rho
#ifdef __MPI
    this->pw_rho->initmpi(GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL, POOL_WORLD);
#endif
    if (this->classname == "ESolver_OF")
    {
        this->pw_rho->setfullpw(inp.of_full_pw, inp.of_full_pw_dim);
    }

    if (inp.nx * inp.ny * inp.nz == 0)
    {
        this->pw_rho->initgrids(inp.ref_cell_factor * cell.lat0, cell.latvec, 4.0 * inp.ecutwfc);
    }
    else
    {
        this->pw_rho->initgrids(inp.ref_cell_factor * cell.lat0, cell.latvec, inp.nx, inp.ny, inp.nz);
    }

    this->pw_rho->initparameters(false, 4.0 * inp.ecutwfc);
    this->pw_rho->ft.fft_mode = inp.fft_mode;
    this->pw_rho->setuptransform();
    this->pw_rho->collect_local_pw();
    this->pw_rho->collect_uniqgg();

    //! 3) initialize the double grid (for uspp) if necessary
    if (GlobalV::double_grid)
    {
        ModulePW::PW_Basis_Sup* pw_rhod_sup = static_cast<ModulePW::PW_Basis_Sup*>(pw_rhod);
#ifdef __MPI
        this->pw_rhod->initmpi(GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL, POOL_WORLD);
#endif
        if (this->classname == "ESolver_OF")
        {
            this->pw_rhod->setfullpw(inp.of_full_pw, inp.of_full_pw_dim);
        }
        if (inp.ndx * inp.ndy * inp.ndz == 0)
        {
            this->pw_rhod->initgrids(inp.ref_cell_factor * cell.lat0, cell.latvec, inp.ecutrho);
        }
        else
        {
            this->pw_rhod->initgrids(inp.ref_cell_factor * cell.lat0, cell.latvec, inp.ndx, inp.ndy, inp.ndz);
        }
        this->pw_rhod->initparameters(false, inp.ecutrho);
        this->pw_rhod->ft.fft_mode = inp.fft_mode;
        pw_rhod_sup->setuptransform(this->pw_rho);
        this->pw_rhod->collect_local_pw();
        this->pw_rhod->collect_uniqgg();
    }

    //! 4) print some information
    this->print_rhofft(inp, GlobalV::ofs_running);

    //! 5) initialize the charge extrapolation method if necessary
    this->CE.Init_CE(GlobalV::NSPIN, GlobalC::ucell.nat, GlobalC::ucell.omega, this->pw_rhod->nrxx, inp.chg_extrap);

    return;
}

//! Something to do after SCF iterations when SCF is converged or comes to the max iter step.
void ESolver_FP::after_scf(const int istep)
{
    // 0) output convergence information
    ModuleIO::output_convergence_after_scf(this->conv_elec, this->pelec->f_en.etot);

    // 1) write fermi energy
    ModuleIO::output_efermi(this->conv_elec, this->pelec->eferm.ef);

    if (istep % PARAM.inp.out_interval == 0)
    {
        // 3) write charge density
        if (PARAM.inp.out_chg[0])
        {
            for (int is = 0; is < GlobalV::NSPIN; is++)
            {
                double* data = nullptr;
                if (PARAM.inp.dm_to_rho)
                {
                    data = this->pelec->charge->rho[is];
                }
                else
                {
                    data = this->pelec->charge->rho_save[is];
                }
                std::string fn = GlobalV::global_out_dir + "/SPIN" + std::to_string(is + 1) + "_CHG.cube";
                ModuleIO::write_cube(
#ifdef __MPI
                    this->pw_big->bz,
                    this->pw_big->nbz,
                    this->pw_rhod->nplane,
                    this->pw_rhod->startz_current,
#endif
                    data,
                    is,
                    GlobalV::NSPIN,
                    istep,
                    fn,
                    this->pw_rhod->nx,
                    this->pw_rhod->ny,
                    this->pw_rhod->nz,
                    this->pelec->eferm.get_efval(is),
                    &(GlobalC::ucell),
                    PARAM.inp.out_chg[1],
                    1);
                if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
                {
                    fn = GlobalV::global_out_dir + "/SPIN" + std::to_string(is + 1) + "_TAU.cube";
                    ModuleIO::write_cube(
#ifdef __MPI
                        this->pw_big->bz,
                        this->pw_big->nbz,
                        this->pw_rhod->nplane,
                        this->pw_rhod->startz_current,
#endif
                        this->pelec->charge->kin_r_save[is],
                        is,
                        GlobalV::NSPIN,
                        istep,
                        fn,
                        this->pw_rhod->nx,
                        this->pw_rhod->ny,
                        this->pw_rhod->nz,
                        this->pelec->eferm.get_efval(is),
                        &(GlobalC::ucell));
                }
            }
        }

        // 4) write potential
        if (PARAM.inp.out_pot == 1 || PARAM.inp.out_pot == 3)
        {
            for (int is = 0; is < GlobalV::NSPIN; is++)
            {
                std::string fn = GlobalV::global_out_dir + "/SPIN" + std::to_string(is + 1) + "_POT.cube";

                ModuleIO::write_cube(
#ifdef __MPI
                    this->pw_big->bz,
                    this->pw_big->nbz,
                    this->pw_rhod->nplane,
                    this->pw_rhod->startz_current,
#endif
                    this->pelec->pot->get_effective_v(is),
                    is,
                    GlobalV::NSPIN,
                    istep,
                    fn,
                    this->pw_rhod->nx,
                    this->pw_rhod->ny,
                    this->pw_rhod->nz,
                    0.0, // efermi
                    &(GlobalC::ucell),
                    3,  // precision
                    0); // out_fermi
            }
        }
        else if (PARAM.inp.out_pot == 2)
        {
            std::string fn = GlobalV::global_out_dir + "/ElecStaticPot.cube";
            ModuleIO::write_elecstat_pot(
#ifdef __MPI
                this->pw_big->bz,
                this->pw_big->nbz,
#endif
                fn,
                istep,
                this->pw_rhod,
                this->pelec->charge,
                &(GlobalC::ucell),
                this->pelec->pot->get_fixed_v());
        }
    }
}

void ESolver_FP::init_after_vc(const Input_para& inp, UnitCell& cell)
{
    ModuleBase::TITLE("ESolver_FP", "init_after_vc");

    if (inp.mdp.md_prec_level == 2)
    {
        if (inp.nx * inp.ny * inp.nz == 0)
        {
            this->pw_rho->initgrids(cell.lat0, cell.latvec, 4.0 * inp.ecutwfc);
        }
        else
        {
            this->pw_rho->initgrids(cell.lat0, cell.latvec, inp.nx, inp.ny, inp.nz);
        }

        this->pw_rho->initparameters(false, 4.0 * inp.ecutwfc);
        this->pw_rho->setuptransform();
        this->pw_rho->collect_local_pw();
        this->pw_rho->collect_uniqgg();

        if (GlobalV::double_grid)
        {
            ModulePW::PW_Basis_Sup* pw_rhod_sup = static_cast<ModulePW::PW_Basis_Sup*>(pw_rhod);
            if (inp.ndx * inp.ndy * inp.ndz == 0)
            {
                this->pw_rhod->initgrids(cell.lat0, cell.latvec, inp.ecutrho);
            }
            else
            {
                this->pw_rhod->initgrids(cell.lat0, cell.latvec, inp.ndx, inp.ndy, inp.ndz);
            }
            this->pw_rhod->initparameters(false, inp.ecutrho);
            pw_rhod_sup->setuptransform(this->pw_rho);
            this->pw_rhod->collect_local_pw();
            this->pw_rhod->collect_uniqgg();
        }
    }
    else
    {
        // only G-vector and K-vector are changed due to the change of lattice
        // vector FFT grids do not change!!
        pw_rho->initgrids(cell.lat0, cell.latvec, pw_rho->nx, pw_rho->ny, pw_rho->nz);
        pw_rho->collect_local_pw();
        pw_rho->collect_uniqgg();

        if (GlobalV::double_grid)
        {
            this->pw_rhod->initgrids(cell.lat0, cell.latvec, pw_rhod->nx, pw_rhod->ny, pw_rhod->nz);
            this->pw_rhod->collect_local_pw();
            this->pw_rhod->collect_uniqgg();
        }

        GlobalC::ppcell.init_vloc(GlobalC::ppcell.vloc, pw_rhod);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "LOCAL POTENTIAL");
    }
    this->pelec->omega = GlobalC::ucell.omega;

    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        cell.symm.analy_sys(cell.lat, cell.st, cell.atoms, GlobalV::ofs_running);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SYMMETRY");
    }

    kv.set_after_vc(GlobalV::NSPIN, cell.G, cell.latvec);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT K-POINTS");

    return;
}

void ESolver_FP::print_rhofft(const Input_para& inp, std::ofstream& ofs)
{
    std::cout << " UNIFORM GRID DIM        : " << pw_rho->nx << " * " << pw_rho->ny << " * " << pw_rho->nz << std::endl;
    std::cout << " UNIFORM GRID DIM(BIG)   : " << pw_big->nbx << " * " << pw_big->nby << " * " << pw_big->nbz
              << std::endl;
    if (GlobalV::double_grid)
    {
        std::cout << " UNIFORM GRID DIM(DENSE) : " << pw_rhod->nx << " * " << pw_rhod->ny << " * " << pw_rhod->nz
                  << std::endl;
    }

    ofs << "\n\n\n\n";
    ofs << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
           ">>>>"
        << std::endl;
    ofs << " |                                                                 "
           "   |"
        << std::endl;
    ofs << " | Setup plane waves of charge/potential:                          "
           "   |"
        << std::endl;
    ofs << " | Use the energy cutoff and the lattice vectors to generate the   "
           "   |"
        << std::endl;
    ofs << " | dimensions of FFT grid. The number of FFT grid on each "
           "processor   |"
        << std::endl;
    ofs << " | is 'nrxx'. The number of plane wave basis in reciprocal space "
           "is   |"
        << std::endl;
    ofs << " | different for charege/potential and wave functions. We also set "
           "   |"
        << std::endl;
    ofs << " | the 'sticks' for the parallel of FFT. The number of plane waves "
           "   |"
        << std::endl;
    ofs << " | is 'npw' in each processor.                                     "
           "   |"
        << std::endl;
    ofs << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
           "<<<<"
        << std::endl;
    ofs << "\n\n\n\n";
    ofs << "\n SETUP THE PLANE WAVE BASIS" << std::endl;

    double ecut = 4 * PARAM.inp.ecutwfc;
    if (inp.nx * inp.ny * inp.nz > 0)
    {
        ecut = this->pw_rho->gridecut_lat * this->pw_rho->tpiba2;
        ofs << "use input fft dimensions for wave functions." << std::endl;
        ofs << "calculate energy cutoff from nx, ny, nz:" << std::endl;
    }

    ModuleBase::GlobalFunc::OUT(ofs, "energy cutoff for charge/potential (unit:Ry)", ecut);

    ModuleBase::GlobalFunc::OUT(ofs,
                                "fft grid for charge/potential",
                                this->pw_rho->nx,
                                this->pw_rho->ny,
                                this->pw_rho->nz);
    ModuleBase::GlobalFunc::OUT(ofs, "fft grid division", pw_big->bx, pw_big->by, pw_big->bz);
    ModuleBase::GlobalFunc::OUT(ofs, "big fft grid for charge/potential", pw_big->nbx, pw_big->nby, pw_big->nbz);
    ModuleBase::GlobalFunc::OUT(ofs, "nbxx", pw_big->nbxx);
    ModuleBase::GlobalFunc::OUT(ofs, "nrxx", this->pw_rho->nrxx);

    ofs << "\n SETUP PLANE WAVES FOR CHARGE/POTENTIAL" << std::endl;
    ModuleBase::GlobalFunc::OUT(ofs, "number of plane waves", this->pw_rho->npwtot);
    ModuleBase::GlobalFunc::OUT(ofs, "number of sticks", this->pw_rho->nstot);

    ofs << "\n PARALLEL PW FOR CHARGE/POTENTIAL" << std::endl;
    ofs << " " << std::setw(8) << "PROC" << std::setw(15) << "COLUMNS(POT)" << std::setw(15) << "PW" << std::endl;

    for (int i = 0; i < GlobalV::NPROC_IN_POOL; ++i)
    {
        ofs << " " << std::setw(8) << i + 1 << std::setw(15) << this->pw_rho->nst_per[i] << std::setw(15)
            << this->pw_rho->npw_per[i] << std::endl;
    }
    ofs << " --------------- sum -------------------" << std::endl;
    ofs << " " << std::setw(8) << GlobalV::NPROC_IN_POOL << std::setw(15) << this->pw_rho->nstot << std::setw(15)
        << this->pw_rho->npwtot << std::endl;

    ModuleBase::GlobalFunc::OUT(ofs, "number of |g|", this->pw_rho->ngg);
    ModuleBase::GlobalFunc::OUT(ofs, "max |g|", this->pw_rho->gg_uniq[this->pw_rho->ngg - 1]);
    ModuleBase::GlobalFunc::OUT(ofs, "min |g|", this->pw_rho->gg_uniq[0]);

    if (GlobalV::double_grid)
    {
        ofs << std::endl;
        ofs << std::endl;
        ofs << std::endl;
        double ecut = PARAM.inp.ecutrho;
        if (inp.ndx * inp.ndy * inp.ndz > 0)
        {
            ecut = this->pw_rhod->gridecut_lat * this->pw_rhod->tpiba2;
            ofs << "use input fft dimensions for the dense part of charge "
                   "density."
                << std::endl;
            ofs << "calculate energy cutoff from ndx, ndy, ndz:" << std::endl;
        }
        ModuleBase::GlobalFunc::OUT(ofs, "energy cutoff for dense charge/potential (unit:Ry)", ecut);

        ModuleBase::GlobalFunc::OUT(ofs,
                                    "fft grid for dense charge/potential",
                                    this->pw_rhod->nx,
                                    this->pw_rhod->ny,
                                    this->pw_rhod->nz);

        ModuleBase::GlobalFunc::OUT(ofs, "nrxx", this->pw_rhod->nrxx);

        ofs << "\n SETUP PLANE WAVES FOR dense CHARGE/POTENTIAL" << std::endl;
        ModuleBase::GlobalFunc::OUT(ofs, "number of plane waves", this->pw_rhod->npwtot);
        ModuleBase::GlobalFunc::OUT(ofs, "number of sticks", this->pw_rhod->nstot);

        ofs << "\n PARALLEL PW FOR dense CHARGE/POTENTIAL" << std::endl;
        ofs << " " << std::setw(8) << "PROC" << std::setw(15) << "COLUMNS(POT)" << std::setw(15) << "PW" << std::endl;

        for (int i = 0; i < GlobalV::NPROC_IN_POOL; ++i)
        {
            ofs << " " << std::setw(8) << i + 1 << std::setw(15) << this->pw_rhod->nst_per[i] << std::setw(15)
                << this->pw_rhod->npw_per[i] << std::endl;
        }
        ofs << " --------------- sum -------------------" << std::endl;
        ofs << " " << std::setw(8) << GlobalV::NPROC_IN_POOL << std::setw(15) << this->pw_rhod->nstot << std::setw(15)
            << this->pw_rhod->npwtot << std::endl;

        ModuleBase::GlobalFunc::OUT(ofs, "number of |g|", this->pw_rhod->ngg);
        ModuleBase::GlobalFunc::OUT(ofs, "max |g|", this->pw_rhod->gg_uniq[this->pw_rhod->ngg - 1]);
        ModuleBase::GlobalFunc::OUT(ofs, "min |g|", this->pw_rhod->gg_uniq[0]);
    }
}

} // namespace ModuleESolver
