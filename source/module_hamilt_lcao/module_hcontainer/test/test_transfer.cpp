#include "gtest/gtest.h"
#include "../transfer.h"
#include "../hcontainer.h"
#include <chrono>
#ifdef __MPI
#include <mpi.h>
#include "../hcontainer_funcs.h"
#endif

// test_size is the number of atoms in the unitcell
// modify test_size to test different size of unitcell
int test_size = 2;
int test_nw = 10;

class TransferTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
#ifdef __MPI
        // MPI parallel settings
        MPI_Comm_size(MPI_COMM_WORLD, &dsize);
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

        // set up a unitcell, with one element and test_size atoms, each atom has test_nw orbitals
        ucell.ntype = 1;
        ucell.nat = test_size;
        ucell.atoms = new Atom[ucell.ntype];
        ucell.iat2it = new int[ucell.nat];
        ucell.iat2ia = new int[ucell.nat];
        ucell.itia2iat.create(ucell.ntype, ucell.nat);
        for (int iat = 0; iat < ucell.nat; iat++)
        {
            ucell.iat2it[iat] = 0;
            ucell.iat2ia[iat] = iat;
            ucell.itia2iat(0, iat) = iat;
        }
        ucell.atoms[0].na = test_size;
        ucell.atoms[0].nw = test_nw;
        ucell.set_iat2iwt(1);
        init_parav();
        // set up a HContainer with ucell
        HR_para = new hamilt::HContainer<double>(ucell, paraV);
    }

    void TearDown() override
    {
        delete HR_para;
        delete paraV;
        delete[] ucell.atoms;
        delete[] ucell.iat2it;
        delete[] ucell.iat2ia;
    }

#ifdef __MPI
    void init_parav()
    {
        int global_row = test_size * test_nw;
        int global_col = test_size * test_nw;
        std::ofstream ofs_running;
        paraV = new Parallel_Orbitals();
        paraV->set_block_size(2/* nb_2d set to be 2*/);
        paraV->set_proc_dim(dsize, 0);
        paraV->mpi_create_cart(MPI_COMM_WORLD);
        paraV->set_local2global(global_row, global_col, ofs_running, ofs_running);
        int lr = paraV->get_row_size();
        int lc = paraV->get_col_size();
        paraV->set_desc(global_row, global_col, lr);
        paraV->set_global2local(global_row, global_col, true, ofs_running);
        paraV->set_atomic_trace(ucell.get_iat2iwt(), test_size, global_row);
    }
#else
    void init_parav()
    {}
#endif

    UnitCell ucell;
    hamilt::HContainer<double>* HR_para = nullptr;
    Parallel_Orbitals *paraV;

    int dsize;
    int my_rank = 0;
};

TEST_F(TransferTest, serialToPara)
{
// get rank of process
#ifdef __MPI

    hamilt::HContainer<double>* HR_serial = nullptr;

// initialize HR_serial
    // if the master process, calculate the value of HR_serial and send to other processes
    if (my_rank == 0)
    {
        HR_serial = new hamilt::HContainer<double>(ucell);
        for(int i = 0; i < HR_serial->size_atom_pairs(); i++)
        {
            hamilt::AtomPair<double>& atom_pair = HR_serial->get_atom_pair(i);
            int atom_i = atom_pair.get_atom_i();
            int atom_j = atom_pair.get_atom_j();
            //lambda function to calculate value of array: (atom_i*test_size+atom_j+k)*test_nw + l
            auto value = [&](int k, int l) -> double {return (((atom_i*test_nw+k)*test_size+atom_j)*test_nw + l);};
            double* data = atom_pair.get_pointer(0);
            for(int k = 0; k < test_nw; k++)
            {
                for(int l = 0; l < test_nw; l++)
                {
                    *data = value(k, l);
                    ++data;
                }
            }
        }
    }

    /*
    hamilt::HTransSerial<double> trans_s(dsize, HR_serial);
    hamilt::HTransPara<double> trans_p(dsize, HR_para);
    // plan indexes
    for(int i=0; i<dsize;++i)
    {
        if(i != 0 && my_rank != 0) continue;
        if(i != my_rank)
        {
            if(my_rank == 0) 
            {
                //trans_s.cal_ap_indexes(i);
                trans_s.send_ap_indexes(i);
            }
            else
            {
                trans_p.receive_ap_indexes(0);
                //trans_p.cal_orb_indexes(0);
                trans_p.send_orb_indexes(0);
            }
            if(my_rank == 0) 
            {
                trans_s.receive_orb_indexes(i);
            }
        }
        else
        {
            std::vector<int> tmp_indexes;
            trans_s.cal_ap_indexes(i, &tmp_indexes);
            trans_p.receive_ap_indexes(i, &tmp_indexes);
            trans_p.cal_orb_indexes(i, &tmp_indexes);
            trans_s.receive_orb_indexes(i, &tmp_indexes);
        }
    }

    // send data
    for(int i=0; i<dsize;++i)
    {
        if(i != 0 && my_rank != 0) continue;
        if(i != my_rank)
        {
            if(my_rank == 0) 
            {
                trans_s.send_data(i);
            }
            else
            {
                trans_p.receive_data(0);
            }
        }
        else
        {
            std::vector<double> tmp_values;
            trans_s.pack_data(i, &tmp_values);
            trans_p.receive_data(i, &tmp_values);
        }
    }*/
    hamilt::transferSerial2Parallel(*HR_serial, HR_para, 0);

    // check data in HR_para
    for(int i = 0; i < HR_para->size_atom_pairs(); i++)
    {
        hamilt::AtomPair<double>& atom_pair = HR_para->get_atom_pair(i);
        int atom_i = atom_pair.get_atom_i();
        int atom_j = atom_pair.get_atom_j();
        //lambda function to calculate value of array: (atom_i*test_size+atom_j+k)*test_nw + l
        auto value = [&](int k, int l) -> double {return (((atom_i*test_nw+k)*test_size+atom_j)*test_nw + l);};
        double* data = atom_pair.get_pointer(0);
        auto row_indexes = paraV->get_indexes_row(atom_i);
        auto col_indexes = paraV->get_indexes_col(atom_j);
        for(int k = 0; k < row_indexes.size(); k++)
        {
            for(int l = 0; l < col_indexes.size(); l++)
            {
                EXPECT_NEAR(*data, value(row_indexes[k], col_indexes[l]), 1e-10);
                ++data;
            }
        }
    }
    if(my_rank == 0)
    {
        delete HR_serial;
    }
#endif
}

TEST_F(TransferTest, paraToSerial)
{
    // get rank of process
#ifdef __MPI

    hamilt::HContainer<double>* HR_serial = nullptr;

    // initialize HR_para
    for(int i = 0; i < HR_para->size_atom_pairs(); i++)
    {
        hamilt::AtomPair<double>& atom_pair = HR_para->get_atom_pair(i);
        int atom_i = atom_pair.get_atom_i();
        int atom_j = atom_pair.get_atom_j();
        //lambda function to calculate value of array: (atom_i*test_size+atom_j+k)*test_nw + l
        auto value = [&](int k, int l) -> double {return (((atom_i*test_nw+k)*test_size+atom_j)*test_nw + l);};
        double* data = atom_pair.get_pointer(0);
        auto row_indexes = paraV->get_indexes_row(atom_i);
        auto col_indexes = paraV->get_indexes_col(atom_j);
        for(int k = 0; k < row_indexes.size(); k++)
        {
            for(int l = 0; l < col_indexes.size(); l++)
            {
                *data = value(row_indexes[k], col_indexes[l]);
                ++data;
            }
        }
    }

    // initialize HR_serial
    // if the master process, calculate the value of HR_serial and send to other processes
    if (my_rank == 0)
    {
        HR_serial = new hamilt::HContainer<double>(ucell);
    }

    /*hamilt::HTransSerial<double> trans_s(dsize, HR_serial);
    hamilt::HTransPara<double> trans_p(dsize, HR_para);

    // plan indexes
    for(int i=0; i<dsize;++i)
    {
        if(i != 0 && my_rank != 0) continue;
        if(i != my_rank)
        {
            if(my_rank == 0) 
            {
                trans_s.send_ap_indexes(i);
            }
            else
            {
                trans_p.receive_ap_indexes(0);
                trans_p.send_orb_indexes(0);
            }
            if(my_rank == 0) 
            {
                trans_s.receive_orb_indexes(i);
            }
        }
        else
        {
            std::vector<int> tmp_indexes; 
            trans_s.cal_ap_indexes(i, &tmp_indexes);
            trans_p.receive_ap_indexes(i, &tmp_indexes);
            trans_p.cal_orb_indexes(i, &tmp_indexes);
            trans_s.receive_orb_indexes(i, &tmp_indexes);
        }
    }

    // send data
    for(int i=0; i<dsize;++i)
    {
        if(i != 0 && my_rank != 0) continue;
        if(i != my_rank)
        {
            if(my_rank != 0) 
            {
                trans_p.send_data(0);
            }
            else
            {
                trans_s.receive_data(i);
            }
        }
        else
        {
            std::vector<double> tmp_values;
            trans_p.pack_data(i, &tmp_values);
            trans_s.receive_data(i, &tmp_values);
        }
    }*/
    hamilt::transferParallel2Serial(*HR_para, HR_serial, 0);

    // check data in HR_serial
    if(my_rank == 0)
    {
        for(int i = 0; i < HR_serial->size_atom_pairs(); i++)
        {
            hamilt::AtomPair<double>& atom_pair = HR_serial->get_atom_pair(i);
            int atom_i = atom_pair.get_atom_i();
            int atom_j = atom_pair.get_atom_j();
            //lambda function to calculate value of array: (atom_i*test_size+atom_j+k)*test_nw + l
            auto value = [&](int k, int l) -> double {return (((atom_i*test_nw+k)*test_size+atom_j)*test_nw + l);};
            double* data = atom_pair.get_pointer(0);
            for(int k = 0; k < test_nw; k++)
            {
                for(int l = 0; l < test_nw; l++)
                {
                    EXPECT_NEAR(*data , value(k, l), 1e-10);
                    ++data;
                }
            }
        }
    }
    delete HR_serial;
#endif
}

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}