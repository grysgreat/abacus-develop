#include "./transfer.h"
#include "./hcontainer_funcs.h"

#ifdef __MPI
#include <mpi.h>

namespace hamilt
{
// transfer the HContainer from serial object to parallel object
template<typename TR>
void transferSerial2Parallel(const hamilt::HContainer<TR>& hR_s,
                             hamilt::HContainer<TR>* hR_p,
                             const int serial_rank)
{
    int my_rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    hamilt::HTransSerial<TR>* trans_s = nullptr;
    if(my_rank == serial_rank)
    {
        trans_s = new hamilt::HTransSerial<TR>(size, const_cast<hamilt::HContainer<TR>*>(&hR_s));
    }
    hamilt::HTransPara<TR> trans_p(size, hR_p);
    // plan indexes
    if(my_rank == serial_rank)
    {
        {// transfer serial_rank to serial_rank 
            std::vector<int> tmp_indexes;
            trans_s->cal_ap_indexes(serial_rank, &tmp_indexes);
            trans_p.receive_ap_indexes(serial_rank, &tmp_indexes);
            trans_p.cal_orb_indexes(serial_rank, &tmp_indexes);
            trans_s->receive_orb_indexes(serial_rank, &tmp_indexes);
        }
        // send indexes to other ranks
        for(int i=0; i<size;++i)
        {
            if(i == serial_rank) continue;
            trans_s->send_ap_indexes(i);
        }
        // receive indexes from other ranks
        for(int i=0; i<size;++i)
        {
            if(i == serial_rank) continue;
            trans_s->receive_orb_indexes(i);
        }
    }
    else // my_rank != serial_rank
    {
        // receive ap_indexes from serial_rank and then send orb_indexes to serial_rank
        trans_p.receive_ap_indexes(serial_rank);
        trans_p.send_orb_indexes(serial_rank);
    }

    // send data
    if(my_rank == serial_rank)
    {
        {// transfer serial_rank to serial_rank 
            std::vector<TR> tmp_values;
            trans_s->pack_data(serial_rank, &tmp_values);
            trans_p.receive_data(serial_rank, &tmp_values);
        }
        // send data to other ranks
        for(int i=0; i<size;++i)
        {
            if(i == serial_rank) continue;
            trans_s->send_data(i);
        }
    }
    else // my_rank != serial_rank
    {
        trans_p.receive_data(serial_rank);
    }
    
    if(my_rank == serial_rank)
    {
        delete trans_s;
    }
}

//transfer the HContainer from parallel object to serial object
template<typename TR>
void transferParallel2Serial(const hamilt::HContainer<TR>& hR_p,
                             hamilt::HContainer<TR>* hR_s,
                             const int serial_rank)
{
    int my_rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    hamilt::HTransSerial<TR>* trans_s = nullptr;
    if(my_rank == serial_rank)
    {
        trans_s = new hamilt::HTransSerial<TR>(size, hR_s);
    }
    hamilt::HTransPara<TR> trans_p(size, const_cast<hamilt::HContainer<TR>*>(&hR_p));

    // plan indexes
    if(my_rank == serial_rank)
    {
        {// transfer serial_rank to serial_rank 
            std::vector<int> tmp_indexes; 
            trans_s->cal_ap_indexes(serial_rank, &tmp_indexes);
            trans_p.receive_ap_indexes(serial_rank, &tmp_indexes);
            trans_p.cal_orb_indexes(serial_rank, &tmp_indexes);
            trans_s->receive_orb_indexes(serial_rank, &tmp_indexes);
        }
        // send indexes to other ranks
        for(int i=0; i<size;++i)
        {
            if(i == serial_rank) continue;
            trans_s->send_ap_indexes(i);
        }
        // receive indexes from other ranks
        for(int i=0; i<size;++i)
        {
            if(i == serial_rank) continue;
            trans_s->receive_orb_indexes(i);
        }
    }
    else // my_rank != serial_rank
    {
        trans_p.receive_ap_indexes(serial_rank);
        trans_p.send_orb_indexes(serial_rank);
    }

    // send data
    if(my_rank != serial_rank)
    {
        trans_p.send_data(serial_rank);
    }
    else // my_rank == serial_rank
    {
        {// transfer serial_rank itself 
            std::vector<TR> tmp_values;
            trans_p.pack_data(serial_rank, &tmp_values);
            trans_s->receive_data(serial_rank, &tmp_values);
        }
        // receive data from other ranks
        for(int i=0; i<size;++i)
        {
            if(i == serial_rank) continue;
            trans_s->receive_data(i);
        }
    }

    if(my_rank == serial_rank)
    {
        delete trans_s;
    }
}

// specialize for double and std::complex<double>
template void transferSerial2Parallel(const hamilt::HContainer<double>& hR_s,
                                      hamilt::HContainer<double>* hR_p,
                                      const int serial_rank);
template void transferSerial2Parallel(const hamilt::HContainer<std::complex<double>>& hR_s,
                                        hamilt::HContainer<std::complex<double>>* hR_p,
                                        const int serial_rank);
template void transferParallel2Serial(const hamilt::HContainer<double>& hR_p,
                                        hamilt::HContainer<double>* hR_s,
                                        const int serial_rank);
template void transferParallel2Serial(const hamilt::HContainer<std::complex<double>>& hR_p,
                                        hamilt::HContainer<std::complex<double>>* hR_s,
                                        const int serial_rank);

} // namespace hamilt   

#endif // __MPI