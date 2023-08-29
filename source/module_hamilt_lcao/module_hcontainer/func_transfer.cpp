#include "./hcontainer_funcs.h"
#include "./transfer.h"

#ifdef __MPI
#include <mpi.h>

//#include <chrono>

namespace hamilt
{
// transfer the HContainer from serial object to parallel object
template <typename TR>
void transferSerial2Parallel(const hamilt::HContainer<TR>& hR_s, hamilt::HContainer<TR>* hR_p, const int serial_rank)
{
    int my_rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    hamilt::HTransSerial<TR>* trans_s = nullptr;
    if (my_rank == serial_rank)
    {
        trans_s = new hamilt::HTransSerial<TR>(size, const_cast<hamilt::HContainer<TR>*>(&hR_s));
    }
    hamilt::HTransPara<TR> trans_p(size, hR_p);
    // plan indexes
    //std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
    if (my_rank == serial_rank)
    {
        // send indexes to other ranks
        for (int i = 0; i < size; ++i)
        {
            if (i == serial_rank)
                continue;
            trans_s->send_ap_indexes(i);
        }

        { // transfer serial_rank to serial_rank
            std::vector<int> tmp_indexes;
            trans_s->cal_ap_indexes(serial_rank, &tmp_indexes);
            trans_p.receive_ap_indexes(serial_rank, &tmp_indexes);
            trans_p.cal_orb_indexes(serial_rank, &tmp_indexes);
            trans_s->receive_orb_indexes(serial_rank, &tmp_indexes);
        }
        
        // receive indexes from other ranks
        for (int i = 0; i < size; ++i)
        {
            if (i == serial_rank)
                continue;
            trans_s->receive_orb_indexes(i);
        }
    }
    else // my_rank != serial_rank
    {
        // receive ap_indexes from serial_rank and then send orb_indexes to serial_rank
        trans_p.receive_ap_indexes(serial_rank);
        trans_p.send_orb_indexes(serial_rank);
    }
    std::vector<TR> all_values;
    long max_size = 0;
    if (my_rank == serial_rank)
    {
        // calculate max size of values
        max_size = trans_s->get_max_size();
        all_values.resize(max_size * size);
    }
    MPI_Bcast(&max_size, 1, MPI_LONG, serial_rank, MPI_COMM_WORLD);
    std::vector<TR> receive_values(max_size);
    /*std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time0
        = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    start_time = std::chrono::high_resolution_clock::now();*/

    // send data
    if (my_rank == serial_rank)
    {
        for (int i = 0; i < size; ++i)
        {
            trans_s->pack_data(i, (all_values.data() + i * max_size));
        }
    }

    /*end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> pre_scatter_time
        = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    start_time = std::chrono::high_resolution_clock::now();*/

    // MPI_scatter to send values
    MPI_Scatter(all_values.data(),
                max_size,
                MPITraits<TR>::datatype(),
                receive_values.data(),
                max_size,
                MPITraits<TR>::datatype(),
                serial_rank,
                MPI_COMM_WORLD);
    /*end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> scatter_time
        = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    start_time = std::chrono::high_resolution_clock::now();*/

    // receive data
    trans_p.receive_data(serial_rank, receive_values.data());

    /*end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> post_scatter_time
        = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    std::cout << " S2P: my_rank = " << my_rank << " indexes_time = " << elapsed_time0.count()
              << " data_trans_time = " << pre_scatter_time.count()<<" "<<scatter_time.count()
              <<" "<<post_scatter_time.count() << std::endl;*/

    if (my_rank == serial_rank)
    {
        delete trans_s;
    }
}

// transfer the HContainer from parallel object to serial object
template <typename TR>
void transferParallel2Serial(const hamilt::HContainer<TR>& hR_p, hamilt::HContainer<TR>* hR_s, const int serial_rank)
{
    int my_rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    hamilt::HTransSerial<TR>* trans_s = nullptr;
    if (my_rank == serial_rank)
    {
        trans_s = new hamilt::HTransSerial<TR>(size, hR_s);
    }
    hamilt::HTransPara<TR> trans_p(size, const_cast<hamilt::HContainer<TR>*>(&hR_p));

    // plan indexes
    //std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
    if (my_rank == serial_rank)
    {
        { // transfer serial_rank to serial_rank
            std::vector<int> tmp_indexes;
            trans_s->cal_ap_indexes(serial_rank, &tmp_indexes);
            trans_p.receive_ap_indexes(serial_rank, &tmp_indexes);
            trans_p.cal_orb_indexes(serial_rank, &tmp_indexes);
            trans_s->receive_orb_indexes(serial_rank, &tmp_indexes);
        }

        // send indexes to other ranks
        for (int i = 0; i < size; ++i)
        {
            if (i == serial_rank)
                continue;
            trans_s->send_ap_indexes(i);
        }
        // receive indexes from other ranks
        for (int i = 0; i < size; ++i)
        {
            if (i == serial_rank)
                continue;
            trans_s->receive_orb_indexes(i);
        }
    }
    else // my_rank != serial_rank
    {
        trans_p.receive_ap_indexes(serial_rank);
        trans_p.send_orb_indexes(serial_rank);
    }
    /*std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time0
        = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    start_time = std::chrono::high_resolution_clock::now();*/

    // send data
    std::vector<TR> receive_values;
    long max_size;
    if (my_rank == serial_rank)
    {
        max_size = trans_s->get_max_size();
        receive_values.resize(max_size * size);
    }
    MPI_Bcast(&max_size, 1, MPI_LONG, serial_rank, MPI_COMM_WORLD);

    // MPI_gather to receive values
    std::vector<TR> send_values;
    send_values.resize(max_size);
    trans_p.pack_data(serial_rank, send_values.data());
    /*end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> pre_gather_time
        = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    start_time = std::chrono::high_resolution_clock::now();*/
    MPI_Gather(send_values.data(),
               max_size,
               MPITraits<TR>::datatype(),
               receive_values.data(),
               max_size,
               MPITraits<TR>::datatype(),
               serial_rank,
               MPI_COMM_WORLD);
    /*end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> gather_time
        = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    start_time = std::chrono::high_resolution_clock::now();*/
    if (my_rank == serial_rank)
    {
        for (int i = 0; i < size; ++i)
        {
            trans_s->receive_data(i, &receive_values[i * max_size]);
        }
    }
    
    /*end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> post_gather_time
        = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    std::cout << " P2S: my_rank = " << my_rank << " elapsed_time0 = " << elapsed_time0.count()
              << " data_trans_time = " << pre_gather_time.count()<<" "<<gather_time.count()
              <<" "<<post_gather_time.count() << std::endl;*/

    if (my_rank == serial_rank)
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