mpicxx -DENABLE_MPI -O3 -Wall -std=c++11   -c -o mpi_ensemble_tcf_rpmd.o ensemble_tcf_rpmd.cpp -Itools -Ipotentials
mpicxx -o mpi_ensemble_tcf_rpmd mpi_ensemble_tcf_rpmd.o libtools.a -lfftw3 -larmadillo
rm -f mpi_ensemble_tcf_rpmd.o
