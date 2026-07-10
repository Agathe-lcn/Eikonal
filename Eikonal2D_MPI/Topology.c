#include <mpi.h>

int main(int argc, char **argv){

    MPI_Init(&argc, &argv);

    int ndims = 2;
    int mesh_size[2] = {200,200};
    int nproc_per_dim[2] = {0,0};
    int periods_per_dim[2] = {0,0};
    int rang, nproc;

    MPI_Comm_rank(MPI_COMM_WORLD, &(rang));
    MPI_Comm_size(MPI_COMM_WORLD, &(nproc));

    int err = 0;

    // Création de la topologie cartésienne
    err = MPI_Dims_create(nproc, ndims, nproc_per_dim);
    MPI_Comm GRID_COMM;
    err = MPI_Cart_create(MPI_COMM_WORLD, ndims, nproc_per_dim, periods_per_dim, 0, &GRID_COMM);

    // Récupération des coordonnées du rang dans la topologie
    int proc_coords[2];
}