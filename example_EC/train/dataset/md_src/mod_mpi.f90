module mod_mpi
    include 'mpif.h'
    integer inode,inode_tot    ! inode=1
    integer nnodes,nnodes_tot  ! nnodes=1
    integer MPI_COMM_MOL   ! within each molecule, actually size 1 here
    integer status(MPI_STATUS_SIZE)
    integer natom_n   ! divided number of atom
end module mod_mpi

