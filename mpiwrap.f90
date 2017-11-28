function mpiworld()
    use mpi
    integer :: mpiworld
    mpiworld=MPI_COMM_WORLD
end function mpiworld

subroutine mpirank(rank)
    use mpi
    implicit none
    integer :: rank
    integer :: ierr
    call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr) 
end subroutine mpirank

function commSplit(comm,color,key)
    use mpi
    implicit none
    integer :: commSplit
    integer :: comm
    integer :: color,key,ierr
    integer :: newcomm
    call MPI_Comm_split(comm,color,key,newcomm,ierr)
    commSplit=newcomm
end function commSplit


subroutine MPI_zAllReduceSum(sendbuf,recvbuf,lecount,comm)
    use mpi
    implicit none
    complex*16, dimension(*) :: sendbuf,recvbuf
    integer :: lecount,ierr,comm

    call MPI_AllReduce(sendbuf,recvbuf,lecount,MPI_DOUBLE_COMPLEX, MPI_SUM,comm,ierr)
    
end

!subroutine mpiBarrier(comm)
!    use mpi
!    implicit none
!   integer :: comm
!    call MPI_Barrier(comm)
!end

function getMPIsize(comm)
    use mpi
    implicit none
    integer :: getMPIsize
    integer :: comm
    integer :: ierr

    call MPI_Comm_size(comm, getMPIsize,ierr)

end function getMPIsize

function getrank(comm)
    use mpi
    implicit none
    integer :: getrank
    integer :: comm
    integer :: ierr
    call MPI_Comm_rank(comm,getrank,ierr)
end function getrank

subroutine printRank(comm)
    use mpi
    implicit none
    integer :: comm
    integer :: ierr
    integer :: rank

    call MPI_Comm_rank(comm,rank,ierr)
    print *,rank
end subroutine printRank

