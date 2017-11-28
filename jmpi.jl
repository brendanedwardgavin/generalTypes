function MPI_Init()
    ierr=0
    #ccall((:MPI_Init,"/local/gavin/OpenMPI/lib/libmpi.so"),Void,(Ref{Int32},),Ref{Int32}(ierr))
    ccall((:mpi_init_,"/local/gavin/OpenMPI/lib/libmpi_mpifh.so"),Void,(Ref{Int32},),Ref{Int32}(ierr))

    return ierr
end

function mpiworld()
    return ccall((:mpiworld_,"./mpiwrap.so"),Int32,())
end

#function MPI_Comm_rank(commworld)
#    myrank=convert(Cint,0)
#    ierr=0
#    ccall((:mpi_comm_rank_,"/local/gavin/OpenMPI/lib/libmpi_mpifh.so"),Void,(Ref{Cint},Ref{Cint},Ref{Cint}),Ref{Cint}(commworld),Ref{Cint}(myrank),Ref{Cint}(ierr))
#    return myrank
#end

function MPI_Comm_rank(comm)
#    rank=convert(Cint,0)
#    ccall((:mpirank_,"./mpiwrap.so"),Void,(Ref{Cint},),Ref{Cint}(rank))
    #println(rank)
#    return rank
    return ccall((:getrank_,"./mpiwrap.so"),Int32,(Ref{Cint},),Ref{Cint}(comm))
end

function MPI_Comm_size(comm)
#    rank=convert(Cint,0)
#    ccall((:mpirank_,"./mpiwrap.so"),Void,(Ref{Cint},),Ref{Cint}(rank))
    #println(rank)
#    return rank
    return ccall((:getmpisize_,"./mpiwrap.so"),Int32,(Ref{Cint},),Ref{Cint}(comm))
end

function MPI_Comm_split(comm,color,key)
    #ccall((:mpi_comm_split,"/local/gavin/OpenMPI/lib/libmpi_mpifh.so"),Void,(Ref{Cint},Ref{Cint},Ref{Cint},Ref{Cint}),Ref{Cint}(comm),Ref{Cint}(color),Ref{Cint}(key),Ref{Cint}(newcomm))
    newcomm=ccall((:commsplit_,"./mpiwrap.so"),Cint,(Ref{Cint},Ref{Cint},Ref{Cint}),Ref{Cint}(comm),Ref{Cint}(color),Ref{Cint}(key))
    return newcomm
end

function MPI_Finalize()
    ierr=0
    ccall((:MPI_Finalize,"/local/gavin/OpenMPI/lib/libmpi.so"),Void,(Ref{Int32},),Ref{Int32}(ierr))
    return ierr
end
