include("jmpi.jl")

ierr=MPI_Init()

include("bicgstab.jl")
include("feast_ns.jl")
include("gT2.jl")

srand(231)

T=Complex128

#generate matrix:
n=500
A=rand(T,n,n)
(Q,R)=qr(A)
for i in 1:n
    A[:,i]=A[:,i]/norm(A[:,i])
end
Q=copy(A)
L=diagm(1.0*collect(1:n))
#L=diagm(rand(n))
A=Q*L*inv(Q)
#A=rand(n,n)
B=eye(T,n)




#FEAST parameters:
ncontour=1
emin=9.5
emax=50.5
m0t=50
m0=convert(Int64,ceil(m0t/ncontour))
nc=4
eps=1e-14
maxit=100
x0=rand(T,n,m0)
y0=rand(T,n,m0)




#Custom Type wrappers:
#A=matrixOp(A)
options=Dict("n" => n, "nvctr" => m0, "nvctrp" => m0, "ptr" => A)

A=matrixOp(A)
#A=userOp(options)
B=matrixOp(B)
#B=userOp(Dict("n" => n, "nvctr" => m0, "nvctrp" => m0, "ptr" => B))
x0=generalVec(x0)




#Parallelization:
MPI_COMM_WORLD=mpiworld()
nprocs=MPI_Comm_size(MPI_COMM_WORLD)
myrank=MPI_Comm_rank(MPI_COMM_WORLD)

color=mod(myrank,ncontour)
key=1
contourcomm=MPI_Comm_split(MPI_COMM_WORLD,color,key)

r=(emax-emin)/(2.0*ncontour)
delta=(emax-emin)/ncontour

emid=emin+color*delta+r




#Run FEAST:
(l,x)=feast_nsR_MPI(contourcomm,color,A,B,x0,nc,emid,r,eps,maxit)

println(l)

ierr=MPI_Finalize()
