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

emid = 25.5
r=5
m0=12
nc=6
eps=1e-14
maxit=100
x0=rand(T,n,m0)
y0=rand(T,n,m0)

B=eye(T,n)


#A=matrixOp(A)
options=Dict("n" => n, "nvctr" => m0, "nvctrp" => m0, "ptr" => A)

#A=matrixOp(A)
A=userOp(options)

#B=matrixOp(B)
B=userOp(Dict("n" => n, "nvctr" => m0, "nvctrp" => m0, "ptr" => B))
x0=generalVec(x0)

(l,x)=feast_nsR(A,B,x0,nc,emid,r,eps,maxit)

println(l)
