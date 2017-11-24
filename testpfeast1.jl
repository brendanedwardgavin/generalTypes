@everywhere include("bicgstab.jl")
@everywhere include("feast_ns.jl")
@everywhere include("gT2.jl")

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

emin=10.5
emax=110.5
#emid = 25.5
#r=5
m0=110
nc=6
eps=1e-14
maxit=100

y0=rand(T,n,m0)

B=eye(T,n)

ncontours=nprocs()-1
m0=convert(Int64,ceil(m0/ncontours))
x0=rand(T,n,m0)
results=Vector{Future}(ncontours)
emids=zeros(T,ncontours)
len=(emax-emin)/ncontours
r=len/2.0
for i in 1:ncontours
    emids[i]=emin+(i-1)*len+r
end

A=matrixOp(A)
B=matrixOp(B)
x0=generalVec(x0)

for i=1:ncontours
    #results[i]=remotecall(feast_nsR,i+1,A,B,x0,nc,emids[i],r,eps,maxit)
    results[i]=@spawnat i+1 feast_nsR(A,B,x0,nc,emids[i],r,eps,maxit) 
end

for i in 1:ncontours
    wait(results[i])
end

for i in 1:ncontours
    println(fetch(results[i])[1])
end
