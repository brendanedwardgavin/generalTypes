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

emin=20.5
emax=30.5
#emid = 25.5
#r=5
m0=8
nc=6
eps=1e-14
maxit=100
x0=rand(T,n,m0)
y0=rand(T,n,m0)

B=eye(T,n)


A=matrixOp(A)
B=matrixOp(B)
x0=generalVec(x0)

ncontours=2
results=Vector{Future}(ncontours)
emids=zeros(T,ncontours)
len=(emax-emin)/ncontours
r=len/2.0
for i in 1:ncontours
    emids[i]=emin+(i-1)*len+r
end

for i=1:ncontours
    results[i]=remotecall(feast_nsR,i+1,A,B,x0,nc,emids[i],r,eps,maxit)
end

for i in 1:ncontours
    wait(results[i])
end

#for i in 1:ncontours
#    println((results[i])[1])
#end
