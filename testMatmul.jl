include("gT2.jl")

srand(231)

T=Complex128

#generate matrix:
n=500
m0=1
for i in 1:10
A=rand(T,n,n)
B=rand(T,n,n)
z=rand(T)

x=rand(T,n,m0)

Am=matrixOp(A)
Bm=matrixOp(B)
xm=generalVec(x)

b=zeros(T,n,m0)

ccall((:zmymatmul_,"./mymatmul.so"),Void,(Ptr{Cuchar},Ref{Complex128},Ref{Complex128},Ref{Complex128},Ref{Int32},Ref{Int32},Ref{Int32}),Ref{Cuchar}('N'),A,x,b,Ref{Int32}(n),Ref{Int32}(n),Ref{Int32}(m0))

println("Error=",norm((A)*x-((Am)*xm).vec))
println("       ",norm(A*x-b))
println("      ",norm((B-A)*x-((Bm-Am)*xm).vec))
println("      ",norm((z*B-A)*x-((z*Bm-Am)*xm).vec))
end
