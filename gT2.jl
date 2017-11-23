import Base: +, -, *,/,transpose,print,norm,copy,conj,Ac_mul_B,zero,size,getindex,setindex!,ctranspose,zeros

#=
TODO:
---Allocate new memory when doing Y=X[:,j]?? Or when doing arithmetic operations and memory copies? Probably best to do it in getindex; need to change this
=#

abstract type generalOp end

struct compositeOp <: generalOp #subtype of generalOp for building composite linear maps, e.g. A=B*C+z*D
    #carries two functions:  opmult(x)=A*x and adjopmult(x)=A'*x
    opmult::Function
    adjopmult::Function
    isAdjoint::Bool
end

struct matrixOp <: generalOp #subtype of generalOp that stores a matrix to use for its functions
    A::Array{Complex128,2}
    opmult::Function
    adjopmult::Function
    isAdjoint::Bool
end

mutable struct metadata_handle
   placehold::NTuple{8,Int64}
end

struct userOp <: generalOp
  ptr::Array{Complex128,2}
  handle:: metadata_handle
  opmult::Function
  adjopmult::Function
  isAdjoint::Bool
end

struct userVec
  ptr::Array{Complex128,2}
  handle:: metadata_handle
  isCovector::Bool
end


#functions for accessing the same fields from both compositeOp and matrixOp in the same format
#necessary because abstract types in Julia do not have member fields
adjOpMult(A::userOp)=A.adjopmult
opMult(A::userOp)=A.opmult
isOpAdjoint(A::userOp)=A.isAdjoint


adjOpMult(A::compositeOp)=A.adjopmult
opMult(A::compositeOp)=A.opmult
adjOpMult(A::matrixOp)=A.adjopmult
opMult(A::matrixOp)=A.opmult
isOpAdjoint(A::matrixOp)=A.isAdjoint
isOpAdjoint(A::compositeOp)=A.isAdjoint

struct vec_metadata
    n::Int32 #use Int32's for Fortran compatibility
    nvctr::Int32
    nvctrp::Int32
end

struct generalVec
    vec::Array{Complex128,2}
    isCovector::Bool
    md::vec_metadata
end

#constructor for converting true vectors (Array{T,1}) into the matrix format that generalVec expects (Array{T,2})
function generalVec(vec::Array{Complex128,1},isCovector::Bool)
    #md=vec_metadata(length(vec),1,1)
    #ccall((:set_metadata_,"./mymatmul.so"),Void,(Ref{Int32},Ref{Int32},Ref{Int32},Ref{vec_metadata_handle}),length(vec),1,1,Ref{vec_metadata_handle}(md_ptr))
    return generalVec(reshape(vec,length(vec),1),isCovector,md)
end

#convenience constructors for generalVec:
function generalVec(vec::Array{Complex128,2})
    md=vec_metadata(size(vec,1),size(vec,2),size(vec,2))
    return generalVec(vec,false,md)
end

function generalVec(vec::Array{Complex128,2},isCovector::Bool)
    md=vec_metadata(size(vec,1),size(vec,2),size(vec,2))
    return generalVec(vec,isCovector,md)
end

function generalVec(vec::Array{Complex128,1})
    return generalVec(reshape(vec,length(vec),1))
end

function userOp(options)

  mh=metadata_handle((0,0,0,0,0,0,0,0))

  ptr=options["ptr"]
  ccall((:operator_definition_,"./mymatmul.so"),Void,
  (Ref{Int32},Ref{Int32},Ref{Int32},Ref{Complex128},Ref{metadata_handle}),
  Ref{Int32}(options["n"]),
  Ref{Int32}(options["nvctr"]),
  Ref{Int32}(options["nvctrp"]),
  ptr,Ref{metadata_handle}(mh))
  #ccall((:op_direct_,"./mymatmul.so"),Void,(Ref{metadata_handle},),Ref{metadata_handle}(mh))

  function apply_to_userVec(x::generalVec)
  	  newvecs=zeros(Complex128,size(x.vec,1),size(x.vec,2))
          ccall((:apply_op_to_vec_,"./mymatmul_futile.so"),Void,(Ref{metadata_handle},Ref{Complex128},Ref{Complex128}),
          Ref{metadata_handle}(mh),x.vec,newvecs)

        return generalVec(newvecs,x.isCovector)
  end

  function apply_adjoint_to_userVec(x::generalVec)
  	  newvecs=zeros(Complex128,size(x.vec,1),size(x.vec,2))
	  ccall((:op_dagger,"./mymatmul_futile.so"),Void,(Ref{metadata_handle},),Ref{metadata_handle}(mh))
      ccall((:apply_op_to_vec_,"./mymatmul_futile.so"),Void,(Ref{metadata_handle},Ref{Complex128},Ref{Complex128}),
      Ref{metadata_handle}(mh),x.vec,newvecs)
	  ccall((:op_direct_,"./mymatmul_futile.so"),Void,(Ref{metadata_handle},),Ref{metadata_handle}(mh))

        return generalVec(newvecs,x.isCovector)
  end

  return userOp(ptr,mh,apply_to_userVec,apply_adjoint_to_userVec,false)
end


#constructor for matrixOp; define matrix vector multiplication here
#PLACEHOLDER; insert your routine here
function matrixOp(A::Array{Complex128,2})
    function myMatVec(x::generalVec)
        #Matrix multiply happens here! Replace the next lines with your own code
        #newvecs=A*x.vec[x.indices1,x.indices2]
	newvecs=zeros(Complex128,size(x.vec,1),size(x.vec,2))

        ccall((:mymatvec_,"./mymatmul.so"),Void,(Ref{Complex128},Ref{Complex128},Ref{Complex128},Ref{vec_metadata}),A,x.vec,newvecs,Ref{vec_metadata}(x.md))

        return generalVec(newvecs,x.isCovector)
    end

    function myAdjMatVec(x::generalVec)
        #Adjoint matrix multiply happens here! Replace the next lines with your own code
        #newvecs=A'*x.vec[x.indices1,x.indices2]
        newvecs=zeros(Complex128,size(x.vec,1),size(x.vec,2))

	ccall((:myadjmatvec_,"./mymatmul.so"),Void,(Ref{Complex128},Ref{Complex128},Ref{Complex128},Ref{vec_metadata}),A,x.vec,newvecs,Ref{vec_metadata}(x.md))

        return generalVec(newvecs,x.isCovector)
    end

    return matrixOp(A,myMatVec,myAdjMatVec,false)
end

#function for printing generalVecs
function Base.show(io::IO, x::generalVec)
    if(x.isCovector)
        print(io,x.vec')
    else
        print(io,x.vec)
    end
end


#~~~~~~~~~ Memory methods ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
function getindex(x::generalVec,k,l)

    if(k==Colon()) #then we're extracting a whole generalVec, or several of them
        ll=l
        if(l==Colon())
            return x
        else
            return generalVec(x.vec[k,l],x.isCovector)
        end
    else #otherwise we're extracting some sort of submatrix...disallowed for now
        error("Extracting coordinate submatrices from generalVec's is currently disallowed.")
    end
end

#vector of all zeros
#PLACEHOLDER; insert your own routine here
function zeros(x::generalVec)
    newvec=zeros(Complex128,size(x.vec,1),size(x.vec,2))
    return generalVec(newvec,x.isCovector)
end

#PLACEHOLDER; do your own memory allocations here
function copy(x::generalVec)
    #allocate new vector and copy coordinates from the old one
    newvec=Array{Complex128,2}(size(x.vec,1),size(x.vec,2))
    y=generalVec(newvec,x.isCovector)
    y[:,:]=x #all of the actual copying is done here; this calls setindex!()
    return y
end

#set all coordinates equal to a constant
#PLACEHOLDER; insert your own routine here
function setindex!(y::generalVec, a::Number,k)
    if(k==Colon())
        for i in 1:size(y,1)
            for j in 1:size(y,2)
                y.vec[i,j]=a
            end
        end
    else
        error("Indexing not fully implemented yet.")
    end
end

#assign array values, i.e. y[k,l]=x
#k and l can be ranges or values
#may need to be changed for custom data structures
function setindex!(y::generalVec,x::generalVec,k,l)
    kk=k
    ll=l
    if(k==Colon())
        kk=1:size(y.vec,1)
    end
    if(ll==Colon())
        ll=1:size(y.vec,2)
    end

    if(k!=Colon())
        error("Addressing submatrices of subspaces not supported yet.")
    end

    if(x.isCovector==y.isCovector)

        if(size(x.vec,1)==length(kk) && size(x.vec,2)==length(ll))
            for i in 1:length(kk)
                for j in 1:length(ll)
                    y.vec[kk[i],ll[j]]=x.vec[i,j]
                end
            end
        else
            error("Shape mismatch between input and output")
        end
    else
        error("Tried assigning a vector to a covector, or vice versa")
    end
end


#~~~~~~~~~~ Vector operations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#PLACEHOLDER; define your own function here
function size(x::generalVec)
    return size(x.vec)
end

function size(x::generalVec,k::Int)
    return size(x.vec,k)
end

function norm(x::generalVec)
    #PLACEHOLDER; INSERT YOUR ROUTINE HERE
    return norm(x.vec)
end

#perform scalar multiplications y=a*x
#allocate new vector, then multiply each of the coordinates of the old vector and store them in the new one
#important: explicitly only affects the submatrix of x.vec given by the indexes in x.indices1 and x.indices2
function *(a::Number, x::generalVec)
    #println("scalar-vector mul")
    #PLACEHOLDER; INSERT YOUR ROUTINE HERE
    newvecs=a*x.vec

    y=generalVec(newvecs,x.isCovector)
    return y
end

#scalar multiplication y=x*a
function *(x::generalVec,a::Number)
    return a*x
end

function /(x::generalVec,a::Number)
    #println("scalar-vector mul")
    #PLACEHOLDER; INSERT YOUR ROUTINE HERE
    newvecs=x.vec/a

    y=generalVec(newvecs,x.isCovector)
    return y
end

#conjugate transpose of vectors, i.e. X'
function ctranspose(x::generalVec)
    return generalVec(x.vec,!x.isCovector)
end


#block postmultiplication of subspace by a small dense matrix, i.e. y=x*A
#PROBABLY A PLACEHOLDER; insert your routine here
function *{T}(x::generalVec, A::Array{T,2})
    if(!x.isCovector)
        newvecs=x.vec*A
        return generalVec(newvecs,false)
    else #don't allow the user to do something that should be a generalOp multiplication
        error("Dense matrix multiplication is not supported-1.")
    end
end

#block postmultiplication of subspace by a small dense matrix, i.e. y=A*x'
#PROBABLY A PLACEHOLDER; insert your routine here
function *{T}(A::Array{T,2},x::generalVec)
    if(x.isCovector)
        newvecs=x.vec*A'
        return generalVec(newvecs,true)
    else #don't allow the user to do something that should be a generalOp multiplication
        error("Dense matrix multiplication is not supported-2.")
    end
end


#~~~~~~~~~~~ Vector - Vector operations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`

#perform vector addition z=x+y
function +(x::generalVec, y::generalVec)
    #make sure they're both vectors or both covectors, and that the size is the same
    if(size(x)==size(y) && x.isCovector == y.isCovector)
        #PLACEHOLDER; INSERT YOUR ROUTINE HERE
        newvecs=y.vec+x.vec

        return generalVec(newvecs,x.isCovector)
    else
        error("Tried to add two vectors of non-conforming sizes.")
    end
end


#perform vector subtraction z=x-y
function -(x::generalVec, y::generalVec)
    #make sure they're both vectors or both covectors, and that the size is the same
    if(size(x)==size(y) && x.isCovector == y.isCovector)
        #PLACEHOLDER; INSERT YOUR ROUTINE HERE
        newvecs=x.vec-y.vec

        return generalVec(newvecs,x.isCovector)
    else
        error("Tried to add two vectors of non-conforming sizes.")
    end
end


#perform block inner products x'*y
function *(x::generalVec, y::generalVec)
    if(x.isCovector && !y.isCovector)
        #PLACEHOLDER; INSERT YOUR ROUTINE HERE
        #retval=x.vec[x.indices1,x.indices2]'*y.vec[y.indices1,y.indices2]

	(n,m)=size(x)
        retval=zeros(Complex128,m,m)
        ccall((:myinnerproduct_,"./mymatmul.so"),Void,(Ref{Complex128},Ref{Complex128},Ref{Complex128},Ref{vec_metadata}),x.vec,y.vec,retval,Ref{vec_metadata}(x.md))

        if(size(retval)==(1,1))
            return retval[1] #return a scalar instead of a 1x1 julia matrix
        else
            return retval
        end
    elseif(x.isCovector == y.isCovector)
        error("Direct products not implemented.")
    else
        error("Outer products not implemented.")
    end
end

#~~~~~~~~~~Operator operations~~~~~~~~~~~~~~~~
#Represent linear combinations of operators through functional composition; makes it possible to write expressions like y=(z*A+B)*x and ensure that everything will be evaluated in terms of matrix-vector multiplication

#operator-operator multiplication, C=A*B
#probably don't have to change this
function *(A::generalOp, B::generalOp)
    am=opMult(A)
    aAdjm=adjOpMult(A)
    bm=opmult(B)
    bAdjm=adjOpMult(B)

    #if either operator is adjoint, then its multiplication routines are swapped:
    if(isOpAdjoint(A))
        am=adjOpMult(A)
        aAdjm=opMult(A)
    end
    if(isOpAdjoint(B))
        bm=adjOpMult(B)
        bAdjm=opmult(B)
    end

    #new left and right multiply functions:
    newmult(x)=am(bm(x))
    newAdjmult(x)=aAdjm(bAdjm(x))
    return compositeOp(newmult,newAdjmult,false)
end

#operator-operator addition, C=A+B
#probably don't have to change this
function +(A::generalOp, B::generalOp)
    am=opMult(A)
    aAdjm=adjOpMult(A)
    bm=opMult(B)
    bAdjm=adjOpMult(B)

    #if either operator is adjoint, then its multiplication routines are swapped:
    if(isOpAdjoint(A))
        am=adjOpMult(A)
        aAdjm=opMult(A)
    end
    if(isOpAdjoint(B))
        bm=adjOpMult(B)
        bAdjm=opmult(B)
    end

    newmult(x)=am(x)+bm(x)
    newAdjmult(x)=aAdjm(x)+bAdjm(x)

    return compositeOp(newmult,newAdjmult,false)
end

function -(A::generalOp, B::generalOp)
    am=opMult(A)
    aAdjm=adjOpMult(A)
    bm=opMult(B)
    bAdjm=adjOpMult(B)

    #if either operator is adjoint, then its multiplication routines are swapped:
    if(isOpAdjoint(A))
        am=adjOpMult(A)
        aAdjm=opMult(A)
    end
    if(isOpAdjoint(B))
        bm=adjOpMult(B)
        bAdjm=opmult(B)
    end

    newmult(x)=am(x)-bm(x)
    newAdjmult(x)=aAdjm(x)-bAdjm(x)

    return compositeOp(newmult,newAdjmult,false)
end


#scalar multiplication times operator, i.e. B=a*A
#probably don't need to change this; it will ultimately use the scalar-vector multiplication routine during function evaluation
function *(a::Number, A::generalOp)
    am=opMult(A)
    aAdjm=adjOpMult(A)
    if(isOpAdjoint(A))
        am=adjOpMult(A)
        aAdjm=opMult(A)
    end

    newmult(x)=a*am(x)
    newAdjmult(x)=a*aAdjm(x)

    return compositeOp(newmult,newAdjmult,false)
end

#more scalar multiplication; B=A*a
function *(A::generalOp, a::Number)
    return a*A
end

#scalar plus operator addition, i.e. B=a+A
#we'll implicitly put the identity operator next to a, i.e. B=aI+A
#probably don't need to change this
function +(a::Number, A::generalOp)
    am=opMult(A)
    aAdjm=adjOpMult(A)
    if(isOpAdjoint(A))
        am=adjOpMult(A)
        aAdjm=opMult(A)
    end

    newmult(x)=a*x+am(x)
    newAdjmult(x)=a*x+aAdjm(x)

    return compositeOp(newmult,newAdjmult,false)
end

#more scalar-operator addition:
function +(A::generalOp, a::Number)
    return a+A
end

#~~~~~~~~~~~~~~Operator-vector operations~~~~~~~~~~~~
#Operator-vector multiplication, i.e. y=A*x
function *(A::generalOp, x::generalVec)
    if(x.isCovector) #don't let user try to do dense matrix multiplication with generalOp
        error("Dense matrix multiplication disabled for generalOps")
    end

    am=opMult(A)
    aAdjm=adjOpMult(A)

    if(isOpAdjoint(A))
        return aAdjm(x)
    else
        return am(x)
    end
end

#covector-operator multiplication, i.e. y = x'*A
function *(x::generalVec,A::generalOp)
    if(!x.isCovector)#don't let user try to do dense matrix multiplication with generalOp
        error("Dense matrix multiplication disabled for generalOps2")
    end

    am=opMult(A)
    aAdjm=adjOpMult(A)

    if(isOpAdjoint(A))
        return am(x)
    else
        return aAdjm(x)
    end

end
