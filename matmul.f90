module blaswrapper_setup
  implicit none

  type, public :: metaparams
     integer :: n=0 !<size of matrix A and B (subspace)
     integer :: nvctr=0 !< number of componenets for each of the vector of the subspace
     integer :: nvctrp=0
  end type metaparams

end module blaswrapper_setup

!matrix multiplication y=A*x
subroutine mymatvec(A,x,y,vec_metadata)
use blaswrapper_setup
implicit none
    complex*16, dimension(*) :: A
    complex*16, dimension(*) :: x
    complex*16, dimension(*) :: y
    type(metaparams) :: vec_metadata
    integer :: n,m
    !the vector_metadata shoud of course be initialized
    n=vec_metadata%n
    m=vec_metadata%nvctr

    !example of parallel implementantion
    !just use dgemm for now
    call zgemm('N','N',n,m,n,1.0d0,A,n,x,n,0.0d0,y,n)

    !call mpiallred(y(1),nvctrp,FMPI_SUM)

end subroutine mymatvec

!matrix multiplication y=A'*x
subroutine myadjmatvec(A,x,y,vec_metadata)
use blaswrapper_setup
implicit none
    complex*16, dimension(*) :: A
    complex*16, dimension(*) :: x
    complex*16, dimension(*) :: y
    type(metaparams) :: vec_metadata
    integer :: n,m
    !the vector_metadata shoud of course be initialized
    n=vec_metadata%n
    m=vec_metadata%nvctr

    !example of parallel implementantion
    !just use dgemm for now
    call zgemm('C','N',n,m,n,1.0d0,A,n,x,n,0.0d0,y,n)

    !call mpiallred(y(1),nvctrp,FMPI_SUM)

end subroutine myadjmatvec


!a=x'*y
subroutine myInnerProduct(x,y,a,vec_metadata)
use blaswrapper_setup
implicit none
    complex*16, dimension(*) :: x
    complex*16, dimension(*) :: y !n x m block vectors
    complex*16, dimension(*) :: a !m x m matrix of inner products
    type(metaparams) :: vec_metadata
    integer :: n,m

    n=vec_metadata%n
    m=vec_metadata%nvctr

    !just use dgemm for now
    call zgemm('C','N',m,m,n,1.0d0,x,n,y,n,0.0d0,a,m)

end subroutine myInnerProduct
