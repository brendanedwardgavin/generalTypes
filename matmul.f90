module blaswrapper_setup
  implicit none
  type, public :: feastarray_params
     integer :: n=0 !<size of matrix A and B (subspace)
     integer :: nvctr=0 !< number of componenets for each of the vector of the subspace
     integer :: nvctrp=0 !, local amount of componenet stored by the processor (for parallel implmentation)
  end type

  type(feastarray_params), save, public ::  vector_metadata
end module

subroutine set_blaswrapper_objects(n,nvctr,nvctrp)
use blaswrapper_setup, only: vector_metadata
implicit none
    integer :: n,nvctr,nvctrp
    vector_metadata%n=n
    vector_metadata%nvctr=nvctr
    vector_metadata%nvctrp=nvctrp
end subroutine


!matrix multiplication y=A*x
subroutine mymatvec(A,x,y)
use blaswrapper_setup, only: vector_metadata
implicit none
    complex*16, dimension(*) :: A
    complex*16, dimension(*) :: x
    complex*16, dimension(*) :: y
    integer :: n,m
    !the vector_metadata shoud of course be initialized
    n=vector_metadata%n
    m=vector_metadata%nvctrp

    !example of parallel implementantion
    !just use dgemm for now
    call zgemm('N','N',n,m,n,1.0d0,A,n,x,n,0.0d0,y,n)

    !call mpiallred(y(1),nvctrp,FMPI_SUM)

end subroutine mymatvec

!matrix multiplication y=A'*x
subroutine myadjmatvec(A,x,y)
use blaswrapper_setup, only: vector_metadata
implicit none
    complex*16, dimension(*) :: A
    complex*16, dimension(*) :: x
    complex*16, dimension(*) :: y
    integer :: n,m
    !the vector_metadata shoud of course be initialized
    n=vector_metadata%n
    m=vector_metadata%nvctrp

    !example of parallel implementantion
    !just use dgemm for now
    call zgemm('C','N',n,m,n,1.0d0,A,n,x,n,0.0d0,y,n)

    !call mpiallred(y(1),nvctrp,FMPI_SUM)

end subroutine myadjmatvec


!a=x'*y
subroutine myInnerProduct(x,y,a)
use blaswrapper_setup, only: vector_metadata
implicit none
    complex*16, dimension(*) :: x
    complex*16, dimension(*) :: y !n x m block vectors
    complex*16, dimension(*) :: a !m x m matrix of inner products
    integer :: n,m

    n=vector_metadata%n
    m=vector_metadata%nvctrp

    !just use dgemm for now
    call zgemm('C','N',m,m,n,1.0d0,x,n,y,n,0.0d0,a,m)

end subroutine myInnerProduct
