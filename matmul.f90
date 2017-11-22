module blaswrapper_setup
  implicit none

  type, public :: metaparams
     integer :: n=0 !<size of matrix A and B (subspace)
     integer :: nvctr=0 !< number of componenets for each of the vector of the subspace
     integer :: nvctrp=0
  end type metaparams


  type, public :: hamiltonian
     integer :: n=0 !<size of matrix A and B (subspace)
     integer :: nvctr=0 !< number of componenets for each of the vector of the subspace
     integer :: nvctrp=0 !, local amount of componenet stored by the processor (for parallel implmentation)
     character(len=1) :: transa='N'
     complex*16, dimension(:,:), pointer :: ptr => null()
  end type hamiltonian

  type, public :: operator_container
     type(hamiltonian), pointer :: h
  end type operator_container

contains

  function hamiltonian_null() result(H)
    implicit none
    type(hamiltonian) :: H
  end function hamiltonian_null

  subroutine hamiltonian_init(n,nvctr,nvctrp,ptr,H)
    implicit none
    integer, intent(in) :: n,nvctr,nvctrp
    complex*16, dimension(n,n), intent(in) :: ptr
    type(hamiltonian), intent(out) :: H

    H%n=n
    H%nvctr=nvctr
    H%nvctrp=nvctrp
    allocate(H%ptr(n,n))
    H%ptr=ptr

  end subroutine hamiltonian_init

  subroutine hamiltonian_free(H)
    implicit none
    type(hamiltonian), intent(inout) :: H

    deallocate(H%ptr)
    H=hamiltonian_null()

  end subroutine hamiltonian_free

  subroutine hamiltonian_application(H,x,y)
    implicit none
    type(hamiltonian), intent(in) :: H
    complex*16, dimension(*) :: x
    complex*16, dimension(*) :: y

    !example of parallel implementantion
    !just use zgemm for now
    call zgemm(H%transa,'N',H%n,H%nvctr,H%n,&
    1.0d0,H%ptr,H%n,x,H%n,0.0d0,y,H%n)
  end subroutine hamiltonian_application

  subroutine conjugate_hamiltonian(H)
    implicit none
    type(hamiltonian), intent(inout) :: H

    H%transa='C'
  end subroutine conjugate_hamiltonian

  subroutine direct_hamiltonian(H)
    implicit none
    type(hamiltonian), intent(inout) :: H
    H%transa='N'
  end subroutine direct_hamiltonian

end module blaswrapper_setup

subroutine operator_definition(n,nvctr,nvctrp,ptr,A)
  use blaswrapper_setup
  implicit none
  integer, intent(in) :: n,nvctr,nvctrp
  complex*16, dimension(n,n), intent(in) :: ptr
  type(operator_container), intent(inout) :: A

  allocate(A%H)
  call hamiltonian_init(n,nvctr,nvctrp,ptr,A%H)
  print *,'operator allocated, size',n,associated(A%H)
end subroutine operator_definition

subroutine op_dagger(A)
  use blaswrapper_setup
  implicit none
  type(operator_container), intent(inout) :: A

  call conjugate_hamiltonian(A%H)
end subroutine op_dagger

subroutine op_direct(A)
  use blaswrapper_setup
  implicit none
  type(operator_container), intent(inout) :: A
  print *,'op direct',associated(A%H)
  call direct_hamiltonian(A%H)
end subroutine op_direct


subroutine apply_op_to_vec(A,x,y)
  use blaswrapper_setup
  implicit none
  type(operator_container) :: A
  complex*16, dimension(*) :: x
  complex*16, dimension(*) :: y

  !here on can control if the containers are allocated
  if (.not. associated(A%H)) then
     print *,'Error, the operator is not associated'
     stop
  end if
  print *,'calling the hamiltonian'
  call hamiltonian_application(A%H,x,y)

end subroutine apply_op_to_vec

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
