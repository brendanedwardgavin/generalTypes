module wrapper
  use f_precisions, dp=>f_double
  use dynamic_memory, only: f_routine,f_release_routine
  implicit none
  type, public :: hamiltonian
     integer :: n=0 !<size of matrix A and B (subspace)
     integer :: nvctr=0 !< number of componenets for each of the vector of the subspace
     integer :: nvctrp=0 !, local amount of componenet stored by the processor (for parallel implmentation)
     character(len=1) :: transa='N'
     complex(dp), dimension(:,:), pointer :: ptr => null()
  end type hamiltonian

  type, public :: operator_container
     type(hamiltonian), pointer :: h
  end type operator_container

  logical, save :: initialized=.false.

contains

  subroutine ensure_initialization()
    implicit none

    if (.not. initialized) then
      call f_lib_initialize()
      initialized=.true.
    end if
  end subroutine ensure_initialization

  function hamiltonian_null() result(H)
    implicit none
    type(hamiltonian) :: H
  end function hamiltonian_null

  subroutine hamiltonian_init(n,nvctr,nvctrp,ptr,H)
    use dynamic_memory
    implicit none
    integer, intent(in) :: n,nvctr,nvctrp
    complex(dp), dimension(n,n), intent(in) :: ptr
    type(hamiltonian), intent(out) :: H

    call ensure_initialization()

    call f_routine(id='hamiltonian_init')
    H%n=n
    H%nvctr=nvctr
    H%nvctrp=nvctrp
    H%ptr=f_malloc_ptr([n,n],id='ptr')
    call f_memcpy(src=ptr,dest=H%ptr)
    call f_release_routine()
  end subroutine hamiltonian_init

  subroutine hamiltonian_free(H)
    use dynamic_memory
    implicit none
    type(hamiltonian), intent(inout) :: H

    call ensure_initialization()
    call f_free_ptr(H%ptr)
    H=hamiltonian_null()

  end subroutine hamiltonian_free

  subroutine hamiltonian_application(H,x,y)
    implicit none
    type(hamiltonian), intent(in) :: H
    complex(dp), dimension(H%n,H%nvctr), intent(in) :: x
    complex(dp), dimension(H%n,H%nvctr), intent(inout) :: y

    call ensure_initialization()
    call f_routine(id='hamiltonian_application')
    !example of parallel implementantion
    !just use zgemm for now
    call zgemm(H%transa,'N',H%n,H%nvctr,H%n,&
    1.0d0,H%ptr,H%n,x,H%n,0.0d0,y,H%n)
    call f_release_routine()
  end subroutine hamiltonian_application

  subroutine conjugate_hamiltonian(H)
    implicit none
    type(hamiltonian), intent(inout) :: H
    call ensure_initialization()
    call f_routine(id='conjugate_hamiltonian')
    H%transa='C'
    call f_release_routine()
  end subroutine conjugate_hamiltonian

  subroutine direct_hamiltonian(H)
    implicit none
    type(hamiltonian), intent(inout) :: H
    call ensure_initialization()
    call f_routine(id='direct_hamiltonian')
    H%transa='N'
    call f_release_routine()
  end subroutine direct_hamiltonian

end module wrapper

subroutine operator_definition(n,nvctr,nvctrp,ptr,A)
  use wrapper
  implicit none
  integer, intent(in) :: n,nvctr,nvctrp
  complex(dp), dimension(n,n), intent(in) :: ptr
  type(operator_container), intent(inout) :: A

  allocate(A%H)
  call hamiltonian_init(n,nvctr,nvctrp,ptr,A%H)
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
  call direct_hamiltonian(A%H)
end subroutine op_direct

subroutine apply_op_to_vec(A,x,y)
  use wrapper
  implicit none
  type(operator_container), intent(in) :: A
  complex(dp), dimension(*) :: x
  complex(dp), dimension(*) :: y

  !here on can control if the containers are allocated
  if (.not. associated(A%H)) then
     call f_err_throw('The Hamiltionian operator is not correctly associated')
     return
  end if
  call hamiltonian_application(A%H,x,y)

end subroutine apply_op_to_vec

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
