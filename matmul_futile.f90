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
    use wrapper_MPI
    implicit none

    if (.not. initialized) then
      call f_lib_initialize()
      initialized=.true.
      call mpiinit()
      print *,'Hello world, I am',mpirank()
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

    !call f_routine(id='hamiltonian_init')
    H%n=n
    H%nvctr=nvctr
    H%nvctrp=nvctrp
    H%ptr=f_malloc_ptr([n,n],id='ptr')
    call f_memcpy(src=ptr,dest=H%ptr)
    !!call f_release_routine()
  end subroutine hamiltonian_init

  subroutine hamiltonian_free(H)
    use dynamic_memory
    implicit none
    type(hamiltonian), intent(inout) :: H
    !!call f_routine(id='hamiltonian_free')
    call ensure_initialization()
    call f_free_ptr(H%ptr)
    H=hamiltonian_null()
    !!call f_release_routine()

  end subroutine hamiltonian_free

  subroutine hamiltonian_application(H,x,y,m)
    implicit none
    integer, intent(in) :: m
    type(hamiltonian), intent(in) :: H
    complex(dp), dimension(H%n,H%nvctr), intent(in) :: x
    !complex(dp), dimension(H%n,H%m), intent(inout) :: y
    complex(dp), dimension(H%n,H%nvctr), intent(inout) :: y

    call ensure_initialization()
    !!call f_routine(id='hamiltonian_application')
    !example of parallel implementantion
    !just use zgemm for now
    !y_vector=f_malloc([H%n,H%nvctrp], id='y_vector')
    !call fmpi_alltoall(y,sendcount=H%np*H%nvctr/mpisize(),recvbuf=y_vector)

    !call zgemm(H%transa,'N',H%n,m_parallel ,H%n,&
    !1.0d0,H%ptr,H%n,x,H%n,0.0d0,y,H%n)

    !call fmpi_alltoall(recvbuf=y,sendcount=H%np*H%nvctr/mpisize(),sendbuf=y_vector)


    call zgemm(H%transa,'N',H%n,m,H%n,&
    1.0d0,H%ptr,H%n,x,H%n,0.0d0,y,H%n)


    !!call f_release_routine()
  end subroutine hamiltonian_application

  subroutine conjugate_hamiltonian(H)
    implicit none
    type(hamiltonian), intent(inout) :: H
    call ensure_initialization()
    !!call f_routine(id='conjugate_hamiltonian')
    !print *,'inside conjugate hamiltonian'
    H%transa='C'
    !!call f_release_routine()
  end subroutine conjugate_hamiltonian

  subroutine direct_hamiltonian(H)
    implicit none
    type(hamiltonian), intent(inout) :: H
    call ensure_initialization()
    !!call f_routine(id='direct_hamiltonian')
    !print *,'inside direct hamiltonian'
    H%transa='N'
    !!call f_release_routine()
  end subroutine direct_hamiltonian

end module wrapper

subroutine op_create(n,nvctr,nvctrp,ptr,A)
  use wrapper
  implicit none
  integer, intent(in) :: n,nvctr,nvctrp
  complex(dp), dimension(n,n), intent(in) :: ptr
  type(operator_container), intent(inout) :: A
  !local variables
  type(operator_container), dimension(12) :: Avec

  print *,f_loc(Avec(4))-f_loc(Avec(3)),'size'

  allocate(A%H)
  call hamiltonian_init(n,nvctr,nvctrp,ptr,A%H)
end subroutine op_create

subroutine op_dagger(A)
  use wrapper
  implicit none
  type(operator_container), intent(inout) :: A

  call conjugate_hamiltonian(A%H)
end subroutine op_dagger

subroutine op_direct(A)
  use  wrapper
  implicit none
  type(operator_container), intent(inout) :: A
  call direct_hamiltonian(A%H)
end subroutine op_direct

subroutine finalize_wrapper()
  use wrapper_mpi
  implicit none
  call mpifinalize()
  call f_lib_finalize()
end subroutine finalize_wrapper

subroutine apply_op_to_vec(A,x,y,m)
  use wrapper
  implicit none
  integer, intent(in) :: m
  type(operator_container), intent(in) :: A
  complex(dp), dimension(*) :: x
  complex(dp), dimension(*) :: y

  !here on can control if the containers are allocated
  if (.not. associated(A%H)) then
     call f_err_throw('The Hamiltionian operator is not correctly associated')
     return
  end if
  call hamiltonian_application(A%H,x,y,m)

end subroutine apply_op_to_vec

subroutine op_release(A)
  use wrapper
  implicit none
  type(operator_container), intent(inout) :: A
  call hamiltonian_free(A%H)
  deallocate(A%H)
end subroutine op_release

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
