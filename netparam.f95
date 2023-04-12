module netparam
  use defined_types
  private
  type, public :: net_param
     type(mat2), allocatable, dimension(:) :: weight
     type(vec), allocatable, dimension(:) :: bias

     type(mat2), allocatable, dimension(:) :: de2_dw
     type(vec), allocatable, dimension(:) :: de2_db
     integer, dimension(:), allocatable :: size
     integer nlayers

     contains 
       procedure, public :: init => init
       procedure, public :: dealloc => dealloc
  end type net_param

  contains
    subroutine init(this,n_layers,size_)
      implicit none
      class(net_param) :: this
      integer :: n_layers,n,mw,m
      integer, dimension(0:n_layers) :: size_

      !make sure the size of the last layer = 1
      size_(n_layers) = 1

      this%nlayers = n_layers
      this%size = size_

      allocate(this%weight(this%nlayers))
      allocate(this%de2_dw(this%nlayers))
      allocate(this%bias(this%nlayers))
      allocate(this%de2_db(this%nlayers))

      do mw = 1,this%nlayers
         allocate(this%weight(mw)%mat2(this%size(mw),this%size(mw-1)))
         allocate(this%bias(mw)%vec(this%size(mw)))

         allocate(this%de2_dw(mw)%mat2(this%size(mw),this%size(mw-1)))
         allocate(this%de2_db(mw)%vec(this%size(mw)))
      end do
      
    end subroutine init

    subroutine dealloc(this)
      implicit none
      class(net_param) :: this
      integer n,mw,m


      do mw = 1,this%nlayers
         deallocate(this%weight(mw)%mat2)
         deallocate(this%de2_dw(mw)%mat2)
      end do

      deallocate(this%weight)

      do n = 1,this%nlayers
         deallocate(this%bias(n)%vec)
         deallocate(this%de2_db(n)%vec)
      end do
      deallocate(this%bias)


    end subroutine dealloc


  end module netparam

