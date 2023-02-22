  use defined_types
  use neural
  implicit none
  integer, allocatable, dimension(:) :: size
  integer :: nlayers,n0,n1,n,m,mw,m1,m2,m0
  real(8) :: rand
  real(8) :: y,y0
  type(net) :: model
  real(8) :: tiny,numeric,analyt
  real(8), allocatable, dimension(:) :: dy1_dx,dy0_dx
  logical calc_second

  calc_second = .false.

  write(*,*) 'TESTING NEURAL NETWORK DERIVATIVES'
  write(*,*) 
  
  nlayers = 5
  
  !set the dimensions of the layers. 
  ! the 0th layer is the input layer
    allocate(size(0:nlayers))
    
    size(0) = 5
    size(1) = 3
    size(2) = 4
    size(3) = 5
    size(4) = 4
    size(5) = 1

    ! initialize model    

    call model%init(nlayers,size,calc_second)
    

    write(*,*) 'NETWORK CREATED WITH ',nlayers,' LAYERS'
    write(*,*)

    !set initial weights, biases, input vector with random numbers
    
    do n1 = 1,size(0)
       call random_number(rand)
       model%x0(n1) = rand - 0.5d0
    end do
    
    do n = 1,nlayers
       do n0 = 1,size(n-1)
          do n1 = 1,size(n)
             call random_number(rand)
             model%weight(n)%mat2(n1,n0) = rand - 0.5d0
          end do
       end do
       model%weight(n)%mat2 = model%weight(n)%mat2
       
       do n1 = 1,size(n) 
          call random_number(rand)
          model%bias(n)%vec(n1) = rand - 0.5d0
       end do
    end do

    allocate(dy0_dx(size(0)))
    allocate(dy1_dx(size(0)))
    
    call model%run_net
    y0 = model%y
    dy0_dx = model%dy_dx
    
    tiny = 1.d-5
    
    write(*,*) 'NUMERICAL/ANALYTICAL/RATIO'
    write(*,*) 

    do m = 1,size(0)
       model%x0(m) = model%x0(m) + tiny
       call model%run_net
       y = model%y
       model%x0(m) = model%x0(m) - tiny
       numeric = (y - y0)/tiny
       analyt = model%dy_dx(m)
       write(*,*) 'dy_dx',m,numeric,analyt,numeric/analyt
    end do
    write(*,*) 
    
    do mw = 1,nlayers
       do m2 = 1,size(mw)
          model%bias(mw)%vec(m2) = model%bias(mw)%vec(m2) + tiny
          call model%run_net
          
          y = model%y
          model%bias(mw)%vec(m2) = model%bias(mw)%vec(m2) - tiny
          numeric = (y - y0)/tiny
          
          analyt = model%dy_dbias(mw)%vec(m2)
          write(*,*) 'dy_db',mw,m2,numeric,analyt,numeric/analyt
       end do
    end do
    
    write(*,*) 

    do mw = 1,nlayers
       do m1 = 1,size(mw-1)
          do m2 = 1,size(mw)
             model%weight(mw)%mat2(m2,m1) = model%weight(mw)%mat2(m2,m1) + tiny
             call model%run_net
             y = model%y
             model%weight(mw)%mat2(m2,m1) = model%weight(mw)%mat2(m2,m1) - tiny
             numeric = (y - y0)/tiny
             
             analyt = model%dy_dweight(mw)%mat2(m2,m1)
             write(*,*) 'dy_dw',mw,m2,m1,numeric,analyt,numeric/analyt
          end do
       end do
    end do

    if(.not.calc_second) stop
    write(*,*) 

    do mw = 1,nlayers
       do m2 = 1,size(mw)
        model%bias(mw)%vec(m2) = model%bias(mw)%vec(m2) + tiny
        
        call model%run_net
        dy1_dx = model%dy_dx
        
        model%bias(mw)%vec(m2) = model%bias(mw)%vec(m2) - tiny
        
        do m0 = 1,size(0)
           numeric = (dy1_dx(m0) - dy0_dx(m0))/tiny
           
           analyt = model%d2y_dxdb(m0)%vecvec(mw)%vec(m2)
           write(*,*) 'd2y_dxdb',mw,m2,numeric,analyt,numeric/analyt
        end do
     end do
  end do
  
  write(*,*) 

  do mw = 1,nlayers
     do m1 = 1,size(mw-1)
        do m2 = 1,size(mw)
           model%weight(mw)%mat2(m2,m1) = model%weight(mw)%mat2(m2,m1) + tiny
           
           call model%run_net
           dy1_dx = model%dy_dx
           
           model%weight(mw)%mat2(m2,m1) = model%weight(mw)%mat2(m2,m1) - tiny
           
           do m0 = 1,size(0)
              numeric = (dy1_dx(m0) - dy0_dx(m0))/tiny
              
              analyt = model%d2y_dxdw(m0)%vecmat2(mw)%mat2(m2,m1)
              write(*,*) 'd2y_dxdw',mw,m2,m1,numeric,analyt,numeric/analyt
           end do
        end do
     end do
  end do


end program
