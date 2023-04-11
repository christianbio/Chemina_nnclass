module data
  use defined_types
  use neural
  integer nlayers
  type(net) :: model
  integer, allocatable, dimension(:) :: size
  type(mat2), allocatable, dimension(:) :: de2_dw
  type(vec), allocatable, dimension(:) :: de2_db
  real(8) :: e2
  integer nminsteps
end module data

  use defined_types
  use neural
  use data
  implicit none
  logical calc_second
  real(8) :: rand,error2
  integer n,n1,n0
  
  calc_second = .true.

  nlayers = 4

  allocate(size(0:nlayers))
    
  size(0) = 1
  size(1) = 8
  size(2) = 5
  size(3) = 2
  size(4) = 1

   call init_random_seed()
   

  allocate(de2_dw(nlayers))
  allocate(de2_db(nlayers))
  do n = 1,nlayers
     allocate(de2_dw(n)%mat2(size(n),size(n-1)))
     allocate(de2_db(n)%vec(size(n)))
  end do
  

  call model%init(nlayers,size,calc_second)

  do n1 = 1,size(0)
     call random_number(rand)
     model%x0(n1) = rand - 0.5d0
  end do
  
  do n = 1,nlayers
     do n0 = 1,size(n-1)
        do n1 = 1,size(n)
           call random_number(rand)
           model%weight(n)%mat2(n1,n0) = (rand - 0.5d0)*0.1d0
        end do
     end do
     model%weight(n)%mat2 = model%weight(n)%mat2
     
     do n1 = 1,size(n) 
        call random_number(rand)
        model%bias(n)%vec(n1) = (rand - 0.5d0)*0.1d0
     end do
  end do
  
  call model%run_net

!  call get_error(error2)
!  write(*,*) 'error2 = ',error2
!  call test
  call minimize

end program


subroutine get_error(error2)
  use data
  implicit none
  integer i,imax,n
  real(8) :: xmin,xmax,x,ytarget,ymodel,error2,delta

  do n = 1,nlayers
     de2_dw(n)%mat2 = 0.0d0 
     de2_db(n)%vec = 0.0d0 
  end do 

  imax = 200
  xmin = 0.0d0 
  xmax = 1.0d0 

  error2 = 0.0d0 
  do i = 0,imax-1
     x = xmin + (xmax - xmin) * dble(i)/dble(imax-1)
     model%x0(1) = x
     call get_func(x,ytarget)
     call model%run_net
     ymodel = model%y
     delta = ymodel - ytarget
     error2 = error2 + delta**2
     do n = 1,nlayers
        de2_dw(n)%mat2 = de2_dw(n)%mat2 + 2.0d0 * delta * model%dy_dweight(n)%mat2
        de2_db(n)%vec = de2_db(n)%vec + 2.0d0 * delta * model%dy_dbias(n)%vec
     end do
     
  end do

  e2 = error2 
end subroutine get_error

subroutine test
  use data
  implicit none
  integer n,i,j,i1,i2
  real(8) :: e0,e1,tiny
  real(8) :: numeric,analyt

  tiny = 1.0d-4

  call get_error(e0)
  do n = 1,nlayers
     do i1 = 1,size(n)
        do i2 = 1,size(n-1)
           model%weight(n)%mat2(i1,i2) = model%weight(n)%mat2(i1,i2) + tiny
           call get_error(e1)           
           numeric = (e1 - e0)/tiny
           analyt = de2_dw(n)%mat2(i1,i2)
           model%weight(n)%mat2(i1,i2) = model%weight(n)%mat2(i1,i2) - tiny
           write(*,*) 'www',n,i1,i2,numeric,analyt,numeric/analyt
        end do
     end do 

     do i1 = 1,size(n) 
        model%bias(n)%vec(i1) = model%bias(n)%vec(i1) + tiny
        call get_error(e1)           
        numeric = (e1 - e0)/tiny
        analyt = de2_db(n)%vec(i1)
        model%bias(n)%vec(i1) = model%bias(n)%vec(i1) - tiny
        write(*,*) 'bbb',n,i1,i2,numeric,analyt,numeric/analyt
     end do
  end do

end subroutine test

subroutine get_func(x,y)
  implicit none
  real(8) :: x,y
  y = dsin(25.0d0 * x)
end subroutine get_func

subroutine minimize
  use data
  use nr_common_data
  use nr_mod
  implicit none
  real(8) :: xvec(60000)
  real(8) :: ftol,fret
  integer n,i1,i2,ll,iter
  

  n =  0
  do ll = 1,nlayers
     do i1 = 1,size(ll)
        do i2 = 1,size(ll-1)
           n = n + 1
           xvec(n) = model%weight(ll)%mat2(i1,i2)
        end do
     end do
     do i1 = 1,size(ll)
        n = n + 1
        xvec(n) = model%bias(ll)%vec(i1)
     end do 
  end do

  ftol = 0.0d0
  linmin_param = 1.0d-6 * 1.d+3
  nprintdata = 1
  
  frprmn_itmax = 8192*8192*8
  
  nminsteps = 0 
  call frprmn(xvec, n, ftol, iter, fret)
  
end subroutine minimize


real(8) function func(xvec)
  use data
  use nr_common_data
  use nr_mod
  implicit none
  real(8) :: xvec(60000)
  real(8) :: error2
  integer n,i1,i2,ll

  n =  0
  do ll = 1,nlayers
     do i1 = 1,size(ll)
        do i2 = 1,size(ll-1)
           n = n + 1
           model%weight(ll)%mat2(i1,i2) = xvec(n)
        end do
     end do
     do i1 = 1,size(ll)
        n = n + 1
        model%bias(ll)%vec(i1) = xvec(n) 
     end do 
  end do

  call get_error(error2) 
  func = error2

end function func

subroutine dfunc(xvec, gvec)
  use data
  use nr_common_data
  use nr_mod
  implicit none
  real(8) :: xvec(60000),gvec(60000) 
  integer n,i1,i2,ll

  n =  0
  do ll = 1,nlayers
     do i1 = 1,size(ll)
        do i2 = 1,size(ll-1)
           n = n + 1
           gvec(n) = de2_dw(ll)%mat2(i1,i2)
        end do
     end do
     do i1 = 1,size(ll)
        n = n + 1
        gvec(n) = de2_db(ll)%vec(i1) 
     end do 
  end do


end subroutine dfunc

subroutine print_data
  use data
  implicit none

  nminsteps = nminsteps + 1
  write(*,*) 'net',nminsteps,e2
  call flush(6)
end subroutine print_data

       subroutine init_random_seed()
         use iso_fortran_env, only: int64
         implicit none
         integer, allocatable :: seed(:)
         integer :: i, n, un, istat, dt(8), pid
         integer(int64) :: t

         call random_seed(size=n)
         allocate (seed(n))
         ! First try if the OS provides a random number generator
         open(newunit=un, file="/dev/urandom", access="stream", &
               form="unformatted", action="read", status="old", iostat=istat)
         if (istat == 0) then
            read(un) seed
            close (un)
         else
            ! Fallback to XOR:ing the current time and pid. The PID is
            ! useful in case one launches multiple instances of the same
            ! program in parallel.
            call system_clock(t)
            if (t == 0) then
               call date_and_time(values=dt)
               t = (dt(1) - 1970)*365_int64*24*60*60*1000 &
                   + dt(2)*31_int64*24*60*60*1000 &
                   + dt(3)*24_int64*60*60*1000 &
                   + dt(5)*60*60*1000 &
                   + dt(6)*60*1000 + dt(7)*1000 &
                   + dt(8)
            end if
            pid = getpid()
            t = ieor(t, int(pid, kind(t)))
            do i = 1, n
               seed(i) = lcg(t)
            end do
         end if
         call random_seed(put=seed)
      contains
         ! This simple PRNG might not be good enough for real work, but is
         ! sufficient for seeding a better PRNG.
         function lcg(s)
            integer :: lcg
            integer(int64) :: s
            if (s == 0) then
               s = 104729
            else
               s = mod(s, 4294967296_int64)
            end if
            s = mod(s*279470273_int64, 4294967291_int64)
            lcg = int(mod(s, int(huge(0), int64)), kind(0))
         end function lcg
      end subroutine init_random_seed
