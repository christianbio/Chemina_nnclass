module neural
  use defined_types
  private
  type, public :: net
     type(mat2), allocatable, dimension(:) :: weight
     type(vec), allocatable, dimension(:) :: bias
     type(mat2), allocatable, dimension(:) :: dy_dweight !  (w1/w2/w3..) / (m1,m2)
     type(vec), allocatable, dimension(:) :: dy_dbias !  (w1/w2/w3..) / (m1)
     type(vecmat2), allocatable, dimension(:) :: d2y_dxdw, d2y_dxdw_ ! (m0) / (w1/w2..) / (m1,m2)
     type(vecvec), allocatable, dimension(:) :: d2y_dxdb, d2y_dxdb_ ! (m0) / (b1/b2..) / (m1)

     integer, dimension(:), allocatable :: size
     integer nlayers
     logical calc_second
     real(8) :: y
     real(8), allocatable, dimension(:) :: x0 
     real(8), allocatable, dimension(:) :: dy_dx

     contains 
       procedure, public :: init => init
       procedure, public :: run_net => run_net
       procedure, public :: dealloc => dealloc
  end type net

  contains
    subroutine init(this,n_layers,size_,calc_second)
      implicit none
      class(net) :: this
      integer :: n_layers,n,ll,m
      integer, dimension(0:n_layers) :: size_
      logical calc_second

      this%calc_second = calc_second

      !make sure the size of the last layer = 1
      size_(n_layers) = 1

      this%nlayers = n_layers
      this%size = size_

      allocate(this%x0(this%size(0)))
      allocate(this%dy_dx(this%size(0)))
      allocate(this%d2y_dxdw(this%size(0)))
      allocate(this%d2y_dxdb(this%size(0)))

      allocate(this%d2y_dxdw_(this%size(0)))
      allocate(this%d2y_dxdb_(this%size(0)))

      do m = 1,this%size(0)
         allocate(this%d2y_dxdw(m)%vecmat2(this%nlayers))
         allocate(this%d2y_dxdb(m)%vecvec(this%nlayers))

         allocate(this%d2y_dxdw_(m)%vecmat2(this%nlayers))
         allocate(this%d2y_dxdb_(m)%vecvec(this%nlayers))
         do n = 1,this%nlayers
            allocate(this%d2y_dxdw(m)%vecmat2(n)%mat2(this%size(n),this%size(n-1)))
            allocate(this%d2y_dxdb(m)%vecvec(n)%vec(this%size(n)))

            allocate(this%d2y_dxdw_(m)%vecmat2(n)%mat2(this%size(n),this%size(n-1)))
            allocate(this%d2y_dxdb_(m)%vecvec(n)%vec(this%size(n)))

            this%d2y_dxdw_(m)%vecmat2(n)%mat2 = 0.0d0 
            this%d2y_dxdb_(m)%vecvec(n)%vec = 0.0d0 
         end do
      end do


      allocate(this%weight(this%nlayers))
      allocate(this%bias(this%nlayers))

      do ll = 1,this%nlayers
         allocate(this%weight(ll)%mat2(this%size(ll),this%size(ll-1)))
      end do
      
      allocate(this%dy_dweight(this%nlayers))
      allocate(this%dy_dbias(this%nlayers))
      do n = 1,this%nlayers
         allocate(this%dy_dweight(n)%mat2(this%size(n),this%size(n-1)))
         allocate(this%bias(n)%vec(this%size(n)))
         allocate(this%dy_dbias(n)%vec(this%size(n)))
      end do

    end subroutine init

    subroutine dealloc(this)
      implicit none
      class(net) :: this
      integer n,ll,m

      do m = 1,this%size(0)
         do n = 1,this%nlayers
            deallocate(this%d2y_dxdw(m)%vecmat2(n)%mat2)
            deallocate(this%d2y_dxdb(m)%vecvec(n)%vec)

            deallocate(this%d2y_dxdw_(m)%vecmat2(n)%mat2)
            deallocate(this%d2y_dxdb_(m)%vecvec(n)%vec)
         end do

         deallocate(this%d2y_dxdw(m)%vecmat2)
         deallocate(this%d2y_dxdb(m)%vecvec)

         deallocate(this%d2y_dxdw_(m)%vecmat2)
         deallocate(this%d2y_dxdb_(m)%vecvec)
      end do

      deallocate(this%x0)
      deallocate(this%dy_dx)
      deallocate(this%d2y_dxdw)
      deallocate(this%d2y_dxdb)

      deallocate(this%d2y_dxdw_)
      deallocate(this%d2y_dxdb_)

      do ll = 1,this%nlayers
         deallocate(this%weight(ll)%mat2)
      end do

      deallocate(this%weight)

      do n = 1,this%nlayers
         deallocate(this%dy_dweight(n)%mat2)
         deallocate(this%bias(n)%vec)
         deallocate(this%dy_dbias(n)%vec)
      end do
      deallocate(this%bias)

      deallocate(this%dy_dweight)
      deallocate(this%dy_dbias)


    end subroutine dealloc


    subroutine run_net(this)
      implicit none
      class(net) :: this
      type(vec), dimension(0:this%nlayers) :: zz
      type(mat2), dimension(this%nlayers) :: dz_dzminus,d2z_dzminus2
      type(vecmat3), dimension(0:this%nlayers) :: dz_dweight    ! (layer) / (w1/w2/w3..) / (z_i,m1,m2)
      type(vecmat2), dimension(0:this%nlayers) :: dz_dbias    ! (layer) / (w1/w2/w3..) / (m1)
      type(vec), dimension(0:this%nlayers) :: dy_dz
      type(vec), dimension(this%nlayers) :: dsig
      real(8), dimension(this%size(0)) :: z0
      real(8), allocatable, dimension(:) :: z,zminus
      type(vecvecmat2), dimension(0:this%nlayers) :: d2y_dzdw ! (layer) / (mz) / (w1/w2/w3..) / (m1,m2)
      type(vecvecvec), dimension(0:this%nlayers) :: d2y_dzdb ! (layer) / (mz) / (b1/b2/b3..) / (m1)
      integer layer
      integer n,nminus,m,m2,m0,ll,m1      

      z0 = this%x0

      do layer = 0,this%nlayers
         allocate(zz(layer)%vec(this%size(layer)))
         allocate(dy_dz(layer)%vec(this%size(layer)))
         allocate(d2y_dzdw(layer)%vecvecmat2(this%nlayers))
         allocate(d2y_dzdb(layer)%vecvecvec(this%nlayers))
         do m = 1,this%nlayers
            allocate(d2y_dzdw(layer)%vecvecmat2(m)%vecmat2(this%size(layer)))
            allocate(d2y_dzdb(layer)%vecvecvec(m)%vecvec(this%size(layer)))
            do m2 = 1,this%size(layer)
               allocate(d2y_dzdw(layer)%vecvecmat2(m)%vecmat2(m2)%mat2(this%size(m),this%size(m-1)))
               allocate(d2y_dzdb(layer)%vecvecvec(m)%vecvec(m2)%vec(this%size(m)))
            end do
         end do
         
         if(layer.gt.0) then 
            allocate(dz_dweight(layer)%vecmat3(this%nlayers))
            allocate(dz_dbias(layer)%vecmat2(this%nlayers))
            allocate(dsig(layer)%vec(this%size(layer-1)))
            do m = 1,this%nlayers
               allocate(dz_dweight(layer)%vecmat3(m)%mat3(this%size(layer),this%size(m),this%size(m-1)))
               allocate(dz_dbias(layer)%vecmat2(m)%mat2(this%size(layer),this%size(m)))
            end do

            allocate(dz_dzminus(layer)%mat2(this%size(layer),this%size(layer-1)))
            allocate(d2z_dzminus2(layer)%mat2(this%size(layer),this%size(layer-1)))
         endif
      end do

      zz(0)%vec = z0

      do layer = 1,this%nlayers
         !forward propagation
         n = this%size(layer) 
         nminus = this%size(layer-1) 
         allocate(z(n))
         
         call forward(this%nlayers,layer,this%size,n,nminus,&
              zz(layer-1)%vec,this%weight(layer)%mat2,this%bias(layer)%vec,&
              dz_dweight(layer-1)%vecmat3,&
              dz_dbias(layer-1)%vecmat2,&
              z,dz_dzminus(layer)%mat2,d2z_dzminus2(layer)%mat2, &
              dsig(layer)%vec,&
              dz_dweight(layer)%vecmat3,&
              dz_dbias(layer)%vecmat2)
         
         zz(layer)%vec = z
         
         if(layer.eq.this%nlayers) then 
            this%y = z(1)
         endif
         
         deallocate(z)

      end do

      do n = 1,this%nlayers
         this%dy_dweight(n)%mat2(:,:) = dz_dweight(this%nlayers)%vecmat3(n)%mat3(1,:,:)
         this%dy_dbias(n)%vec(:) = dz_dbias(this%nlayers)%vecmat2(n)%mat2(1,:)
      end do

      dy_dz(this%nlayers)%vec(1) = 1.0d0

      do ll = 1,this%nlayers
         d2y_dzdw(this%nlayers)%vecvecmat2(ll)%vecmat2(1)%mat2 = 0.0d0 
         d2y_dzdb(this%nlayers)%vecvecvec(ll)%vecvec(1)%vec = 0.0d0 
      end do

      do layer = this%nlayers,1,-1
         n = this%size(layer) 
         nminus = this%size(layer-1) 
         
         call backward(layer,this%size,n,nminus,&
              this%nlayers,dz_dzminus(layer)%mat2,&
              dz_dweight(layer-1)%vecmat3,&
              dz_dbias(layer-1)%vecmat2,&
              dsig(layer)%vec,&
              d2y_dzdw(layer)%vecvecmat2,&
              d2y_dzdw(layer-1)%vecvecmat2,&
              d2y_dzdb(layer)%vecvecvec,&
              d2y_dzdb(layer-1)%vecvecvec,&
              d2z_dzminus2(layer)%mat2, &
              dy_dz(layer)%vec,dy_dz(layer-1)%vec,this%calc_second)
      end do


      this%dy_dx = dy_dz(0)%vec

      do ll = 1,this%nlayers
         do m0 = 1,this%size(0)
            this%d2y_dxdw(m0)%vecmat2(ll)%mat2 = d2y_dzdw(0)%vecvecmat2(ll)%vecmat2(m0)%mat2
            this%d2y_dxdb(m0)%vecvec(ll)%vec = d2y_dzdb(0)%vecvecvec(ll)%vecvec(m0)%vec
         end do
      end do

    end subroutine run_net


    subroutine forward(nlayers,layer,size,n,nminus,&
         zminus,w,b,&
         dzminus_dw,&
         dzminus_db,&
         z,dz_dzminus,d2z_dzminus2, &
         dsig,&
         dz_dw,&
         dz_db)
      use defined_types
      implicit none
      integer layer,n,nminus,nlayers
      real(8), dimension(n,nminus) :: w,dz_dzminus,d2z_dzminus2
      real(8), dimension(n) :: b
      real(8), dimension(nminus) :: zminus,sig,dsig,d2sig
      real(8), dimension(n) :: z
      integer, dimension(0:nlayers) :: size  
      integer m1,m2,m0,k0,k1,ll,m
      real(8) :: z0,dz

      type(mat3), dimension(nlayers) :: dzminus_dw
      type(mat3), dimension(nlayers) :: dz_dw
      type(mat2), dimension(nlayers) :: dz_db
      type(mat2), dimension(nlayers) :: dzminus_db

      do m = 1,nlayers
         dz_dw(m)%mat3 = 0.0d0
         dz_db(m)%mat2 = 0.0d0 
      end do
      
      if(layer.eq.1) then 
         z = matmul(w,zminus) + b 
         dz_dzminus = w 
         
         do m0 = 1,n
            dz_db(1)%mat2(m0,m0) = 1.0d0 
         end do

         do k0 = 1,nminus
            z0 = zminus(k0)
            do m0 = 1,n
               dz_dw(1)%mat3(m0,m0,k0) = z0
            end do
         end do
         dsig = 1.0d0 
      else
         call get_sigma(nminus,zminus,sig,dsig,d2sig)
         z = matmul(w,sig) + b
      endif

      if(layer.gt.1) then 
         do m1 = 1,nminus
            do m2 = 1,n
               dz_dzminus(m2,m1) = w(m2,m1) * dsig(m1)
               d2z_dzminus2(m2,m1) = w(m2,m1) * d2sig(m1)
            end do
         end do
      endif

      do ll = 1,nlayers
         if(ll.lt.layer) then 
            do m1 = 1,nminus
               do k1 = 1,size(ll)
                  dz = dzminus_db(ll)%mat2(m1,k1)
                  do m2 = 1,n
                     dz_db(ll)%mat2(m2,k1) = dz_db(ll)%mat2(m2,k1) + dz_dzminus(m2,m1) * dz
                  end do
               end do
            end do
            
            
            do k0 = 1,size(ll-1)
               do k1 = 1,size(ll)
                  do m1 = 1,nminus
                     dz = dzminus_dw(ll)%mat3(m1,k1,k0)
                     do m2 = 1,n
                        dz_dw(ll)%mat3(m2,k1,k0) = dz_dw(ll)%mat3(m2,k1,k0) + dz_dzminus(m2,m1) * dz
                     end do
                  end do
               end do
            end do
               
         endif
      end do
      
      if(layer.ne.1) then 
         do k1 = 1,n
            dz_db(layer)%mat2(k1,k1) = dz_db(layer)%mat2(k1,k1) + 1.0d0 
         end do
         do k0 = 1,nminus
            do k1 = 1,n
               dz_dw(layer)%mat3(k1,k1,k0) = sig(k0)
            end do
         end do
      endif
      
    end subroutine forward


    subroutine get_sigma(n,xvec, sigvec, dsigvec, d2sigvec)
      implicit none
      integer n
      real(8) :: efac, eplus, eminus
      real(8), dimension(n) :: xvec, sigvec, dsigvec, d2sigvec
      !tanh         
      !         sigvec = dtanh(xvec)
      !         dsigvec = 1.0d0 - sigvec**2
      !         d2sigvec = -2.0d0*sigvec*dsigvec
      !logistic
      sigvec = 1.0d0 / (1.0d0 + dexp(-xvec))
      dsigvec = sigvec * (1.0d0 - sigvec)
      d2sigvec = dsigvec * (1.0d0 -2.0d0 * sigvec)
      
    end subroutine get_sigma

    subroutine backward(layer,size,n,nminus,&
         nlayers,dz_dzminus,&
         dzminus_dw,&
         dzminus_db,&
         dsig,&
         d2y_dzdw, & 
         d2y_dzminusdw,&
         d2y_dzdb, & 
         d2y_dzminusdb,&
         d2z_dzminus2, &
         dy_dz,dy_dzminus,calc_second)
      use defined_types
      implicit none
      integer layer,nlayers,n,nminus
      integer, dimension(0:nlayers) :: size  
      
      real(8), dimension(nminus) :: dsig
      type(mat3), dimension(nlayers) :: dzminus_dw
      type(mat2), dimension(nlayers) :: dzminus_db
      real(8), intent(in), dimension(n,nminus) :: dz_dzminus
      
      real(8), dimension(n,nminus) :: d2z_dzminus2

      real(8), dimension(n) :: dy_dz
      real(8), dimension(nminus) :: dy_dzminus
      integer m1,m2,k1,k2,ll
      type(vecmat2), dimension(nlayers) :: d2y_dzdw
      type(vecmat2), dimension(nlayers) :: d2y_dzminusdw
      type(vecvec), intent(in), dimension(nlayers) :: d2y_dzdb
      type(vecvec), dimension(nlayers) :: d2y_dzminusdb
      real(8) :: dsigma,dd,dd1,dd2,dk1,dz1,dz2
      real(8), allocatable, dimension(:,:) :: amat
      logical calc_second

      dy_dzminus = matmul(dy_dz, dz_dzminus)  
      
      if(.not.calc_second) return

      do k1 = 1,nminus
         do ll = 1,nlayers
            d2y_dzminusdb(ll)%vecvec(k1)%vec = 0.0d0 
            d2y_dzminusdw(ll)%vecmat2(k1)%mat2 = 0.0d0 
         end do
      end do
      
      do ll = 1,nlayers
         if(layer.gt.ll) then 

            do k1 = 1,nminus
               do k2 = 1,n
                  dz1 = d2z_dzminus2(k2,k1) * dy_dz(k2)
                  dz2 = dz_dzminus(k2,k1)
                  do m2 = 1,size(ll)
                     d2y_dzminusdb(ll)%vecvec(k1)%vec(m2) = d2y_dzminusdb(ll)%vecvec(k1)%vec(m2) + &
                       dzminus_db(ll)%mat2(k1,m2) * dz1 & 
                      + d2y_dzdb(ll)%vecvec(k2)%vec(m2) * dz2
                  end do
               end do
            end do

            allocate(amat(size(ll),size(ll-1)))

            do k1 = 1,nminus
               amat(:,:) = dzminus_dw(ll)%mat3(k1,:,:)
               do k2 = 1,n
                  d2y_dzminusdw(ll)%vecmat2(k1)%mat2 = d2y_dzminusdw(ll)%vecmat2(k1)%mat2 + &
                       amat * dy_dz(k2)  * d2z_dzminus2(k2,k1)  & 
                       + d2y_dzdw(ll)%vecmat2(k2)%mat2 * dz_dzminus(k2,k1)
               end do
            end do
            
            deallocate(amat)

         else !i.e. if(layer.le.ll)

            !expensive loop!

            if(layer.lt.nlayers) then 
               do k1 = 1,nminus
                  do k2 = 1,n
                     d2y_dzminusdb(ll)%vecvec(k1)%vec = d2y_dzminusdb(ll)%vecvec(k1)%vec + &
                          d2y_dzdb(ll)%vecvec(k2)%vec * dz_dzminus(k2,k1)
                     d2y_dzminusdw(ll)%vecmat2(k1)%mat2 = d2y_dzminusdw(ll)%vecmat2(k1)%mat2 + &
                          d2y_dzdw(ll)%vecmat2(k2)%mat2 * dz_dzminus(k2,k1)
                  end do
               end do
            endif
            
            
         endif
      end do
   
      do k1 = 1,nminus
         dsigma = dsig(k1)
         do m2 = 1,size(layer)
            d2y_dzminusdw(layer)%vecmat2(k1)%mat2(m2,k1) = d2y_dzminusdw(layer)%vecmat2(k1)%mat2(m2,k1) + & 
                 dy_dz(m2) * dsigma
         end do
      end do
      

    end subroutine backward


  end module neural

