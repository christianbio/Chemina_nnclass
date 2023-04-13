program fitnn
   use molneural_data
   implicit none
   character(len=32), dimension(1024) :: control_list
   real(8), dimension(1024, 8) :: param_list
   character(2048), dimension(1024) :: string_list
   character(64) :: filename
   integer fileno
   integer ncontrol
   logical file_exists
   integer :: istart,ilen

   calc_forces = .true.


   minmode = "none"
   use_sumtype = .false.
   
   ii0 = 0
   
   neighbor0 = 0
   mnet_allocated = .false.
   
   min_error_val = 1.d+23
   nminsteps = 0 
   arc_cut = 100.0d0    

   rmin_b = 0.0d0
   rmax_b = 4.0d0
   
   rmin_nb = 0.0d0
   rmax_nb = 16.0d0
   
   nbinmax = 128
   additive = .false.
   
   max_col = 256

   mnet_alpha = 20.0d0 
   calc_valid = .true.
   restart_calc_valid = .false.

   call check_files_exist
   
   open(10,file = 'pf.xyz')
   open(11,file = 'efile')
   
   open(22, file='mnet.dat')
   open(32, file = 'cnet.dat')
   
   call open_append_file(20, 'mnet.out')
   call open_append_file(30, 'cnet.out')
   
   call open_append_file(24, 'mnet.best')
   call open_append_file(34, 'cnet.best')
   
   call open_append_file(20, 'mnet.out')
   
   open(100, file='nncontrol.dat')
   
   first_read = .true.
   
   aufac = 2625.0d0
   aufac_f = 2625.0d0/0.52918d0
   
   evfac = 96.487d0
   evfac_f = 96.487d0

   evfac = 1.0d0
   evfac_f = 1.0d0 

   read_mnet_flag = .false.
   read_cnet_flag = .false.
   
   itmax = 4096
   randval = 0.0d0 
   radial_cut = 3.0d0 

   call init_random_seed()

   call read_control(control_list, param_list, string_list, ncontrol)

   call set_units
   
   print_forces = .false.
   call find_nsteps_total
   call read_natoms

   call read_efile

   call calc_variances
   
   call allocate_initial
   
   call allocate_cols

   call read_train_data
   call count_atoms_and_molecules

   call calc_coulomb
   call shift_energies

   !must come after read_train_data
   call find_atlist
   call find_used_cols
   
   net_type = "mnet"
   call allocate(net_type)
   if(.not.additive.and.(.not.restart)) call randomize(randval, additive)

   if(read_mnet_flag) call read_mnet0
   if(read_cnet_flag) call read_cnet

   if(additive.and.(.not.restart)) call randomize(randval, additive)
   call print_mnet(20)

   if(run_cnet_flag) net_type = 'cnet'
   if(run_mnet_flag) net_type = 'mnet'

   if(restart) call restart_log
   call run
   
   stop

end program fitnn


module networks_mod
   use omp_lib
   use molneural_data
   use nr_common_data
   use neural
   implicit none

   CONTAINS

      subroutine run_train_cnet
         implicit none
         integer ns,m1,m2,id1,id2,nc_used,etype,id_center,cindex,cindex_atom,ncols,id
         integer iatom,itype
         real(8) :: delta_energy,energy_pred
         real(8), allocatable, dimension(:) :: charge_pred

         real(8), allocatable, dimension(:,:) :: force_pred, delta_force

         real(8), allocatable, dimension(:,:) :: denergy_dw1
         real(8), allocatable, dimension(:,:,:) :: denergy_dw2,denergy_dw3
         real(8), allocatable, dimension(:,:) :: denergy_db1,denergy_db2,denergy_db3
         real(8), allocatable, dimension(:) :: denergy_dconst

         real(8), allocatable, dimension(:,:,:,:) :: denergy_dw
         real(8), allocatable, dimension(:,:,:) :: denergy_db

         real(8), allocatable, dimension(:,:,:,:) :: dcharge_dw1
         real(8), allocatable, dimension(:,:,:,:) :: dcharge_dw2,dcharge_dw3
         real(8), allocatable, dimension(:,:,:) :: dcharge_db1,dcharge_db2,dcharge_db3

         real(8), allocatable, dimension(:,:,:,:,:) :: dcharge_dw
         real(8), allocatable, dimension(:,:,:,:) :: dcharge_db

         real(8), allocatable, dimension(:,:,:,:) :: dforce_dw1

         real(8), allocatable, dimension(:,:,:,:,:) :: dforce_dw2, dforce_dw3
         real(8), allocatable, dimension(:,:,:,:) :: dforce_db1,dforce_db2,dforce_db3

         real(8), allocatable, dimension(:,:,:,:,:,:) :: dforce_dw
         real(8), allocatable, dimension(:,:,:,:,:) :: dforce_db

         real(8) :: delta_charge,frac,cfac
         real(8) :: net_timing
         character(3) :: atname

         character(16) :: net_type

!         write(*,*) 'RUN TRAIN CNET'

         cnet_de2_db1 = 0.0d0
         cnet_de2_db2 = 0.0d0
         cnet_de2_dw2 = 0.0d0         

         error2_charge = 0.0d0 
         error2_charge_val = 0.0d0 


         do id_center = 1,natom_list
            do id1 = 0,natom_list
               do id2 = id1,natom_list
                  do etype = 0,3
                     if(col_used(etype,id2,id1,id_center)) then
                        nc_used = ncol_used(etype,id2,id1,id_center)
                        do m1 = 1,nc_used
                           cindex = col_index(etype,m1,id2,id1,id_center)
                           do m2 = 1,mnet_s1
                              cnet_de2_dw1(m2,cindex) = 0.0d0 
                           end do
                        end do
                     endif
                  end do
               end do
            end do
         end do

         if(minmode.ne.'energy') then 
            error2_c_reg = 0.0d0 
            call get_cnet_reg         
         endif

         !$OMP PARALLEL DEFAULT(SHARED) &
         !$OMP PRIVATE(energy_pred, charge_pred, & 
         !$OMP force_pred, delta_energy, delta_charge,delta_force, atname, itype, m2, iatom, cindex,cindex_atom, &
         !$OMP dcharge_dw1, dcharge_dw2, dcharge_db1, dcharge_db2, &
         !$OMP denergy_dw1, denergy_dw2, denergy_db1, denergy_db2, &
         !$OMP denergy_db, denergy_dw, &
         !$OMP dforce_dw1, dforce_dw2, dforce_db1,cfac,frac)


         allocate(dcharge_dw1(mnet_size(1),ncols_atom_max,natom_list,natoms_max))
         allocate(dcharge_db1(mnet_size(1), natom_list, natoms_max))
         allocate(charge_pred(natoms_max))


         !$OMP DO REDUCTION(+:error2_energy, error2_force, error2_charge, &
         !$OMP error2_energy_val, error2_force_val, &
         !$OMP mnet_de2_dw1, & 
         !$OMP mnet_de2_dw, mnet_de2_db,net_timing)

         do ns = 1,nsteps_total

            net_type = 'cnet'

            call get_mnet_energy(ns, net_type, &
                 energy_pred, charge_pred, &
                 denergy_dw1, denergy_dw, denergy_db, & 
                 dcharge_dw1, dcharge_dw, dcharge_db, &
                 dforce_dw1, dforce_dw, dforce_db, & 
                 denergy_dconst,force_pred,net_timing)
            
            !calculate error

            do iatom = 1,natoms_step(ns)

               delta_charge = charge_pred(iatom) - charge(iatom,ns)
               charge0(iatom,ns) = charge_pred(iatom)

               if (train_step(ns)) then
                  
                  atname = atname_step(iatom,ns)
                  call get_atom_id(atname,itype)

                  !weight the contribution to error2_energy by the fraction of atoms of this element
                  frac = dble(charge_count(itype))/dble(ncharges_total)
                  cfac = (1.0d0 / dble(ncharges_train(itype)))/ var_charge(itype)
                  error2_charge = error2_charge + frac * cfac * (delta_charge)**2 
                  
                  
                  do itype = 1,natom_list
                     do id1 = 0,natom_list
                        do id2 = id1,natom_list
                           do etype = 0,3
                              if(col_used(etype,id2,id1,itype)) then
                              nc_used = ncol_used(etype,id2,id1,itype)
                              do m1 = 1,nc_used
                                 cindex = col_index(etype,m1,id2,id1,itype)
                                 cindex_atom = col_atom_index(etype,m1,id2,id1,itype)
                                 do m2 = 1,mnet_s1
                                    cnet_de2_dw1(m2,cindex) = cnet_de2_dw1(m2,cindex) + &
                            2.0d0 * delta_charge * dcharge_dw1(m2,cindex_atom,itype,iatom) * frac * cfac
                                  end do
                               end do
                            endif
                         end do
                      end do
                   end do
                end do
                
                
                do itype = 1,natom_list
                   do m1 = 1,mnet_size(1)
                      cnet_de2_db1(m1,itype) = cnet_de2_db1(m1,itype) &
                             + 2.0d0*delta_charge*dcharge_db1(m1,itype,iatom) * frac * cfac
                      do m2 = 1,mnet_size(2)
                         cnet_de2_dw2(m2,m1,itype) = cnet_de2_dw2(m2,m1,itype) &
                              + 2.0d0*delta_charge*dcharge_dw2(m2,m1,itype,iatom) * frac * cfac
                      end do 
                   end do
 
                   !work on later
!                   cnet_de2_db2(itype) = cnet_de2_db2(itype) &
!                        + 2.0d0 * delta_charge * dcharge_db2(itype,iatom) * frac * cfac
                end do
                
               else !i.e. validation

                  atname = atname_step(iatom,ns)
                  call get_atom_id(atname,itype)

                  frac = dble(charge_count(itype))/dble(ncharges_total)
                  cfac = (1.0d0 / dble(ncharges_valid(itype)))/ var_charge(itype)

                  error2_charge_val = error2_charge_val + frac * cfac * (delta_charge)**2 
               endif !if (train_step(ns)) then
            end do
         end do
         


         !$OMP END DO

         !$OMP END PARALLEL

!         write(*,*) 'error2_charge = ',error2_charge

      end subroutine run_train_cnet

     subroutine run_train_mnet
         implicit none
         integer ns, iatom, ivec, itype, i, j
         real(8) :: ffac, efac, ffac_val,efac_val,delta_energy,nfac
         integer layer,m,m1,m2,m3,id1,id2,id_center,nc_used,cindex,etype,struct_id
         integer frag_id,nfrag,nf,ns2
         integer, allocatable, dimension(:) :: perm_step
         
         real(8) :: energy_pred,frac,var_est
         real(8), allocatable, dimension(:) :: charge_pred
         real(8), allocatable, dimension(:,:) :: force_pred, delta_force

         real(8), allocatable, dimension(:,:) :: denergy_dw1
         real(8), allocatable, dimension(:,:,:) :: denergy_dw2,denergy_dw3
         real(8), allocatable, dimension(:,:) :: denergy_db1,denergy_db2,denergy_db3
         real(8), allocatable, dimension(:) :: denergy_dconst

         real(8), allocatable, dimension(:,:,:,:) :: denergy_dw
         real(8), allocatable, dimension(:,:,:) :: denergy_db
         

         real(8), allocatable, dimension(:,:,:,:) :: dcharge_dw1
         real(8), allocatable, dimension(:,:,:,:) :: dcharge_dw2,dcharge_dw3
         real(8), allocatable, dimension(:,:,:) :: dcharge_db1,dcharge_db2,dcharge_db3

         real(8), allocatable, dimension(:,:,:,:,:) :: dcharge_dw
         real(8), allocatable, dimension(:,:,:,:) :: dcharge_db

         real(8), allocatable, dimension(:,:,:,:) :: dforce_dw1

         real(8), allocatable, dimension(:,:,:,:,:) :: dforce_dw2,dforce_dw3
         real(8), allocatable, dimension(:,:,:,:) :: dforce_db1,dforce_db2,dforce_db3

         real(8), allocatable, dimension(:,:,:,:,:,:) :: dforce_dw
         real(8), allocatable, dimension(:,:,:,:,:) :: dforce_db

         real(8) :: mult
         real(8) :: net_timing
         real(8) :: const,rand
         real(8), dimension(8192) :: energy_store         

         !only calculate the validation error once every 10 steps
         calc_valid = .true.

         if(minmode.eq.'none') then 
            energy_pred_store = 0.0d0 
         endif

         if(minmode.eq.'energy') then 
            do ns = 1,nsteps_total
               energy_store(ns) = energy_pred_store(ns)
            end do
         endif

         if(mod(nminsteps+1,20).ne.0) calc_valid = .false.
         if(nminsteps.eq.0) calc_valid = .true.
         if(restart_calc_valid) then
            calc_valid = .true.
            restart_calc_valid = .false.
         endif

         mnet_de2_db = 0.0d0 
         mnet_de2_dw = 0.0d0 

         mnet_de2_dconst = 0.0d0 

         net_timing = 0.0d0 

         do id_center = 1,natom_list
            do id1 = 0,natom_list
               do id2 = id1,natom_list
                  do etype = 0,3
                     if(col_used(etype,id2,id1,id_center)) then
                        nc_used = ncol_used(etype,id2,id1,id_center)
                        do m1 = 1,nc_used
                           cindex = col_index(etype,m1,id2,id1,id_center)
                           do m2 = 1,mnet_size(1)
                              mnet_de2_dw1(m2,cindex) = 0.0d0 
                           end do
                        end do
                     endif
                  end do
               end do
            end do
         end do
         !charge term
         do itype = 1,natom_list
            cindex = charge_col_index(itype)
            do m2 = 1,mnet_size(1)
               mnet_de2_dw1(m2,cindex) = 0.0d0 
            end do
         end do 
         
         allocate(perm_step(nsteps_total))
         perm_step = 0 
         ns = 1
         do while(.true.)
            call random_number(rand)
            ns2 = int(rand * nsteps_total) + 1
            if(perm_step(ns2).eq.0) then 
               perm_step(ns2) = ns
!               write(*,*) 'aa',ns,ns2
               ns = ns + 1
            endif
            if(ns.eq.nsteps_total+1) exit
         end do

         error2_force = 0.0d0
         error2_energy = 0.0d0

         if(calc_valid) then
            error2_force_val = 0.0d0
            error2_energy_val = 0.0d0
         endif
         
         if(minmode.ne.'energy') then 
            error2_m_reg = 0.0d0 
            call get_mnet_reg         
         endif

         mnet_size_max = 0 
         do m = 1,mnet_nlayers
            if(mnet_size(m).gt.mnet_size_max) mnet_size_max = mnet_size(m)
         end do


         !$OMP PARALLEL DEFAULT(SHARED) &
         !$OMP PRIVATE(energy_pred, force_pred, delta_energy, delta_force, &
         !$OMP         denergy_dw1, denergy_dw2, denergy_dw3,&
         !$OMP         denergy_db1, denergy_db2, denergy_db3, &
         !$OMP         denergy_db, denergy_dw, &  
         !$OMP         denergy_dconst, &
         !$OMP         dforce_dw1, dforce_dw2, dforce_dw3, dforce_db1,dforce_db2, &
         !$OMP         dforce_dw, dforce_db, &
         !$OMP         var_est,mult,const,frag_id,nf,nfrag,struct_id,&
         !$OMP         efac,frac,ns,ns2,nfac)


         allocate(denergy_dw1(mnet_size(1),ncols_wmat1))
         
         allocate(denergy_db(mnet_size_max,natom_list,mnet_nlayers))
         allocate(denergy_dw(mnet_size_max,mnet_size_max,natom_list,mnet_nlayers))

         allocate(denergy_dconst(128))

         allocate(dforce_dw1(3,ncols_wmat1,mnet_size(1),natoms_max))
         
         allocate(force_pred(3, natoms_max))
         allocate(delta_force(3, natoms_max))

         allocate(dforce_dw(3,mnet_size_max,mnet_size_max, natoms_max, natom_list,mnet_nlayers))
         allocate(dforce_db(3,mnet_size_max, natoms_max, natom_list,mnet_nlayers))

         !$OMP DO REDUCTION(+:error2_energy, error2_force, &
         !$OMP error2_energy_val, error2_force_val, mnet_de2_dw1, &
         !$OMP mnet_de2_dw, mnet_de2_db,mnet_de2_dconst,net_timing,energy_pred_store)

         do ns = 1, nsteps_total
            ns2 = perm_step(ns)
            nfrag = nfrag_ns(ns2) 

            if((.not.calc_valid).and.(.not.train_step(ns2))) then
               cycle
            endif


            calc_forces = .false.
            if(minmode.eq.'all'.or.minmode.eq.'none') then 
               calc_forces = .true.
               call get_mnet_energy(ns2, net_type, &
                    energy_pred, charge_pred, &
                    denergy_dw1, denergy_dw, denergy_db, & 
                    dcharge_dw1, dcharge_dw, dcharge_db, &
                    dforce_dw1, dforce_dw, dforce_db, & 
                    denergy_dconst,force_pred,net_timing)
               energy_pred_store(ns2) = energy_pred_store(ns2) + energy_pred
            else
            if(minmode.eq.'energy') energy_pred = energy_store(ns2)
         endif
            


            denergy_dconst = 0.0d0 
            do nf = 1,nfrag_ns(ns2)
               frag_id = frag_id_ns(ns2,nf)
               mult = frag_mult_ns(ns2,nf)
               const = mnet_const(frag_id)
               
               energy_pred = energy_pred + const * mult
               denergy_dconst(frag_id) = denergy_dconst(frag_id) + mult
            end do

            train_energy1(ns2) = energy_pred
            
            ffac = (1.0d0/dble(3*natoms_train))/var_force
            ffac_val = (1.0d0 /dble(3*natoms_val))/var_force

            nfac = 1.0d0 / dble(natoms_step(ns2))
            delta_energy = (energy_pred - train_energy(ns2))*nfac

            var_est = 0.0d0 
               do nf = 1,nfrag
                  mult = frag_mult_ns(ns2,nf)
                  frag_id = frag_id_ns(ns2,nf)
                  var_est = var_est + dble(mult) * var_energy(frag_id)
               end do

            if (train_step(ns2)) then

               struct_id = struct_id_ns(ns2)
               !estimated variance for this structure 

               efac = (1.0d0/dble(nsteps_train_struct_id(struct_id)))/var_est
               !weight the contribution to error2_energy by the fraction of frames containing this molecule.
               frac = dble(struct_id_count(struct_id))/dble(nsteps_total)
               
               do id_center = 1,natom_list
                  do id1 = 0,natom_list
                     do id2 = id1,natom_list
                        do etype = 0,3
                           if(col_used(etype,id2,id1,id_center)) then
                              nc_used = ncol_used(etype,id2,id1,id_center)
                              do m1 = 1,nc_used
                                 cindex = col_index(etype,m1,id2,id1,id_center)                                 
                                 do m2 = 1,mnet_size(1)
                                    mnet_de2_dw1(m2,cindex) = mnet_de2_dw1(m2,cindex) + &
                                         2.0d0 * delta_energy * denergy_dw1(m2,cindex)*efac*nfac*frac
                                 end do
                              end do
                           endif
                        end do
                     end do
                  end do
               end do
               !charge term
               do id_center = 1,natom_list
                  cindex = charge_col_index(id_center)
                  do m2 = 1,mnet_size(1)
                     mnet_de2_dw1(m2,cindex) = mnet_de2_dw1(m2,cindex) + &
                          2.0d0 * delta_energy * denergy_dw1(m2,cindex)*efac*nfac*frac
                  end do
               end do


               !now do the other mnet arrays

               do layer = 1,mnet_nlayers
                  do itype = 1, natom_list
                     if(layer.lt.mnet_nlayers) then 
                        do m1 = 1, mnet_size(layer)
                           mnet_de2_db(m1,itype,layer) = mnet_de2_db(m1,itype,layer) &
                                &            + 2.0d0*delta_energy*denergy_db(m1,itype,layer)*efac*nfac*frac
                        end do
                     endif

                     if(layer.gt.1) then 
                        do m1 = 1,mnet_size(layer-1)
                           do m2 = 1,mnet_size(layer)
                              mnet_de2_dw(m2,m1,itype,layer) = mnet_de2_dw(m2,m1,itype,layer) &
                                   &            + 2.0d0*delta_energy*denergy_dw(m2,m1,itype,layer)*efac*nfac*frac
                           end do
                        end do
                     endif
                  end do
               end do

                  do nf = 1,nfrag_ns(ns2)
                     frag_id = frag_id_ns(ns2,nf)               
                     mnet_de2_dconst(frag_id) = mnet_de2_dconst(frag_id) &
                          &            + 2.0d0*delta_energy*denergy_dconst(frag_id)*efac*nfac*frac
                  end do 
                  
                  error2_energy = error2_energy + efac*delta_energy**2 * frac
               else

               struct_id = struct_id_ns(ns2)
               efac_val = (1.0d0 /dble(nsteps_val_struct_id(struct_id)))/var_est
               frac = dble(struct_id_count(struct_id))/dble(nsteps_total)


               error2_energy_val = error2_energy_val + efac_val*delta_energy**2 * frac
            end if !if(train_step(ns2))


            if(calc_forces) then 

            !*****************************************************'
            !***     FORCES 
            !*****************************************************'

            do iatom = 1, natoms_step(ns2)
               do ivec = 1, 3
                  delta_force(ivec, iatom) = (force_pred(ivec, iatom) &
                       - train_force(ivec, iatom, ns2))
                  train_force1(ivec, iatom, ns2) = force_pred(ivec, iatom)

                  if (train_step(ns2)) then
                     error2_force = error2_force + ffac*delta_force(ivec, iatom)**2
                  else
                     error2_force_val = error2_force_val + ffac_val*delta_force(ivec, iatom)**2
                  end if
               end do
            end do

            if (train_step(ns2)) then

               do iatom = 1,natoms_step(ns2)
                  do ivec = 1,3
                     do id_center = 1,natom_list
                        do id1 = 0,natom_list
                           do id2 = id1,natom_list
                              do etype = 0,3
                                 if(col_used(etype,id2,id1,id_center)) then 
                                    nc_used = ncol_used(etype,id2,id1,id_center)
                                    do m1 = 1,nc_used
                                       cindex = col_index(etype,m1,id2,id1,id_center)
                                       do m2 = 1,mnet_size(1)
                                          mnet_de2_dw1(m2,cindex) = mnet_de2_dw1(m2,cindex) &
                                               + 2.0d0*delta_force(ivec, iatom)*dforce_dw1(ivec,cindex,m2,iatom)*ffac
                                       end do
                                    end do
                                 endif
                              end do
                           end do
                        end do
                     end do

!charges
                     do id_center = 1,natom_list
                        cindex = charge_col_index(id_center)
                        do m2 = 1,mnet_size(1)
                           mnet_de2_dw1(m2,cindex) = mnet_de2_dw1(m2,cindex) &
                                + 2.0d0*delta_force(ivec, iatom)*dforce_dw1(ivec,cindex,m2,iatom)*ffac
                        end do
                     end do

                  end do
               end do
               

               !now do the other arrays

               do layer = 1,mnet_nlayers
                  do itype = 1, natom_list
                     do iatom = 1, natoms_step(ns2)
                        do ivec = 1,3                 
                           if(layer.lt.mnet_nlayers) then 
                              do m1 = 1, mnet_size(layer)
                                 mnet_de2_db(m1,itype,layer) = mnet_de2_db(m1,itype,layer) &
                                      + 2.0d0*delta_force(ivec,iatom)*dforce_db(ivec,m1,iatom,itype,layer)*ffac
                              end do
                           endif

                           if(layer.gt.1) then 
                              do m1 = 1,mnet_size(layer-1)
                                 do m2 = 1,mnet_size(layer)
                                    mnet_de2_dw(m2,m1,itype,layer) = mnet_de2_dw(m2,m1,itype,layer) &
                                         + 2.0d0*delta_force(ivec, iatom)*dforce_dw(ivec,m2,m1,iatom,itype,layer)*ffac
                                 end do
                              end do

                           endif
                           
                        end do
                     end do
                  end do ! iatom
               end do ! itype
               endif

            end if

         end do
         !$OMP END DO

         !$OMP END PARALLEL

         
         error_train = error2_force + error2_energy + error2_m_reg 
         error_val = error2_force_val + error2_energy_val
         m_reg_error_train = error2_m_reg

!         write(*,*) 'TIME',net_timing
       end subroutine run_train_mnet

       
       subroutine  get_mnet_energy(ns, net_type, &
                 energy_pred, charge_pred, &
                 denergy_dw1, denergy_dw, denergy_db, & 
                 dcharge_dw1, dcharge_dw, dcharge_db, &
                 dforce_dw1, dforce_dw, dforce_db, & 
                 denergy_dconst,force_pred,net_timing)

         implicit none
         ! Dummy arguments
         integer, intent(in) :: ns
         real(8), intent(out) ::  energy_pred
         real(8), allocatable, dimension(:) :: charge_pred,charge0_pred
         real(8), intent(out), dimension(:,:,:,:) :: dforce_dw1         
         
         real(8), intent(out), dimension(:,:) :: denergy_dw1
         real(8), allocatable, dimension(:,:,:) :: denergy_dw2,denergy_dw3
         real(8), allocatable, dimension(:,:) :: denergy_db1,denergy_db2,denergy_db3
         real(8), intent(out), dimension(:) :: denergy_dconst

         real(8), allocatable, dimension(:,:,:,:) :: denergy_dw
         real(8), allocatable, dimension(:,:,:) :: denergy_db

         real(8), allocatable, dimension(:,:,:,:) :: dcharge_dw1
         real(8), allocatable, dimension(:,:,:,:) :: dcharge_dw2,dcharge_dw3
         real(8), allocatable, dimension(:,:,:) :: dcharge_db1,dcharge_db2,dcharge_db3

         real(8), allocatable, dimension(:,:,:,:,:) :: dcharge_dw
         real(8), allocatable, dimension(:,:,:,:) :: dcharge_db

         real(8), allocatable, dimension(:,:,:,:) :: dcharge0_dw1
         real(8), allocatable, dimension(:,:,:,:) :: dcharge0_dw2,dcharge0_dw3
         real(8), allocatable, dimension(:,:,:) :: dcharge0_db1,dcharge0_db2,dcharge0_db3

         real(8), intent(out), dimension(:,:) :: force_pred
         real(8), allocatable, dimension(:,:,:,:,:) :: dforce_dw2,dforce_dw3
         real(8), allocatable, dimension(:,:,:,:) :: dforce_db1,dforce_db2,dforce_db3

         real(8), allocatable, dimension(:,:,:,:,:,:) :: dforce_dw
         real(8), allocatable, dimension(:,:,:,:,:) :: dforce_db

         ! real(8), intent(out), dimension(:,:,:) :: dforce_db2
         real(8), allocatable, dimension(:,:) :: evector0,evector1,evector_qr,evector_qrplus,evector_qrminus
         real(8), allocatable, dimension(:) :: evector_data
         real(8), allocatable, dimension(:) :: zz,evalue0,evalue1,evalue_qr,evalue_qrplus,evalue_qrminus
         integer, allocatable, dimension(:) :: eigentype,dis_type,pair_index,kplus_index,kminus_index
         real(8), allocatable, dimension(:, :) :: w1
         real(8), allocatable, dimension(:) :: x1vec, b1vec, w2
         real(8), allocatable, dimension(:) :: x1vec_0,x1vec_1,x1vec_qr,x1vec_qrplus,x1vec_qrminus
         real(8) :: b2vec
         integer n,nc_used,struct_id,frag_id,kk
         real(8) :: couple
         real(8), allocatable, dimension(:) :: dy_dx1,dy_dz1,dy_dm
         real(8), allocatable, dimension(:,:) :: dy_dw1
         real(8), allocatable, dimension(:) :: dy_dw2
         real(8), allocatable, dimension(:) :: dy_db1
         real(8), allocatable, dimension(:,:,:) :: d2y_dx1dw1,d2y_dz1dw1
         real(8), allocatable, dimension(:,:,:) :: d2y_dx1dw2
         real(8), allocatable, dimension(:,:) :: d2y_dx1db1,d2y_dz1db1
         real(8), allocatable, dimension(:) :: d2y_dx1db2,d2y_dz1db2
         ! Local variables
         real(8) :: dy_db2
         real(8) :: yvec
         real(8) :: xdif,ydif,zdif,rdis,rdis2,ev2,arc_cut2
         integer m0,m1,m2,m3,m10,m11,iatom,jatom,layer
         real(8), dimension(0:11,512) :: xyz
         integer :: pair_atno(2,512) 
         real(8) :: rdisi
         integer ndis,ndis0,ndis1,mm,nn,ndis_central,ndis_non_central
         integer iat,jat
         real(8) :: aa,aa0
         character(3), dimension(512) :: atname_sorted
         integer, dimension(512) :: atindex_sorted
         integer, dimension(0:128) :: range
         integer nrange,n1,n2,n1_start,n1_end,n2_start,n2_end,k1,k2,npair,j
         integer n_types,itype,ipair
         integer, dimension(128) :: npairs_of_type
         integer i1,i2,ivec,itype0,jtype0
         integer natoms0,k,atomid
         integer, dimension(512) :: atindex
         character(3), dimension(512) :: atname
         character(3) :: aname
         real(8) :: d_evaluei_d_evalue
         real(8) :: rdis_1,rdis_2,rdis_1i,rdis_2i
         real(8) :: xdif_1,ydif_1,zdif_1,xdif_2,ydif_2,zdif_2
         real(8) :: aa1,aa2,aa_arc
         real(8) :: efac_1,efac_2,efac,efac_arc
         integer, dimension(128) :: atomid_list
         integer iatomid,jatomid,atomid_central
         integer maxid,minid
         integer, dimension(0:2,512) :: pairlist
         integer id1,id2,id_center,cindex,n_types0,n0
         integer iat0,jat0,ipair0,nn1,nn2,iat1,jat1,iat2,jat2
         real(8) :: alpha,beta,beta_arc,ecut
         integer jj,ndis2,ndis3,nd,id,nf,m
         integer itype1,itype2,etype1,etype2,nd1,nd2,p1,p2,p,etype
         integer nev_data,nev_data_max,ev_index,cindex_atom,mult
         integer nplus,nminus,kplus,kminus
         character(16) :: net_type
         real(8) :: charge_total,qi,qj         
         logical inside_cutoff,use_qrtype
         character(16) :: sort_direction
         type(net) :: model   
         integer, allocatable, dimension(:) :: size
         real(8) :: yvec2,tiny
         real(8), dimension(3) :: rr,rr_1,rr_2
         real(8) :: net_timing,start,end

         use_qrtype = .true.

         couple = 0.1d0 
         alpha = mnet_alpha
         beta = radial_cut - 0.5d0 
         beta_arc = arc_cut - 0.5d0 
         ecut = 1.d-4
         energy_pred = 0.0d0 

         arc_cut2 = arc_cut**2

         !zero forces

         !only need these forces for the training stage.

         if(net_type.eq.'mnet') then 
            denergy_dw1 = 0.0d0 
            dforce_dw1 = 0.0d0 

            denergy_db = 0.0d0 
            denergy_dw = 0.0d0 

            dforce_db = 0.0d0 
            dforce_dw = 0.0d0 
            
            do iatom = 1, natoms_step(ns)
               force_pred(:, iatom) = 0.0d0
            end do
         else if(net_type.eq.'cnet') then
            do iatom = 1, natoms_step(ns)
               do itype = 1,natom_list
                  do m1 = 1, mnet_size(1)
                     dcharge_db1(m1, itype, iatom) = 0.0d0
                     do m2 = 1,mnet_size(2)
                        dcharge_dw2(m2, m1,itype, iatom) = 0.0d0
                     end do
                  end do
               end do
            end do
            
            dcharge_db2 = 0.0d0 
            dcharge_dw1 = 0.0d0 
         endif
            
         do iatom = 1,natoms_step(ns)
            
            natoms0 = nbonded(iatom,ns) + 1
            j = 0 
            do k = 0,nbonded(iatom,ns)
               iat = atlist(k,iatom,ns)
               j = j + 1
               atname(j) = atname_step(iat,ns)
               atindex(j) = iat
            end do! do k = 0,nbonded(iatom,ns)
            
            call sort_string_list(atname,atindex,atname_sorted,atindex_sorted,&
                 atomid_list,range,nrange,natoms0)
!!$            write(*,*) 'sorted: ',(atname_sorted(j),j=1,natoms0)
!!$            write(*,*) 'index: ',(atindex_sorted(j),j=1,natoms0)
!!$            write(*,*) 'atomid_list',(atomid_list(j),j=1,nrange)
!!$            write(*,*) 'range',(range(j),j=0,nrange)
!!$            write(*,*) 'nrange ',nrange

            call get_atom_id(atname_sorted(1),atomid_central)
            
         npair = 0 
         n_types = 0
         n_types0 = 0 
         npairs_of_type = 0 
         ndis = 0 
         ndis_central = 0 
         ndis_non_central = 0 

         do n1 = 1,nrange
            n1_start = range(n1-1)
            n1_end = range(n1) - 1
            do n2 = n1,nrange
               n2_start = range(n2-1)
               n2_end = range(n2) - 1
               
               if(n1.eq.n2) then 
                  !same type
                  !check to see if there is more than zero pairs
                  if(n1_end-n1_start.gt.0) then 
                     n_types = n_types + 1
                     ipair = 0 
                     do k1 = 1,n1_end - n1_start + 1
                        do k2 = k1+1,n1_end - n1_start + 1
                           npair = npair + 1
                           ipair = ipair + 1
                           
                           iat = atindex_sorted(n1_start + k1 - 1)
                           jat = atindex_sorted(n2_start + k2 - 1)

                           xdif = coord(1,jat,ns) - coord(1,iat,ns)
                           ydif = coord(2,jat,ns) - coord(2,iat,ns)
                           zdif = coord(3,jat,ns) - coord(3,iat,ns)
                           rdis2 = xdif**2 + ydif**2 + zdif**2 

                           inside_cutoff = .true.
                           if(rdis2.gt.arc_cut2) inside_cutoff = .false.
                           if(inside_cutoff) ndis_non_central = ndis_non_central + 1

                           if(inside_cutoff) then 
                              rdis = dsqrt(rdis2)
                              
                              ndis = ndis + 1
                              npairs_of_type(n_types) = npairs_of_type(n_types) + 1
                              
                              pair_atno(1,ndis) = iat
                              pair_atno(2,ndis) = jat
                              xyz(0,ndis) = rdis
                              xyz(1,ndis) = xdif
                              xyz(2,ndis) = ydif
                              xyz(3,ndis) = zdif
                              xyz(4,ndis) = 1.0d0/rdis
                              xyz(5,ndis) = dexp(alpha*(rdis - beta))
                              xyz(6,ndis) = dexp(alpha*(rdis - beta_arc))
                           endif
                        end do
                     end do
                  endif
               else
!                 i.e. n1 .ne. n2
                  n_types = n_types + 1
                  if(n1.eq.1) n_types0 = n_types0 + 1
                  
                  ipair = 0 
                  do k1 = 1,n1_end - n1_start + 1
                     do k2 = 1,n2_end - n2_start + 1
                        npair = npair + 1
                        ipair = ipair + 1

                        iat = atindex_sorted(n1_start+k1 - 1)
                        jat = atindex_sorted(n2_start +k2 - 1)
                        
                        xdif = coord(1,jat,ns) - coord(1,iat,ns)
                        ydif = coord(2,jat,ns) - coord(2,iat,ns)
                        zdif = coord(3,jat,ns) - coord(3,iat,ns)
                        rdis2 = xdif**2 + ydif**2 + zdif**2

                        inside_cutoff = .true.
                        if(n1.eq.1) then 
                           !central atom
                           ndis_central = ndis_central + 1
                        else
                           if(rdis2.gt.arc_cut2) inside_cutoff = .false.
                           if(inside_cutoff) ndis_non_central = ndis_non_central + 1
                        endif
                        if(inside_cutoff) then 
                           rdis = dsqrt(rdis2)
                           ndis = ndis + 1
                           npairs_of_type(n_types) = npairs_of_type(n_types) + 1
                           
                           pair_atno(1,ndis) = iat
                           pair_atno(2,ndis) = jat
                           xyz(0,ndis) = rdis
                           xyz(1,ndis) = xdif
                           xyz(2,ndis) = ydif
                           xyz(3,ndis) = zdif
                           xyz(4,ndis) = 1.0d0/rdis
                           xyz(5,ndis) = dexp(alpha*(rdis - beta))
                           xyz(6,ndis) = dexp(alpha*(rdis - beta_arc))
                        endif
                     end do
                  end do
               endif
               
            end do
         end do

         ndis2 = ndis
         if(use_sumtype) ndis2 = ndis2 + ndis_non_central
         if(use_qrtype.and.net_type.eq.'mnet') ndis2 = ndis2 + ndis_central

         allocate(eigentype(ndis2))
         allocate(dis_type(ndis2))
         allocate(pair_index(ndis2))                           
         allocate(kplus_index(ndis2))                           
         allocate(kminus_index(ndis2))                           

!make eigentype, which indicates whether the interaction is a distance, (eigentype = 0) or a sumtype (eigentype = 1)         

         nd = 0
         k = 0 
         do itype0 = 1,n_types
            kplus = 0 
            kminus = 0 
            do ipair = 1,npairs_of_type(itype0)
               nd = nd + 1
               k = k + 1
               eigentype(nd) = 0
               dis_type(nd) = itype0
               pair_index(nd) = k
            end do

            if(itype0.le.n_types0.and.net_type.eq.'mnet') then 
               !do the qrplus terms

               k = k - npairs_of_type(itype0)
               do ipair = 1,npairs_of_type(itype0)
                  k = k + 1
                  
                  iat = pair_atno(1,k)
                  jat = pair_atno(2,k)
                  qj = charge(jat,ns)
                  if(qj.eq.0) qj = 1.d-10

                  if(qj.ge.0) then 
                     nd = nd + 1
                     kplus = kplus + 1
                     eigentype(nd) = 2
                     dis_type(nd) = itype0
                     pair_index(nd) = k
                     kplus_index(nd) = kplus
!                     write(*,*) 'kkknd',nd,kplus
                  endif
               end do

               !do the qrminus terms
               k = k - npairs_of_type(itype0)
               do ipair = 1,npairs_of_type(itype0)
                  k = k + 1
                  
                  iat = pair_atno(1,k)
                  jat = pair_atno(2,k)
                  qj = charge(jat,ns)
                  if(qj.eq.0) qj = 1.d-10
!                  write(*,*) 'qq1',qj,itype0

                  if(qj.lt.0) then 
                     nd = nd + 1
                     kminus = kminus + 1
                     eigentype(nd) = 3
                     dis_type(nd) = itype0
                     pair_index(nd) = k
                     kminus_index(nd) = kminus
                  endif
               end do
            endif

            if(use_sumtype.and.itype0.gt.n_types0) then 
               k = k - npairs_of_type(itype0)
               do ipair = 1,npairs_of_type(itype0)
                  nd = nd + 1
                  k = k + 1
                  eigentype(nd) = 1
                  dis_type(nd) = itype0 
                  pair_index(nd) = k
               end do
            endif
         end do

         if(net_type.eq.'mnet') then 
            ndis3 = ndis2 + 1
         else
            ndis3 = ndis2
         endif

         allocate(x1vec(ndis3))
         
         allocate(dy_dx1(ndis3))
         
         allocate(zz(ndis3))
         allocate(size(0:mnet_nlayers))

         size(0) = ndis3
         do m = 1,mnet_nlayers
            size(m) = mnet_size(m)
         end do

         call model%init(mnet_nlayers,size,calc_forces)

         !put in the new w1 here.

!         write(*,*) 'central',atomid_central

         npair = 0 
         n_types = 0 
         do n1 = 1,nrange
            n1_start = range(n1-1)
            n1_end = range(n1) - 1
            iatomid = atomid_list(n1)
            do n2 = n1,nrange
               jatomid = atomid_list(n2)
               if(n1.eq.n2.and.n1_start.eq.n1_end) cycle
               n_types = n_types + 1
               
               do k = 1,npairs_of_type(n_types)
                  npair = npair + 1
                  maxid = max(iatomid,jatomid)
                  minid = min(iatomid,jatomid)
                  pairlist(0,npair) = k
                  pairlist(1,npair) = maxid
                  pairlist(2,npair) = minid
               end do

            end do
         end do

         do m1 = 1,mnet_size(1)
            if(net_type.eq.'mnet') then 
            else if(net_type.eq.'cnet') then 
               !work on later!
!               w2(m1) = cnet_wmat2(m1,atomid_central)
!               b1vec(m1) = cnet_b1vec(m1,atomid_central)
            endif
         end do

         
         if(net_type.eq.'mnet') then 
            do layer = 1,mnet_nlayers
               if(layer.gt.1) then 
                  do m1 = 1,mnet_size(layer-1)
                     do m2 = 1,mnet_size(layer)
                        model%weight(layer)%mat2(m2,m1) = mnet_w(m2,m1,atomid_central,layer)
                     end do
                  end do
               endif
               if(layer.lt.mnet_nlayers) then 
                  do m1 = 1,mnet_size(layer)
                     model%bias(layer)%vec(m1) = mnet_b(m1,atomid_central,layer) 
                  end do 
                  model%bias(mnet_nlayers)%vec = 0.0d0 
               endif
            end do
         else if(net_type.eq.'cnet') then 
         endif
            
         if(net_type.eq.'mnet') then 
            do m0 = 1,ndis3
               if(m0.le.ndis2) then 
                  p = pair_index(m0)
                  k = pairlist(0,p)
                  maxid = pairlist(1,p)
                  minid = pairlist(2,p)
                  etype = eigentype(m0)
                  
                  kk = k
                  if(etype.eq.2) kk = kplus_index(m0)
                  if(etype.eq.3) kk = kminus_index(m0)

                  cindex = col_index(etype,kk,maxid,minid,atomid_central)
               else
                  cindex = charge_col_index(atomid_central)
               endif
               do m1 = 1,mnet_size(1)
                  model%weight(1)%mat2(m1,m0) = mnet_wmat1(m1,cindex)
               end do
            end do
         else if(net_type.eq.'cnet') then 

               do m0 = 1,ndis2
                  p = pair_index(m1)
                  k = pairlist(0,p)
                  maxid = pairlist(1,p)
                  minid = pairlist(2,p)
                  etype = eigentype(m0)
                  cindex = col_index(etype,k,maxid,minid,atomid_central)
                  cindex_atom = col_atom_index(etype,k,maxid,minid,atomid_central)
!                  if(cindex_atom.eq.1) write(*,*) 'wmat1',cnet_wmat1(m2,cindex),m2,cindex,&
!        'x',etype,k,maxid,minid,atomid_central

                  do m1 = 1,mnet_size(1)
                     w1(m1,m0) = cnet_wmat1(m1,cindex)
                  end do
               end do
            endif

         if(net_type.eq.'mnet') then 
            b2vec = 0.0d0 
            model%bias(mnet_nlayers)%vec(:) = 0.0d0 
         else if(net_type.eq.'cnet') then 
            b2vec = cnet_b2vec(atomid_central)            
            model%bias(cnet_nlayers)%vec(:) = cnet_b2vec(atomid_central)
         endif
         
         ndis0 = 0 
         ndis1 = 0 
         nev_data = 0 


         nev_data = 0 
         do itype0 = 1,n_types
!            !factor of 2 if we include qr types
            nev_data = nev_data + 2*npairs_of_type(itype0)**2
            if(use_sumtype.and.itype0.gt.n_types0) then
               nev_data = nev_data + npairs_of_type(itype0)**2
            endif
         end do
         
         allocate(evector_data(nev_data))


         nev_data = 0 
         do itype0 = 1,n_types
            allocate(x1vec_0(npairs_of_type(itype0)))
            allocate(x1vec_1(npairs_of_type(itype0)))            
            allocate(x1vec_qr(npairs_of_type(itype0)))            
            allocate(evalue0(npairs_of_type(itype0)))
            allocate(evalue1(npairs_of_type(itype0)))            
            allocate(evalue_qr(npairs_of_type(itype0)))
            allocate(evector0(npairs_of_type(itype0),npairs_of_type(itype0)))
            allocate(evector1(npairs_of_type(itype0),npairs_of_type(itype0)))            
            allocate(evector_qr(npairs_of_type(itype0),npairs_of_type(itype0)))            

            !count nqr_plus,nqr_minus

            nplus = 0 
            nminus = 0 
            if(itype0.le.n_types0) then 
               nd = ndis0
               do ipair = 1,npairs_of_type(itype0)
                  nd = nd + 1
                  iat = pair_atno(1,nd)
                  jat = pair_atno(2,nd)
                  qj = charge(jat,ns)
                  if(qj.eq.0) qj = 1.d-10
!                  write(*,*) 'qq2',qj
                  if(qj.ge.0) then 
                     nplus = nplus + 1
                  else
                     nminus = nminus + 1
                  endif
               end do
!               write(*,*) 'nplus',nplus,'nminus',nminus,'itype',itype0
            endif

            allocate(x1vec_qrplus(nplus))
            allocate(x1vec_qrminus(nminus))
            allocate(evector_qrplus(nplus,nplus))
            allocate(evector_qrminus(nminus,nminus))
            allocate(evalue_qrplus(nplus))
            allocate(evalue_qrminus(nminus))

            kminus = 0 
            kplus = 0 
            do ipair = 1,npairs_of_type(itype0)
               ndis0 = ndis0 + 1

               rdis = xyz(0,ndis0)
               x1vec_0(ipair) = rdis
               iat = pair_atno(1,ndis0)
               jat = pair_atno(2,ndis0)

               !if this pair does not involve the central atom

               if(iat.ne.iatom.and.jat.ne.iatom) then
                  
               !find the corresponding pairs connecting to the central atom               

                  do n0 = 1,ndis_central
                     iat0 = pair_atno(1,n0)
                     if(iat0.eq.iatom) then 
                        jat0 = pair_atno(2,n0)
                        if(jat0.eq.iat) nn1 = n0
                        if(jat0.eq.jat) nn2 = n0 
                     endif
                  end do

                  rdis_1 = xyz(0,nn1)
                  efac_1 = xyz(5,nn1)

                  rdis_2 = xyz(0,nn2)                  
                  efac_2 = xyz(5,nn2)
                  
                  efac_arc = xyz(6,ndis0)

                  x1vec_0(ipair) = x1vec_0(ipair) + efac_1 + efac_2 + efac_arc
                  !factor of 2, because the sum distances pick up exponentials from both
                  !the arc bond and the two radial bonds
                  x1vec_1(ipair) = rdis_1 + rdis_2 + rdis + 2.0d0*efac_1 + 2.0d0*efac_2 + efac_arc
               else
                  rdis = xyz(0,ndis0)
                  efac = xyz(5,ndis0)
                  x1vec_0(ipair) = x1vec_0(ipair) + efac
                  qj = charge(jat,ns)
                  if(qj.eq.0) qj = 1.d-10

                  if(qj.ge.0) then 
                     kplus = kplus + 1
                     x1vec_qrplus(kplus) = (rdis + efac)/qj
                  else
                     kminus = kminus + 1
                     x1vec_qrminus(kminus) = (rdis + efac)/qj
                  endif
               endif

            end do !do ipair = 1,npairs_of_type(itype0)

            if(npairs_of_type(itype0).gt.0) then 
               sort_direction = "ascending"
               call eigen(x1vec_0,couple,npairs_of_type(itype0),evalue0,evector0,sort_direction)
            endif            
            do i1 = 1,npairs_of_type(itype0)
               ndis1 = ndis1 + 1
               zz(ndis1) = 1.0d0/evalue0(i1)
            end do

            do i1 = 1,npairs_of_type(itype0)
               do i2 = 1,npairs_of_type(itype0)
                  nev_data = nev_data + 1
                  evector_data(nev_data) = evector0(i1,i2)
               end do
            end do
            


            !now do the same for qr. (Maybe simplify loop later.)

            !note: npairs_of_type(itype0) should equal nplus + nminus

            if(itype0.le.n_types0) then 
               if(nplus.gt.0) then 
                  sort_direction = "ascending"
                  call eigen(x1vec_qrplus,couple,nplus,evalue_qrplus,evector_qrplus,sort_direction)
!                  write(*,*) 'plus',1.0d0/evalue_qrplus
                  do i1 = 1,nplus
                     ndis1 = ndis1 + 1
                     zz(ndis1) = 1.0d0/evalue_qrplus(i1)
                  end do
                  do i1 = 1,nplus
                     do i2 = 1,nplus
                        nev_data = nev_data + 1
                        evector_data(nev_data) = evector_qrplus(i1,i2)
                     end do
                  end do
               endif
               if(nminus.gt.0) then 
                  sort_direction = "descending"
                  call eigen(x1vec_qrminus,couple,nminus,evalue_qrminus,evector_qrminus,sort_direction)
!                  write(*,*) 'minus',1.0d0/evalue_qrminus
                  do i1 = 1,nminus
                     ndis1 = ndis1 + 1
                     zz(ndis1) = 1.0d0/evalue_qrminus(i1)
                  end do
                  do i1 = 1,nminus
                     do i2 = 1,nminus
                        nev_data = nev_data + 1
                        evector_data(nev_data) = evector_qrminus(i1,i2)
                     end do
                  end do
               endif
            endif

            if(use_sumtype.and.itype0.gt.n_types0) then

               if(npairs_of_type(itype0).gt.0) then 
                  sort_direction = "ascending"
                  call eigen(x1vec_1,couple,npairs_of_type(itype0),evalue1,evector1,sort_direction)
               endif
               
               do i1 = 1,npairs_of_type(itype0)
                  ndis1 = ndis1 + 1
                  zz(ndis1) = 1.0d0/evalue1(i1)
               end do
            
               do i1 = 1,npairs_of_type(itype0)
                  do i2 = 1,npairs_of_type(itype0)
                     nev_data = nev_data + 1
                     evector_data(nev_data) = evector1(i1,i2)
                  end do
               end do
            endif


            deallocate(x1vec_0)
            deallocate(x1vec_1)            
            deallocate(x1vec_qr)            
            deallocate(evalue0)
            deallocate(evector0)
            deallocate(evalue_qr)
            deallocate(evalue1)
            deallocate(evector1)
            deallocate(evector_qr)            
            deallocate(x1vec_qrplus)
            deallocate(x1vec_qrminus)
            deallocate(evector_qrplus)
            deallocate(evector_qrminus)
            deallocate(evalue_qrplus)
            deallocate(evalue_qrminus)

         end do

         ! include the charge
         if(net_type.eq.'mnet') then
            ndis1 = ndis1 + 1
            zz(ndis1) = charge(iatom,ns)
         endif


         model%x0(:) = zz(:)

         !now run the neural network!!!
         call cpu_time(start)
         call model%run_net
         call cpu_time(end) 
         net_timing =  net_timing + end - start
         yvec = model%y

!         write(*,*) 'yvec = ',yvec,ns,iatom

         if(net_type.eq.'mnet') then 
            energy_pred = energy_pred + yvec
         else if(net_type.eq.'cnet') then 
            charge_pred(iatom) = yvec
         endif
         
         if(net_type.eq.'mnet') then 
            do m1 = 1,ndis3
               if(m1.le.ndis2) then 
                  p = pair_index(m1)
                  k = pairlist(0,p)
                  maxid = pairlist(1,p)
                  minid = pairlist(2,p)
                  etype = eigentype(m1)

                  kk = k 
                  if(etype.eq.2) kk = kplus_index(m1)
                  if(etype.eq.3) kk = kminus_index(m1)
                  cindex = col_index(etype,kk,maxid,minid,atomid_central)
               else
                  cindex = charge_col_index(atomid_central) 
               endif
               do m2 = 1, mnet_size(1)
                  denergy_dw1(m2,cindex) = denergy_dw1(m2,cindex) + model%dy_dweight(1)%mat2(m2,m1)
               end do
            end do
         else if(net_type.eq.'cnet') then 
            do m1 = 1,ndis3
               p = pair_index(m1)
               k = pairlist(0,p)
               maxid = pairlist(1,p)
               minid = pairlist(2,p)
               etype = eigentype(m1)
               do m2 = 1, mnet_size(1)
                  cindex_atom = col_atom_index(etype,k,maxid,minid,atomid_central)
                  dcharge_dw1(m2,cindex_atom,atomid_central,iatom) = & 
                       dcharge_dw1(m2,cindex_atom,atomid_central,iatom) + model%dy_dweight(1)%mat2(m2,m1)
               end do
            end do
         endif

         if(net_type.eq.'mnet') then 

            do layer = 1,mnet_nlayers
               if(layer.lt.mnet_nlayers) then 
                  do m1 = 1,mnet_size(layer)
                     denergy_db(m1,atomid_central,layer) = denergy_db(m1,atomid_central,layer) + model%dy_dbias(layer)%vec(m1)
                  end do
               endif
               if(layer.gt.1) then 
                  do m1 = 1,mnet_size(layer-1)
                     do m2 = 1,mnet_size(layer)
                        denergy_dw(m2,m1,atomid_central,layer) = denergy_dw(m2,m1,atomid_central,layer) + &
                             model%dy_dweight(layer)%mat2(m2,m1)
                     end do
                  end do
               endif
            end do



         else if(net_type.eq.'cnet') then 
            do m1 = 1,mnet_size(1)
               dcharge_db1(m1,atomid_central,iatom) = dcharge_db1(m1,atomid_central,iatom) + model%dy_dbias(1)%vec(m1)
               do m2 = 1,mnet_size(2)
                  dcharge_dw2(m2,m1,atomid_central,iatom) = dcharge_dw2(m2,m1,atomid_central,iatom) + & 
                       model%dy_dweight(2)%mat2(m2,m1)
               end do 
            end do
            do m2 = 1,mnet_size(2)
               dcharge_db2(m2,atomid_central,iatom) = dcharge_db1(m2,atomid_central,iatom) + model%dy_dbias(2)%vec(m1)
               do m3 = 1,mnet_size(3)
                  dcharge_dw2(m3,m2,atomid_central,iatom) = dcharge_dw2(m3,m2,atomid_central,iatom) + & 
                       model%dy_dweight(3)%mat2(m2,m1)
               end do 
            end do



            charge_pred(iatom) = charge_pred(iatom) + b2vec
            

         endif

         if(net_type.eq.'mnet') then 
         
            !convert from zz derivatives to x derivatives
            dy_dx1 = 0.0d0

!            d2y_dx1dw1 = 0.0d0 
!            d2y_dx1dw2 = 0.0d0 
!            d2y_dx1db1 = 0.0d0 
         
            !now use the hellmann-feynman theorem to get the right derivatives. 

            ev_index = 0 
            do nd1 = 1,ndis2
               itype1 = dis_type(nd1)
               etype1 = eigentype(nd1)
               do nd2 = 1,ndis2
                  itype2 = dis_type(nd2)
                  etype2 = eigentype(nd2)
                  if(itype1.eq.itype2.and.etype1.eq.etype2) then
                     
                     ev_index = ev_index + 1
                     ev2 = evector_data(ev_index)**2

                     d_evaluei_d_evalue = -zz(nd2)**2

                     dy_dx1(nd1) = dy_dx1(nd1) + model%dy_dx(nd2)*ev2 * d_evaluei_d_evalue

                     do layer = 1,mnet_nlayers

                        if(layer.lt.mnet_nlayers) then 
                           do m1 = 1,size(layer)
                              model%d2y_dxdb_(nd1)%vecvec(layer)%vec(m1) = &
                                   model%d2y_dxdb_(nd1)%vecvec(layer)%vec(m1) + & 
                                   model%d2y_dxdb(nd2)%vecvec(layer)%vec(m1)*ev2 * d_evaluei_d_evalue
                           end do
                        endif

                        do m1 = 1,size(layer-1)
                           do m2 = 1,size(layer)
                              model%d2y_dxdw_(nd1)%vecmat2(layer)%mat2(m2,m1) = &
                                   model%d2y_dxdw_(nd1)%vecmat2(layer)%mat2(m2,m1) + & 
                                   model%d2y_dxdw(nd2)%vecmat2(layer)%mat2(m2,m1)*ev2 * d_evaluei_d_evalue
                           end do
                        end do

                     end do 
                     
                  endif
               end do
            end do
         endif
            
         deallocate(evector_data)

            !only need the predictive forces for mnet
         if(net_type.eq.'mnet'.and.minmode.ne.'energy') then         

            do n = 1,ndis2
               p = pair_index(n)
               etype = eigentype(n)
               
               iat = pair_atno(1,p)
               jat = pair_atno(2,p) 
               
               rdis = xyz(0,p)
               rr(1) = xyz(1,p)
               rr(2) = xyz(2,p)
               rr(3) = xyz(3,p)
               rdisi = xyz(4,p)
               efac = xyz(5,p)

               aa0 = rdisi
               if(etype.eq.2.or.etype.eq.3) then
                  !iat should be the central charge, so need charge on jat
                  qj = charge(jat,ns)
                  if(qj.eq.0) qj = 1.d-10
                  
                  !multiply by reciprocal charge
                  aa0 = aa0/qj
               endif

               force_pred(:,iat) = force_pred(:,iat) + dy_dx1(n) * rr(:) * aa0
               force_pred(:,jat) = force_pred(:,jat) - dy_dx1(n) * rr(:) * aa0

               
               !if this pair does not involve the central atom
               if(iat.ne.iatom.and.jat.ne.iatom) then 
                  
                     efac_arc = xyz(6,p)
                     if(efac_arc.gt.ecut) then 
                        aa_arc = dy_dx1(n) * alpha * efac_arc * rdisi
                        
                        force_pred(:,jat) = force_pred(:,jat) - aa_arc * rr(:)
                        force_pred(:,iat) = force_pred(:,iat) + aa_arc * rr(:)
                     endif
                     
                  !find the corresponding pairs connecting to the central atom                              
                  do n0 = 1,ndis_central
                     iat0 = pair_atno(1,n0)
                     if(iat0.eq.iatom) then 
                        jat0 = pair_atno(2,n0)
                        if(jat0.eq.iat) nn1 = n0
                        if(jat0.eq.jat) nn2 = n0 
                     endif
                  end do
                  
                  efac_1 = xyz(5,nn1)
                  efac_2 = xyz(5,nn2)

                  if((etype.eq.0.and.efac_1.gt.ecut).or.etype.eq.1) then 
                     rdis_1 = xyz(0,nn1)
                     rr_1(1) = xyz(1,nn1)
                     rr_1(2) = xyz(2,nn1)
                     rr_1(3) = xyz(3,nn1)
                     rdis_1i = xyz(4,nn1)

                     iat1 = pair_atno(1,nn1)
                     jat1 = pair_atno(2,nn1) 

                     
                     if(etype.eq.0) then
                        aa1 = dy_dx1(n) * alpha * efac_1 * rdis_1i
                     else if(etype.eq.1) then 
                        aa1 = dy_dx1(n) * rdis_1i * (1.0d0 + 2.0d0 * alpha * efac_1)
                     endif
                     
                     force_pred(:,jat1) = force_pred(:,jat1) - aa1 * rr_1(:)
                     force_pred(:,iat1) = force_pred(:,iat1) + aa1 * rr_1(:)
                     
                  endif
                  
                  if((etype.eq.0.and.efac_2.gt.ecut).or.etype.eq.1) then 
                     rdis_2 = xyz(0,nn2)                  
                     rr_2(1) = xyz(1,nn2)
                     rr_2(2) = xyz(2,nn2)
                     rr_2(3) = xyz(3,nn2)
                     rdis_2i = xyz(4,nn2)
                     
                     iat2 = pair_atno(1,nn2)
                     jat2 = pair_atno(2,nn2) 


                     if(etype.eq.0) then
                        aa2 = dy_dx1(n) * alpha * efac_2 * rdis_2i
                     else if(etype.eq.1) then 
                        aa2 = dy_dx1(n) * rdis_2i* (1.0d0 + 2.0d0 * alpha * efac_2)
                     endif

                     force_pred(:,jat2) = force_pred(:,jat2) - aa2 * rr_2(:)
                     force_pred(:,iat2) = force_pred(:,iat2) + aa2 * rr_2(:)

                  endif
                  
               else if(efac.gt.ecut) then 
                  
                  aa = dy_dx1(n) * alpha * efac * rdisi
                  if(etype.eq.2.or.etype.eq.3) aa = aa/qj                  
                  
                  force_pred(:,jat) = force_pred(:,jat) - aa * rr(:)
                  force_pred(:,iat) = force_pred(:,iat) + aa * rr(:)

               endif
               
               if(calc_forces) then 
               do m1 = 1,mnet_size(1)
                  do m0 = 1,ndis3
                     aa = model%d2y_dxdw_(n)%vecmat2(1)%mat2(m1,m0)*rdisi

                     if(m0.le.ndis2) then 
                        p = pair_index(m0)
                        
                        k = pairlist(0,p)
                        maxid = pairlist(1,p)
                        minid = pairlist(2,p)
                        etype2 = eigentype(m0)

                        kk = k 
                        if(etype2.eq.2) kk = kplus_index(m0)
                        if(etype2.eq.3) kk = kminus_index(m0)

                        cindex = col_index(etype2,kk,maxid,minid,atomid_central)
                     else
                        cindex = charge_col_index(atomid_central)                         
                     endif
                        
                     if(etype.eq.2.or.etype.eq.3) aa = aa/qj

                     dforce_dw1(:,cindex,m1,iat) = dforce_dw1(:,cindex,m1,iat) + aa * rr(:)
                     dforce_dw1(:,cindex,m1,jat) = dforce_dw1(:,cindex,m1,jat) - aa * rr(:)
                        
                     if(iat.ne.iatom.and.jat.ne.iatom) then

                        if(efac_arc.gt.ecut) then 
                           aa_arc = model%d2y_dxdw_(n)%vecmat2(1)%mat2(m1,m0)*alpha*efac_arc*rdisi

                           dforce_dw1(:,cindex,m1,jat) = dforce_dw1(:,cindex,m1,jat) - aa_arc * rr(:)
                           dforce_dw1(:,cindex,m1,iat) = dforce_dw1(:,cindex,m1,iat) + aa_arc * rr(:)

                        endif

                        if((etype.eq.0.and.efac_1.gt.ecut).or.etype.eq.1) then 

                           if(etype.eq.0) then
                              aa1 = model%d2y_dxdw_(n)%vecmat2(1)%mat2(m1,m0) * alpha * efac_1 * rdis_1i
                           else if(etype.eq.1) then 
                              aa1 = model%d2y_dxdw_(n)%vecmat2(1)%mat2(m1,m0) * rdis_1i * (1.0d0 + 2.0d0 * alpha * efac_1)
                           endif

                           dforce_dw1(:,cindex,m1,jat1) = dforce_dw1(:,cindex,m1,jat1) - aa1 * rr_1(:)
                           dforce_dw1(:,cindex,m1,iat1) = dforce_dw1(:,cindex,m1,iat1) + aa1 * rr_1(:)

                        endif

                        if((etype.eq.0.and.efac_2.gt.ecut).or.etype.eq.1) then 

                           if(etype.eq.0) then
                              aa2 = model%d2y_dxdw_(n)%vecmat2(1)%mat2(m1,m0) * alpha * efac_2 * rdis_2i
                           else if(etype.eq.1) then 
                              aa2 = model%d2y_dxdw_(n)%vecmat2(1)%mat2(m1,m0) * rdis_2i * (1.0d0 + 2.0d0 * alpha * efac_2)
                           endif

                           dforce_dw1(:,cindex,m1,jat2) = dforce_dw1(:,cindex,m1,jat2) - aa2 * rr_2(:)
                           dforce_dw1(:,cindex,m1,iat2) = dforce_dw1(:,cindex,m1,iat2) + aa2 * rr_2(:)

                        endif


                     else if(efac.gt.ecut) then 
                        aa = model%d2y_dxdw_(n)%vecmat2(1)%mat2(m1,m0)*alpha*rdisi*efac
                        if(etype.eq.2.or.etype.eq.3) aa = aa/qj                        

                        dforce_dw1(:,cindex,m1,jat) = dforce_dw1(:,cindex,m1,jat) - aa * rr(:)
                        dforce_dw1(:,cindex,m1,iat) = dforce_dw1(:,cindex,m1,iat) + aa * rr(:)

                     endif
                  end do !do m0 1,ndis3
               end do

               !now do the other arrays


               do layer = 1,mnet_nlayers

                  if(layer.gt.1) then 
                  do m1 = 1,mnet_size(layer-1)
                        do m2 = 1,mnet_size(layer)
                     
                           aa = model%d2y_dxdw_(n)%vecmat2(layer)%mat2(m2,m1)*rdisi
                           dforce_dw(:,m2,m1,iat,atomid_central,layer) = &
                                dforce_dw(:,m2,m1,iat,atomid_central,layer) + aa * rr(:)
                           dforce_dw(:,m2,m1,jat,atomid_central,layer) = &
                                dforce_dw(:,m2,m1,jat,atomid_central,layer) - aa * rr(:)
                           
                           if(iat.ne.iatom.and.jat.ne.iatom) then
                              
                              if(efac_arc.gt.ecut) then 
                                 aa_arc = model%d2y_dxdw_(n)%vecmat2(layer)%mat2(m2,m1)*alpha*efac_arc*rdisi
                                 dforce_dw(:,m2,m1,jat,atomid_central,layer) = &
                                      dforce_dw(:,m2,m1,jat,atomid_central,layer) - aa_arc * rr(:)
                                 dforce_dw(:,m2,m1,iat,atomid_central,layer) = & 
                                      dforce_dw(:,m2,m1,iat,atomid_central,layer) + aa_arc * rr(:)
                              endif
                        
                              if((etype.eq.0.and.efac_1.gt.ecut).or.etype.eq.1) then 
                                 
                                 if(etype.eq.0) then
                                    aa1 = model%d2y_dxdw_(n)%vecmat2(layer)%mat2(m2,m1) * alpha * efac_1 * rdis_1i
                                 else if(etype.eq.1) then 
                                    aa1 = model%d2y_dxdw_(n)%vecmat2(layer)%mat2(m2,m1) * rdis_1i * (1.0d0 + 2.0d0 * alpha * efac_1)
                                 endif
                                 
                                 dforce_dw(:,m2,m1,jat1,atomid_central,layer) = &
                                      dforce_dw(:,m2,m1,jat1,atomid_central,layer) - aa1 * rr_1(:)
                                 dforce_dw(:,m2,m1,iat1,atomid_central,layer) = &
                                      dforce_dw(:,m2,m1,iat1,atomid_central,layer) + aa1 * rr_1(:)
                                 
                              endif
                              if((etype.eq.0.and.efac_2.gt.ecut).or.etype.eq.1) then 
                                 
                                 if(etype.eq.0) then
                                    aa2 = model%d2y_dxdw_(n)%vecmat2(layer)%mat2(m2,m1) * alpha * efac_2 * rdis_2i
                                 else if(etype.eq.1) then 
                                    aa2 = model%d2y_dxdw_(n)%vecmat2(layer)%mat2(m2,m1) * rdis_2i * (1.0d0 + 2.0d0 * alpha * efac_2)
                                 endif
                                 
                                 dforce_dw(:,m2,m1,jat2,atomid_central,layer) = &
                                      dforce_dw(:,m2,m1,jat2,atomid_central,layer) - aa2 * rr_2(:)
                                 dforce_dw(:,m2,m1,iat2,atomid_central,layer) = &
                                      dforce_dw(:,m2,m1,iat2,atomid_central,layer) + aa2 * rr_2(:)
                                 
                              endif
                              
                           else if(efac.gt.ecut) then 
                              
                              aa = model%d2y_dxdw_(n)%vecmat2(layer)%mat2(m2,m1)*rdisi * alpha * efac
                              
                              dforce_dw(:,m2,m1,jat,atomid_central,layer) = dforce_dw(:,m2,m1,jat,atomid_central,layer) - aa * rr(:)
                              dforce_dw(:,m2,m1,iat,atomid_central,layer) = dforce_dw(:,m2,m1,iat,atomid_central,layer) + aa * rr(:)
                              
                           endif
                        end do !do m2 = 1,mnet_size(layer)
                     end do !do m1 = 1,mnet_size(layer-1)
                     endif

                     if(layer.lt.mnet_nlayers) then 
                        do m1 = 1,mnet_size(layer)

                           aa = model%d2y_dxdb_(n)%vecvec(layer)%vec(m1) * rdisi
                           
                           dforce_db(:,m1,iat,atomid_central,layer) = dforce_db(:,m1,iat,atomid_central,layer) + aa * rr(:)
                           dforce_db(:,m1,jat,atomid_central,layer) = dforce_db(:,m1,jat,atomid_central,layer) - aa * rr(:)
                        
                     if(iat.ne.iatom.and.jat.ne.iatom) then
                        
                        if(efac_arc.gt.ecut) then 
                           aa_arc = model%d2y_dxdb_(n)%vecvec(layer)%vec(m1) * alpha * efac_arc * rdisi
                           
                           dforce_db(:,m1,jat,atomid_central,layer) = dforce_db(:,m1,jat,atomid_central,layer) - aa_arc * rr(:)
                           dforce_db(:,m1,iat,atomid_central,layer) = dforce_db(:,m1,iat,atomid_central,layer) + aa_arc * rr(:)
                           
                        endif
                        
                        if((etype.eq.0.and.efac_1.gt.ecut).or.etype.eq.1) then 
                           
                           if(etype.eq.0) then
                              aa1 = model%d2y_dxdb_(n)%vecvec(layer)%vec(m1) * alpha * efac_1 * rdis_1i
                           else if(etype.eq.1) then 
                              aa1 = model%d2y_dxdb_(n)%vecvec(layer)%vec(m1) * rdis_1i * (1.0d0 + 2.0d0 * alpha * efac_1)
                           endif
                     
                           dforce_db(:,m1,jat1,atomid_central,layer) = dforce_db(:,m1,jat1,atomid_central,layer) - aa1 * rr_1(:)
                           dforce_db(:,m1,iat1,atomid_central,layer) = dforce_db(:,m1,iat1,atomid_central,layer) + aa1 * rr_1(:)
                           
                        endif
                        if((etype.eq.0.and.efac_2.gt.ecut).or.etype.eq.1) then 
                           if(etype.eq.0) then
                              aa2 = model%d2y_dxdb_(n)%vecvec(layer)%vec(m1) * alpha * efac_2 * rdis_2i
                           else if(etype.eq.1) then 
                              aa2 = model%d2y_dxdb_(n)%vecvec(layer)%vec(m1) * rdis_2i * (1.0d0 + 2.0d0 * alpha * efac_2)
                           endif
                           
                           dforce_db(:,m1,jat2,atomid_central,layer) = dforce_db(:,m1,jat2,atomid_central,layer) - aa2 * rr_2(:)
                           dforce_db(:,m1,iat2,atomid_central,layer) = dforce_db(:,m1,iat2,atomid_central,layer) + aa2 * rr_2(:)
                        endif
                        
                     else if(efac.gt.ecut) then 
                        aa = model%d2y_dxdb_(n)%vecvec(layer)%vec(m1) * rdisi * alpha * efac
                        
                        dforce_db(:,m1,jat,atomid_central,layer) = dforce_db(:,m1,jat,atomid_central,layer) - aa * rr(:)
                        dforce_db(:,m1,iat,atomid_central,layer) = dforce_db(:,m1,iat,atomid_central,layer) + aa * rr(:)
                     endif
                     
                  end do !do m1 = 1,mnet_size(1)
                  endif
               end do !do layer = 1,mnet_nlayers
               
            endif !if(calc_forces)
            end do !do n = 1,ndis2
            
         endif
            
            deallocate(x1vec)
            
            deallocate(dy_dx1)
            
            deallocate(zz)
            deallocate(eigentype)
            deallocate(dis_type)
            deallocate(pair_index)                                    
            deallocate(kplus_index)                                    
            deallocate(kminus_index)                                    

            call model%dealloc()
            deallocate(size)

         end do !do iatom = 1,natoms_step(ns)


         !charge neutrality

         if(net_type.eq.'cnet') then 
            allocate(charge0_pred(natoms_max))
            
            charge_total = 0.0d0 
            charge0_pred = charge_pred

            !xxx warning, possibly inefficient
            dcharge0_dw1 = dcharge_dw1
            dcharge0_dw2 = dcharge_dw2
            dcharge0_db1 = dcharge_db1
            dcharge0_db2 = dcharge_db2

            do iatom = 1,natoms_step(ns)
               do jatom = 1,natoms_step(ns)
                  charge_pred(iatom) = charge_pred(iatom) - charge0_pred(jatom)/dble(natoms_step(ns))
 
                  dcharge_dw2(:,:,:,iatom) = dcharge_dw2(:,:,:,iatom) - dcharge0_dw2(:,:,:,jatom)/dble(natoms_step(ns))
                  dcharge_db1(:,:,iatom) = dcharge_db1(:,:,iatom) - dcharge0_db1(:,:,jatom)/dble(natoms_step(ns))
                  dcharge_db2(:,:,iatom) = dcharge_db2(:,:,iatom) - dcharge0_db2(:,:,jatom)/dble(natoms_step(ns))

                  do id1 = 0,natom_list
                     do id2 = id1,natom_list
                        do etype = 0,3
                           do id = 1,natom_list
                              
                              
                              if(col_used(etype,id2,id1,id)) then
                                 
                                 nc_used = ncol_used(etype,id2,id1,id)
                                 do m1 = 1,nc_used

                                    cindex_atom = col_atom_index(etype,m1,id2,id1,id)
                                    
                                    if(atomid_step(jatom,ns).eq.id) then 
                                       do m2 = 1,mnet_size(1)
                                          dcharge_dw1(m2,cindex_atom,id,iatom) = &
                                               dcharge_dw1(m2,cindex_atom,id,iatom) & 
                                               - dcharge0_dw1(m2,cindex_atom,id,jatom)/dble(natoms_step(ns))
                                       end do 
                                    endif

                                 end do
                                 
                              endif
                              
                                 
                           end do
                           
                        end do
                     end do
                  end do
                  
               end do
            end do
         endif

!         write(*,*) 'energy_pred = ',energy_pred,ns
       end subroutine get_mnet_energy

      subroutine mnetwork(n1, n2, w1, w2, x1vec, b1vec, b2vec, &
                          yvec, dy_dw1, dy_dw2, dy_dx1, dy_db1, dy_db2, &
                          d2y_dx1dw1, d2y_dx1dw2, d2y_dx1db1 )
        implicit none
        ! Dummy arguments
        integer n1, n2
        real(8), intent(in), dimension(n1, n2) :: w1
        real(8), intent(in), dimension(n2) :: w2
        real(8), intent(in), dimension(n1) :: x1vec
        real(8), intent(in), dimension(n2) :: b1vec
        real(8), intent(in) :: b2vec
        real(8), intent(out) :: yvec
        real(8), intent(out), dimension(:,:) :: dy_dw1
        real(8), intent(out), dimension(:) :: dy_dw2, dy_dx1, dy_db1
        real(8), intent(out), dimension(:,:,:) :: d2y_dx1dw1
        real(8), intent(out), dimension(:,:) :: d2y_dx1dw2, d2y_dx1db1
        real(8) :: dy_db2
        ! NOTE: Put allocatable above result in Segmentation fault
        
        ! Local variables
        integer j
        real(8), dimension(n2) :: x2vec, x2vec0
        real(8), dimension(n2) :: dsigma, d2sigma
        
        ! Forward pass
        x2vec0 = MATMUL(transpose(w1), x1vec) + b1vec
        call get_sigma(x2vec0, x2vec, dsigma, d2sigma)
        yvec = DOT_PRODUCT(w2, x2vec)
        ! yvec = yvec + b2vec
        
        ! Backward pass
        dy_dw2 = x2vec
        dy_db1 = w2 * dsigma
        dy_dw1 = MATMUL(SPREAD(x1vec,2,1), SPREAD(dy_db1,1,1))  ![n1,n2]=[n1,1] matmul [1,n2]
        dy_dx1 = MATMUL(w1, dy_db1)
        d2y_dx1dw2 = TRANSPOSE(w1 * SPREAD(dsigma,1,n1)) ! [n2,n1]=([n1,n2] * [n1,n2])^T
        d2sigma    = w2 * d2sigma  ! d2sigma now cache a temporary value
        d2y_dx1db1 = TRANSPOSE(w1 * SPREAD(d2sigma,1,n1)) ! [n2,n1]=([n1,n2] * [n1,n2])
        do j = 1, n1
           ! [n1,n2] = [n1,1] matmul ([1, n2])
           d2y_dx1dw1(:,:,j) = MATMUL(SPREAD(x1vec,2,1), SPREAD(w1(j,:)*d2sigma,1,1))
           d2y_dx1dw1(j,:,j) = d2y_dx1dw1(j, :, j) + dy_db1
        end do
        
!        dy_db2 = 1.0d0 
!        yvec = yvec + b2vec

      end subroutine mnetwork
      
      subroutine test_mnet
         implicit none
         integer itype, i, m, n,frag_id
         real(8) error2_0, error2_1,tiny, diff,analyt
         integer id1,id2,id_center,m1,m2,nc_used,cindex,etype,layer
         real(8) :: error2_energy_0,error2_force_0,error2_m_reg_0
         real(8) :: error2_energy_1,error2_force_1,error2_m_reg_1

         write(*, *) '*********************'
         write(*, *) '***   MTEST GRADIENTS'
         write(*, *) '*********************'
         call flush (6)

         tiny = 1d-4
         call run_train_mnet

         error2_energy_0 = error2_energy
         error2_force_0 = error2_force
         error2_m_reg_0 = error2_m_reg
         error2_0 = error2_force + error2_energy + error2_m_reg 

!!$         do itype = 1,natom_list
!!$            if(use_type(itype)) then 
!!$            cindex = charge_col_index(itype) 
!!$            do m2 = 1,mnet_s1
!!$               mnet_wmat1(m2,cindex) = mnet_wmat1(m2,cindex) + tiny
!!$               call run_train_mnet
!!$               
!!$               error2_1 = error2_force + error2_energy + error2_m_reg
!!$               diff = (error2_1 - error2_0)/tiny
!!$               
!!$               analyt = mnet_de2_dw1(m2,cindex)
!!$               write(*,*) 'w1diff charge',itype,m2,diff,analyt,diff/analyt
!!$!               write(888,*) 'w1diff charge',itype,m2,diff,analyt,diff/analyt
!!$!               call flush(888)
!!$               mnet_wmat1(m2,cindex) = mnet_wmat1(m2,cindex) - tiny                           
!!$               
!!$            end do
!!$            endif
!!$         end do
!!$         write(*,*)
!!$
!!$         do id_center = 1,natom_list
!!$            do id1 = 0,natom_list
!!$               do id2 = id1,natom_list
!!$                  do etype = 3,3
!!$                     if(etype.eq.1.and.(.not.use_sumtype)) exit
!!$                     if(col_used(etype,id2,id1,id_center)) then 
!!$                        nc_used = ncol_used(etype,id2,id1,id_center)
!!$                        do m1 = 1,nc_used
!!$                           cindex = col_index(etype,m1,id2,id1,id_center)
!!$                           do m2 = 1,mnet_s1
!!$                              mnet_wmat1(m2,cindex) = mnet_wmat1(m2,cindex) + tiny
!!$                              call run_train_mnet
!!$
!!$                              error2_1 = error2_force + error2_energy + error2_m_reg
!!$                              diff = (error2_1 - error2_0)/tiny
!!$
!!$                              analyt = mnet_de2_dw1(m2,cindex)
!!$                              write(*, *) 'w1diff',etype,m1,m2,id2,id1,id_center,diff,analyt,diff/analyt
!!$                              mnet_wmat1(m2,cindex) = mnet_wmat1(m2,cindex) - tiny                           
!!$
!!$                           end do
!!$                        end do
!!$                     endif
!!$                  end do
!!$               end do
!!$            end do
!!$         end do
!!$


         do itype = 1, natom_list
            if(use_type(itype)) then 
            !     test derivatives wrt mnet_w
            do layer = 1,mnet_nlayers
               
               if(layer.gt.1) then 
                  write(*,*) 
                  do m1 = 1,mnet_size(layer-1)
                     do m2 = 1, mnet_size(layer)
                        mnet_w(m2,m1,itype,layer) =  mnet_w(m2,m1,itype,layer) + tiny
                        call run_train_mnet
                        
                        error2_energy_1 = error2_energy
                        error2_force_1 = error2_force
                        error2_m_reg_1 = error2_m_reg
                        error2_1 = error2_force + error2_energy + error2_m_reg
                        diff = (error2_1 - error2_0)/tiny
                        
                        write(*, *) 'wdiff', layer, m2, m1, itype, diff, &
                             mnet_de2_dw(m2,m1,itype,layer),diff/mnet_de2_dw(m2,m1,itype,layer),&
                        error2_1,error2_0,error2_energy_0,error2_energy_1,&
                        error2_force_0,error2_force_1,&
                        error2_m_reg_0,error2_m_reg_1
                        
                        call flush (6)
                        mnet_w(m2,m1,itype,layer) =  mnet_w(m2,m1,itype,layer) - tiny
                     end do
                  end do
                  write(*, *)
                  call flush (6)
               endif

!!$               !     test derivatives wrt bias1
!!$
!!$               if(layer.lt.mnet_nlayers) then 
!!$                  do m1 = 1, mnet_size(layer)
!!$                     mnet_b(m1,itype,layer) = mnet_b(m1,itype,layer) + tiny
!!$                     call run_train_mnet
!!$                     
!!$                     error2_1 = error2_force + error2_energy + error2_m_reg
!!$                     diff = (error2_1 - error2_0)/tiny
!!$                     
!!$                     write(*, *) 'bdiff', layer, 'itype',itype,m1, diff, &
!!$                          mnet_de2_db(m1,itype,layer), diff/mnet_de2_db(m1,itype,layer)
!!$                     call flush (6)
!!$                     mnet_b(m1,itype,layer) = mnet_b(m1,itype,layer) - tiny
!!$                  end do
!!$               endif
            end do
         endif
         end do
         !     test derivatives wrt const
               write(*, *)

               do frag_id = 1,nfrag_id

                  mnet_const(frag_id) = mnet_const(frag_id) + tiny
                  call run_train_mnet
               
                  error2_1 = error2_force + error2_energy + error2_m_reg
                  diff = (error2_1 - error2_0)/tiny
                  
                  write(*, *) 'constdiff', frag_id, diff, mnet_de2_dconst(frag_id), diff/mnet_de2_dconst(frag_id)
                  call flush (6)
                  
                  mnet_const(frag_id) = mnet_const(frag_id) - tiny            
            end do
            
      end subroutine test_mnet

      subroutine test_cnet
         implicit none
         integer itype, i, m, n,frag_id
         real(8) error2_0, error2_1, tiny, diff,analyt
         integer id1,id2,id_center,m1,m2,nc_used,cindex,cindex_atom,etype
         
!!$         write(*, *) '*********************'
!!$         write(*, *) '***   CTEST GRADIENTS'
!!$         write(*, *) '*********************'
!!$         call flush (6)
!!$
!!$         tiny = 1d-7
!!$         call run_train_cnet
!!$
!!$         error2_0 = error2_charge + error2_c_reg
!!$
!!$        do id_center = 1,natom_list
!!$           do id1 = 0,natom_list
!!$              do id2 = id1,natom_list
!!$                 do etype = 0,3
!!$                    if(etype.eq.1.and.(.not.use_sumtype)) exit
!!$                    if(col_used(etype,id2,id1,id_center)) then 
!!$                       nc_used = ncol_used(etype,id2,id1,id_center)
!!$                       do m1 = 1,nc_used
!!$                          cindex = col_index(etype,m1,id2,id1,id_center)
!!$                          cindex_atom = col_atom_index(etype,m1,id2,id1,id_center)
!!$                          do m2 = 1,mnet_s1
!!$                             cnet_wmat1(m2,cindex) = cnet_wmat1(m2,cindex) + tiny
!!$                             call run_train_cnet
!!$                             
!!$                             error2_1 = error2_charge + error2_c_reg                             
!!$                             diff = (error2_1 - error2_0)/tiny
!!$                             
!!$                             analyt = cnet_de2_dw1(m2,cindex)
!!$                             write(*, *) 'w1diff',etype,m1,m2,id2,id1,id_center,diff,analyt,diff/analyt
!!$                             cnet_wmat1(m2,cindex) = cnet_wmat1(m2,cindex) - tiny                           
!!$                             
!!$                          end do
!!$                       end do
!!$                    endif
!!$                 end do 
!!$              end do
!!$           end do
!!$        end do
!!$
!!$         do itype = 1, natom_list
!!$            !     test derivatives wrt cnet_wmat2
!!$            
!!$            write(*,*) 
!!$            do m1 = 1, mnet_size(1)
!!$               do m2 = 1,mnet_size(2)
!!$                  cnet_wmat2(m2,m1,itype) = cnet_wmat2(m2,m1,itype) + tiny
!!$                  call run_train_cnet
!!$                  
!!$                  error2_1 = error2_charge + error2_c_reg
!!$                  diff = (error2_1 - error2_0)/tiny
!!$                  
!!$                  write(*, *) 'w2diff', n, diff, cnet_de2_dw2(m2,m1,itype), diff/cnet_de2_dw2(m2,m1,itype)
!!$                  call flush (6)
!!$                  cnet_wmat2(m2,m1,itype) = cnet_wmat2(m2,m1,itype) - tiny
!!$               end do 
!!$            end do
!!$
!!$            write(*, *)
!!$            call flush (6)
!!$            !     test derivatives wrt cnet_b1vec
!!$            
!!$            do m = 1, mnet_s1
!!$               cnet_b1vec(m, itype) = cnet_b1vec(m, itype) + tiny
!!$               call run_train_cnet
!!$               
!!$               error2_1 = error2_charge + error2_c_reg
!!$               diff = (error2_1 - error2_0)/tiny
!!$               
!!$               write(*, *) 'b1diff', m, diff, cnet_de2_db1(m, itype), diff/cnet_de2_db1(m, itype)
!!$               call flush (6)
!!$               cnet_b1vec(m, itype) = cnet_b1vec(m, itype) - tiny
!!$            end do
!!$
!!$
!!$            write(*,*)
!!$            !test derivatives wrt b2vec
!!$
!!$            cnet_b2vec(itype) = cnet_b2vec(itype) + tiny
!!$            call run_train_cnet
!!$            
!!$            error2_1 = error2_charge + error2_c_reg
!!$            diff = (error2_1 - error2_0)/tiny
!!$            write(*, *) 'b2diff', m, diff, cnet_de2_db2(itype), diff/cnet_de2_db2(itype)
!!$            call flush (6)
!!$            
!!$            cnet_b2vec(itype) = cnet_b2vec(itype) - tiny            
!!$            
!!$
!!$         end do !do itype = 1,natom_list
!!$            
!!$
!!$            
          end subroutine test_cnet


      subroutine test_mnet_forces
        implicit none
        real(8) :: energy_pred
        real(8), allocatable, dimension(:) :: charge_pred
        real(8), allocatable, dimension(:,:) :: force_pred, delta_force

        real(8), allocatable, dimension(:,:) :: denergy_dw1
        real(8), allocatable, dimension(:,:,:) :: denergy_dw2,denergy_dw3
        real(8), allocatable, dimension(:,:) :: denergy_db1,denergy_db2,denergy_db3
        real(8), allocatable, dimension(:) :: denergy_dconst
        
        real(8), allocatable, dimension(:,:,:) :: denergy_db
        real(8), allocatable, dimension(:,:,:,:) :: denergy_dw
         
        real(8), allocatable, dimension(:,:,:,:) :: dcharge_dw1
        real(8), allocatable, dimension(:,:,:,:) :: dcharge_dw2,dcharge_dw3
        real(8), allocatable, dimension(:,:,:) :: dcharge_db1,dcharge_db2,dcharge_db3

        real(8), allocatable, dimension(:,:,:,:,:) :: dcharge_dw
        real(8), allocatable, dimension(:,:,:,:) :: dcharge_db

        real(8), allocatable, dimension(:,:,:,:,:) :: dforce_dw2,dforce_dw3
        real(8), allocatable, dimension(:,:,:,:) :: dforce_db1,dforce_db2,dforce_db3
        
        real(8), allocatable, dimension(:,:,:,:,:,:) :: dforce_dw
        real(8), allocatable, dimension(:,:,:,:,:) :: dforce_db

        real(8), allocatable, dimension(:,:,:,:) :: dforce_dw1
        
        real(8) :: u0,u1,numerical,analytical
        real(8) tiny,net_timing
        integer ns,ivec,iatom

        tiny = 1.d-6

        allocate(denergy_dw1(mnet_size(1),ncols_wmat1))

        allocate(denergy_dw2(mnet_size(2),mnet_size(1), natom_list))
        allocate(denergy_dw3(mnet_size(3),mnet_size(2), natom_list))

        allocate(denergy_db1(mnet_size(1), natom_list))
        allocate(denergy_db2(mnet_size(2), natom_list))
        allocate(denergy_db3(mnet_size(3), natom_list))
        

        allocate(denergy_db(mnet_size_max,natom_list,mnet_nlayers))
        allocate(denergy_dw(mnet_size_max,mnet_size_max,natom_list,mnet_nlayers))
        allocate(dforce_dw(3,mnet_size_max,mnet_size_max, natoms_max, natom_list,mnet_nlayers))
        allocate(dforce_db(3,mnet_size_max, natoms_max, natom_list,mnet_nlayers))
        allocate(denergy_dconst(128))

        allocate(force_pred(3, natoms_max))
        allocate(dforce_dw2(3,mnet_size(2),mnet_size(1),natoms_max,natom_list))

        allocate(dforce_db1(3,mnet_size(1),natoms_max,natom_list))

        allocate(dforce_dw1(3,ncols_wmat1,mnet_size(1),natoms_max))
        
        write(*,*) 'TEST FORCES'

        do ns = 1,nsteps_total        

            call get_mnet_energy(ns, net_type, &
                 energy_pred, charge_pred, &
                 denergy_dw1, denergy_dw, denergy_db, & 
                 dcharge_dw1, dcharge_dw, dcharge_db, &
                 dforce_dw1, dforce_dw, dforce_db, & 
                 denergy_dconst,force_pred,net_timing)

           u0 = energy_pred

           do iatom = 1,natoms_step(ns)
              do ivec = 1,3
                 
                 coord(ivec,iatom,ns) = coord(ivec,iatom,ns) + tiny
                 call get_mnet_energy(ns, net_type, &
                      energy_pred, charge_pred, &
                      denergy_dw1, denergy_dw, denergy_db, & 
                      dcharge_dw1, dcharge_dw, dcharge_db, &
                      dforce_dw1, dforce_dw, dforce_db, & 
                      denergy_dconst,force_pred,net_timing)


                 u1 = energy_pred
                 numerical = (u0-u1)/tiny
                 analytical = force_pred(ivec,iatom)
                 write(*,*) ivec,iatom,'numerical',numerical,'analytical',analytical,numerical/analytical
                 coord(ivec,iatom,ns) = coord(ivec,iatom,ns) - tiny        
              end do
           end do
        end do
        
        stop
      end subroutine test_mnet_forces

      
      subroutine get_sigma(xvec, sigvec, dsigvec, d2sigvec)
         implicit none
         real(8) :: efac, eplus, eminus
         real(8), dimension(:) :: xvec, sigvec, dsigvec, d2sigvec
!tanh         
!         sigvec = dtanh(xvec)
!         dsigvec = 1.0d0 - sigvec**2
!         d2sigvec = -2.0d0*sigvec*dsigvec
         !logistic
         sigvec = 1.0d0 / (1.0d0 + dexp(-xvec))
         dsigvec = sigvec * (1.0d0 - sigvec)
         d2sigvec = dsigvec * (1.0d0 -2.0d0 * sigvec)

      end subroutine get_sigma

end module networks_mod


subroutine allocate(net)
   use molneural_data
   use parse_text
   implicit none
   character(16) :: net
   integer m

!      write(*,*) 'ALLOCATE'
   
   !save the current value of neighbor, so that we can compare if it changes later.
   
   neighbor0 = neighbor
   mnet_allocated = .true.
   
   cnet_nlayers = 2
   
   allocate(cnet_size(0:cnet_nlayers))


   mnet_size_max = 0 
   do m = 1,mnet_nlayers
      if(mnet_size(m).gt.mnet_size_max) mnet_size_max = mnet_size(m)
   end do
   
   
   cnet_size(0) = 1
   cnet_size(1) = mnet_size(1)
   cnet_size(2) = 1
   
   !allocate the neural network parameter arrays
   
   allocate(mnet_w(mnet_size_max,mnet_size_max,natom_list,mnet_nlayers))
   allocate(mnet_b(mnet_size_max,natom_list,mnet_nlayers))
   allocate (mnet_const(1024))
   allocate(mnet_de2_db(mnet_size_max,natom_list,mnet_nlayers))
   allocate(mnet_de2_dw(mnet_size_max,mnet_size_max,natom_list,mnet_nlayers))
   
   mnet_w = 0.0d0 
   mnet_b = 0.0d0 
   mnet_const = 0.0d0 
   
   allocate(cnet_wmat2(mnet_size(2),mnet_size(1),natom_list))
   allocate(cnet_b1vec(mnet_s1max, natom_list))
   allocate(cnet_b2vec(1024))
   
   
   allocate(mnet_de2_dconst(1024))

end subroutine allocate

      real(8) function d2sigma(sigma)
         real(8) :: sigma

         d2sigma = sigma*(1.0d0 - sigma)*(1.0d0 - 2.0d0*sigma)

      end function d2sigma

      subroutine randomize(y, add)
         use molneural_data
         implicit none
         integer i, j, itype, npairtype,nc,frag_id,cindex
         real(8) :: rand
         real(8) :: y
         integer, dimension(43) :: xx
         integer m1,m2,m3,layer
         logical add


         xx = 0.0d0

         do nc = 1,ncols_wmat1
            do m1 = 1,mnet_size(1)
               call random_number(rand)
               if(add) then 
                  mnet_wmat1(m1,nc) = mnet_wmat1(m1,nc) + (rand - 0.5d0)*y
                  cnet_wmat1(m1,nc) = cnet_wmat1(m1,nc) + (rand - 0.5d0)*y
               else
                  mnet_wmat1(m1,nc) = (rand - 0.5d0)*y
                  cnet_wmat1(m1,nc) = (rand - 0.5d0)*y
               endif
            end do
         end do

         do itype = 1, natom_list
            do layer = 1,mnet_nlayers

               if(layer.lt.mnet_nlayers) then 
                  do m1 = 1, mnet_size(layer)
                     call random_number(rand)
                     if (add) then
                        mnet_b(m1,itype,layer) = mnet_b(m1,itype,layer) + (rand - 0.5)*y
                     else
                        mnet_b(m1,itype,layer) = (rand - 0.5d0)*y
                     end if
                  end do
               endif

               if(layer.gt.1) then 
                  do m2 = 1, mnet_size(layer)
                     do m1 = 1,mnet_size(layer-1)
                        call random_number(rand)
                        if (add) then
                           mnet_w(m2,m1,itype,layer) = mnet_w(m2,m1,itype,layer) + (rand - 0.5d0)*y
                           !                  cnet_wmat2(j, itype) = cnet_wmat2(j, itype) + (rand - 0.5d0)*y
                        else
                           mnet_w(m2,m1,itype,layer) = (rand - 0.5d0)*y
                           !                  cnet_wmat2(j, itype) = (rand - 0.5d0)*y
                        end if
                     end do
                  end do
               endif

            end do 

         end do ! itype

         do frag_id = 1,nfrag_id
            call random_number(rand)
            if (add) then
               mnet_const(frag_id) = mnet_const(frag_id) + (rand - 0.5)*y
!               cnet_b2vec(frag_id) = cnet_b2vec(frag_id) + (rand - 0.5)*y
            else
               mnet_const(frag_id) = (rand - 0.5)*y
!               cnet_b2vec(frag_id) = (rand - 0.5)*y               
            end if
         end do



       end subroutine randomize

      subroutine read_train_data
         use molneural_data
         use parse_text
         implicit none
         integer ns, nstep, ios, iatom, iat, itype, jtype, m, n, k, j, nn, nn1, nn2
         integer atno1, atno2,jatom
         character(len=16) :: string, dstring
         real(8) :: fx, fy, fz
         integer i, ntrain, ntest, npairtype
         real(8) :: rand, xdif, ydif, zdif, rdis, energy
         character(len=32) :: name
         character(len=32), dimension(0:10) :: tlist
         integer fileno
         integer get_pairtype, nitems
         character(len=2048) :: buffer
         character(len=2048), dimension(64) :: textlist
         character(64) :: filename
         integer ibin,atomid
         integer, dimension(8192) :: count_b, count_nb
         integer type1, type2
         integer mm, kk
         integer nat,struct_id
         logical pass
         character(3) :: atname
         

         !      write(*,*) 'READ TRAIN DATA'

         ntrain = 0
         rewind(10)
         rewind(11)

         use_type = .false. 
         do ns = 1,nsteps_total

            read(10,*,iostat = ios) nat
            
            if(ios.ne.0) exit

            natoms_step(ns) = nat

            read(11,*) energy
            train_energy(ns) = energy*evfac
            
            do iatom = 1,nat
               read(10,*) atname,coord(:,iatom,ns),fx,fy,fz,charge(iatom,ns)
               call get_atom_id(atname,atomid)
               atomid_step(iatom,ns) = atomid
               use_type(atomid) = .true.
               atname_step(iatom,ns) = atname
               train_force(1, iatom, ns) = fx*evfac_f
               train_force(2, iatom, ns) = fy*evfac_f
               train_force(3, iatom, ns) = fz*evfac_f               
            end do

         end do !do while(.true.)

         call flush (6)

         natoms_train = 0
         natoms_val = 0 


         do ns = 1, nsteps_total
            struct_id = struct_id_ns(ns)
            if(train_step(ns)) then 
               natoms_train = natoms_train + natoms_step(ns)
            else
               natoms_val = natoms_val + natoms_step(ns)
            endif
         end do
         if(.not.restart) then 
            write(*,*) 'nsteps_train',nsteps_train
            write(*,*) 'nsteps_val',nsteps_val
            call flush(6)
         endif
       end subroutine read_train_data


      subroutine minimize
         use nr_common_data
         use nr_mod
         use molneural_data
         implicit none
         integer i, j, n, itype, mtype, frag_id
         real(8) :: xvec(60000)
         real(8) :: ftol, fret, uu, uu0, uu1, vv0, vv1, ww0, ww1, func
         real start, finish
         integer iter, npairtype
         integer id_center,id1,id2,nc_used,m1,m2,layer,cindex,etype

         n = 0 

         if(net_type.eq.'mnet') then 

            if(minmode.eq.'all') then 
               do id_center = 1,natom_list
                  !charge term
                  cindex = charge_col_index(id_center)
                  do m2 = 1,mnet_size(1)
                     n = n + 1
                     xvec(n) = mnet_wmat1(m2,cindex)
                  end do
                  
                  do id1 = 0,natom_list
                     do id2 = id1,natom_list
                        do etype = 0,3
                           if(etype.eq.1.and.(.not.use_sumtype)) exit
                           if(col_used(etype,id2,id1,id_center)) then 
                              nc_used = ncol_used(etype,id2,id1,id_center)
                              do m1 = 1,nc_used
                                 cindex = col_index(etype,m1,id2,id1,id_center)
                                 do m2 = 1,mnet_size(1)
                                    n = n + 1
                                    xvec(n) = mnet_wmat1(m2,cindex)
                                 end do
                              end do
                           endif
                        end do
                     end do
                  end do
               end do
               
               do layer = 1,mnet_nlayers
                  do itype = 1, natom_list
                     if(use_type(itype)) then 

                     if(layer.gt.1) then 
                        do m1 = 1, mnet_size(layer-1)
                           do m2 = 1,mnet_size(layer)
                              n = n + 1
                              xvec(n) = mnet_w(m2,m1,itype,layer)
                           end do
                        end do
                     endif
                     if(layer.lt.mnet_nlayers) then 
                        do m1 = 1,mnet_size(layer)
                           n = n + 1
                           xvec(n) = mnet_b(m1,itype,layer) 
                        end do
                     endif
                  endif
                  end do
               end do
            endif

      do frag_id = 1,nfrag_id
         n = n + 1
         xvec(n) = mnet_const(frag_id)
      end do
      else if(net_type.eq.'cnet') then
!!$         if(minmode.eq.'all') then 
!!$            do id_center = 1,natom_list
!!$               do id1 = 0,natom_list
!!$                  do id2 = id1,natom_list
!!$                     do etype = 0,3
!!$                        if(etype.eq.1.and.(.not.use_sumtype)) exit
!!$                        if(col_used(etype,id2,id1,id_center)) then 
!!$                           nc_used = ncol_used(etype,id2,id1,id_center)
!!$                           do m1 = 1,nc_used
!!$                              cindex = col_index(etype,m1,id2,id1,id_center)
!!$                              do m2 = 1,mnet_size(1)
!!$                                 n = n + 1
!!$                                 xvec(n) = cnet_wmat1(m2,cindex) 
!!$                              end do
!!$                           end do
!!$                        endif
!!$                     end do
!!$                  end do
!!$               end do
!!$            end do
!!$         endif
!!$         do itype = 1, natom_list
!!$            if(minmode.eq.'all') then 
!!$               do i = 1, mnet_size(1)
!!$                  n = n + 1
!!$                  xvec(n) = cnet_wmat2(i, itype) 
!!$               end do
!!$               do i = 1, mnet_size(1)
!!$                  n = n + 1
!!$                  xvec(n) = cnet_b1vec(i, itype) 
!!$               end do
!!$            endif
!!$            n = n + 1
!!$            xvec(n) = cnet_b2vec(itype) 
!!$         end do
!!$         
      endif
         
         !      call cpu_time(start)
         !      do i = 1,1024
         !         uu = func(xvec)
         !         write(*,*) i,'uu',error_train,error_test
         !      end do
         !      call cpu_time(finish)
         !      write(*,*) 'time = ',finish - start
         !      stop
         
      ftol = 0.0d0
      linmin_param = 1.0d-6 * 1.d+3
      nprintdata = 1
      
      frprmn_itmax = itmax
      if(minmode.eq.'energy') frprmn_itmax = 20
      
      if(net_type.eq.'mnet') then 
         write(*, *) 'MINIMIZING MNET WITH ', frprmn_itmax, 'STEPS'
      else if(net_type.eq.'cnet') then
         write(*,*) 'MINIMIZING CNET WITH ', frprmn_itmax, 'STEPS'
      endif
            
      call flush (6)
      
      if(.not.restart) nminsteps = 0
      call frprmn(xvec, n, ftol, iter, fret)
      
      call print_energies_and_forces
      
    end subroutine minimize
    
    real(8) function func(xvec)
      use nr_common_data
      use nr_mod
      use molneural_data
      use networks_mod
      real(8) :: xvec(60000)
      integer i, j, n, itype, npairtype,frag_id
      integer id_center,id1,id2,nc_used,m1,m2,cindex,etype
      
      
      n = 0
      
      if(net_type.eq.'mnet') then 
      if(minmode.eq.'all') then 
         do id_center = 1,natom_list
            !charge term
            cindex = charge_col_index(id_center)
            do m2 = 1,mnet_size(1)
               n = n + 1
               mnet_wmat1(m2,cindex) = xvec(n)
            end do
            
            do id1 = 0,natom_list
               do id2 = id1,natom_list
                  do etype = 0,3
                     if(etype.eq.1.and.(.not.use_sumtype)) exit
                     if(col_used(etype,id2,id1,id_center)) then 
                        nc_used = ncol_used(etype,id2,id1,id_center)
                        do m1 = 1,nc_used
                           cindex = col_index(etype,m1,id2,id1,id_center)
                           do m2 = 1,mnet_size(1)
                              n = n + 1
                              mnet_wmat1(m2,cindex) = xvec(n)
                           end do
                        end do
                     endif
                  end do
               end do
            end do
         end do
         
         do layer = 1,mnet_nlayers
            do itype = 1, natom_list
               if(use_type(itype)) then 
               if(layer.gt.1) then 
                  do m1 = 1, mnet_size(layer-1)
                     do m2 = 1,mnet_size(layer)
                        n = n + 1
                        mnet_w(m2,m1,itype,layer) = xvec(n)
                     end do
                  end do
               endif
               if(layer.lt.mnet_nlayers) then 
                  do m1 = 1,mnet_size(layer)
                     n = n + 1
                     mnet_b(m1,itype,layer) = xvec(n)
                  end do
               endif
            endif
         end do
         end do
      endif

      do frag_id = 1,nfrag_id
         n = n + 1
         mnet_const(frag_id) = xvec(n)
      end do
      else if(net_type.eq.'cnet') then 
!!$            if(minmode.eq.'all') then 
!!$               do id_center = 1,natom_list
!!$                  do id1 = 0,natom_list
!!$                     do id2 = id1,natom_list
!!$                        do etype = 0,3
!!$                           if(etype.eq.1.and.(.not.use_sumtype)) exit
!!$                           if(col_used(etype,id2,id1,id_center)) then 
!!$                              nc_used = ncol_used(etype,id2,id1,id_center)
!!$                              do m1 = 1,nc_used
!!$                                 cindex = col_index(etype,m1,id2,id1,id_center)
!!$                                 do m2 = 1,mnet_size(1)
!!$                                    n = n + 1
!!$                                    cnet_wmat1(m2,cindex) = xvec(n)
!!$                                 end do
!!$                              end do
!!$                           endif
!!$                        end do
!!$                     end do
!!$                  end do
!!$               end do
!!$            endif
!!$            do itype = 1, natom_list
!!$               if(minmode.eq.'all') then 
!!$                  do i = 1, mnet_size(1)
!!$                     n = n + 1
!!$                     cnet_wmat2(i, itype) = xvec(n)
!!$                  end do
!!$                  do i = 1, mnet_size(1)
!!$                     n = n + 1
!!$                     cnet_b1vec(i, itype) = xvec(n)
!!$                  end do
!!$               endif
!!$               n = n + 1
!!$               cnet_b2vec(itype) = xvec(n)
!!$            end do
         endif
      
         if(net_type.eq.'mnet') then 
            call run_train_mnet
            func = error2_force + error2_energy + error2_m_reg
         else if(net_type.eq.'cnet') then 
            call run_train_cnet
            func = error2_charge + error2_c_reg
         endif
      
    end function func
       
    subroutine dfunc(xvec, gvec)
      use nr_common_data
      use nr_mod
      use molneural_data
      implicit none
      real(8) :: xvec(60000), gvec(60000)
      integer i, j, n, itype, npairtype,frag_id
      integer id_center,id1,id2,nc_used,m1,m2,layer,cindex,etype
      
      n = 0
      
      if(net_type.eq.'mnet') then 

      if(minmode.eq.'all') then 
         do id_center = 1,natom_list
            
            !charge term
            cindex = charge_col_index(id_center)
            do m2 = 1,mnet_size(1)
               n = n + 1
               gvec(n) = mnet_de2_dw1(m2,cindex) 
            end do
            
            do id1 = 0,natom_list
               do id2 = id1,natom_list
                  do etype = 0,3
                     if(etype.eq.1.and.(.not.use_sumtype)) exit
                     if(col_used(etype,id2,id1,id_center)) then 
                        nc_used = ncol_used(etype,id2,id1,id_center)
                        do m1 = 1,nc_used
                           cindex = col_index(etype,m1,id2,id1,id_center)
                           do m2 = 1,mnet_size(1)
                              n = n + 1
                              gvec(n) = mnet_de2_dw1(m2,cindex) 
                           end do
                        end do
                     endif
                  end do
               end do
            end do
         end do
         
         do layer = 1,mnet_nlayers
            do itype = 1, natom_list
               if(use_type(itype)) then 
               if(layer.gt.1) then 
                  do m1 = 1, mnet_size(layer-1)
                     do m2 = 1,mnet_size(layer)
                        n = n + 1
                        gvec(n) = mnet_de2_dw(m2,m1,itype,layer) 
                     end do
                  end do
               endif
               if(layer.lt.mnet_nlayers) then 
                  do m1 = 1,mnet_size(layer)
                     n = n + 1
                     gvec(n) = mnet_de2_db(m1,itype,layer) 
                  end do
               endif
               endif
            end do
         end do
      endif

      do frag_id = 1,nfrag_id
         n = n + 1
         gvec(n) = mnet_de2_dconst(frag_id)
      end do
      else if(net_type.eq.'cnet') then 
!!$            if(minmode.eq.'all') then 
!!$               do id_center = 1,natom_list
!!$                  do id1 = 0,natom_list
!!$                     do id2 = id1,natom_list
!!$                        do etype = 0,3
!!$                           if(etype.eq.1.and.(.not.use_sumtype)) exit
!!$                           if(col_used(etype,id2,id1,id_center)) then 
!!$                              nc_used = ncol_used(etype,id2,id1,id_center)
!!$                              do m1 = 1,nc_used
!!$                                 cindex = col_index(etype,m1,id2,id1,id_center)
!!$                                 do m2 = 1,mnet_size(1)
!!$                                    n = n + 1
!!$                                    gvec(n) = cnet_de2_dw1(m2,cindex) 
!!$                                 end do
!!$                              end do
!!$                           endif
!!$                        end do
!!$                     end do
!!$                  end do
!!$               end do
!!$            endif
!!$            
!!$            do itype = 1, natom_list
!!$               if(minmode.eq.'all') then 
!!$                  do i = 1, mnet_size(1)
!!$                     n = n + 1
!!$                     gvec(n) = cnet_de2_dw2(i, itype) 
!!$                  end do
!!$                  do i = 1, mnet_size(1)
!!$                     n = n + 1
!!$                     gvec(n) = cnet_de2_db1(i, itype) 
!!$                  end do
!!$               endif
!!$               n = n + 1
!!$               gvec(n) = cnet_de2_db2(itype) 
!!$            end do
      endif

    end subroutine dfunc
    
    subroutine print_results
      use molneural_data
      use networks_mod
      implicit none
      integer n
      
      print_forces = .false.
      if (mod(nminstep, 1) .eq. 0) then
         write(*, *) 'PRINT FORCES'
         call flush (6)
         print_forces = .true.
         call open_force_energy_files
      end if
      
      if(net_type.eq.'mnet') then 
         call run_train_mnet
      else if(net_type.eq.'cnet') then
         call run_train_cnet
      endif
      
      print_forces = .false.
      
      if(mod(nminstep,1).eq.0) then 
         call print_error(0)
         call flush(6)
      endif
      
      if (mod(nminstep, 1) .eq. 0) then
         
         if(net_type.eq.'mnet') then 
         call print_mnet(20)
         else if(net_type.eq.'cnet') then 
            call print_cnet(30)
         endif

      end if
      
    end subroutine print_results
    
    subroutine print_data
      use molneural_data
      implicit none
      !     dummy subroutine for nr_mod call-back during minimization. But can be used to
      !     print results during a relaxation.
      nminsteps = nminsteps + 1
      
      if(mod(nminsteps,10).eq.0) then
         call flush(6) 

            if(net_type.eq.'mnet') then 
               call print_energies_and_forces
            else if(net_type.eq.'cnet') then 
               call print_charges
            endif
         
            if(net_type.eq.'mnet') then 
               call print_mnet(20)
            else if(net_type.eq.'cnet') then 
               call print_cnet(30)
            endif               

         if (error_val.lt.min_error_val) then
            min_error_val = error_val
               if(net_type.eq.'mnet') then 
                  call print_mnet(24)
               else if(net_type.eq.'cnet') then 
                  call print_cnet(34)
               endif
         endif
      endif
      
      call print_error(nminsteps)
    end subroutine print_data
    

    subroutine read_mnet0
      use molneural_data
      use parse_text
      implicit none
      character(len=32) :: string,pstring
      character(len=2048) :: buffer
      character(len=2048), dimension(64) :: textlist
      character(len=32), dimension(0:10) :: tlist
      integer itype, i, j, k, n1, n2, ios, nitems, ii, molmax0,n
      integer nc,ncols,etype,nfid
      integer atomid_central,id1,id2,cindex,m1,m2,nstruct,frag_id,id
      real(8) :: uu
      character(len=32) :: name, net
      character(len = 3) :: atname,atname1,atname2,atname3
      real(8), allocatable, dimension(:) :: vec
      real(8) :: rand,p
      logical begin_pass
      real(8) :: scale       
      
      begin_pass = .false.
      do while(.true.)
         read(22,*,iostat = ios) string
         if(ios.ne.0) exit
         if(string.eq.'BEGIN') then
            begin_pass = .true.
            exit
         endif
      end do
      if(.not.begin_pass) then
         write(*,*) "ERROR. 'BEGIN DATA' NOT FOUND IN MNET.DAT"
         stop
      endif
         
         net_type = "mnet"

         if(.not.restart) write(*,*) 'READ_MNET'
         do while(.true.)
            read(22,*,iostat = ios) string,atname
            if(ios.ne.0) exit

            if(string.eq.'MNET_BIAS2') then 
               backspace(22)
               read(22,*) string,nfid
               do frag_id = 1,nfid
                  read(22,*) i,mnet_const(frag_id)
               end do
               exit
            endif


            call get_atom_id(atname,atomid_central)
            read(22,*) string
            
!   charge term
            read(22,*) string
            if(atomid_central.ne.0) then 
               cindex = charge_col_index(atomid_central)
            else
               cindex = 0
            endif

            if(cindex.ne.0) then 
               read(22,*) (mnet_wmat1(m2,cindex),m2=1,mnet_size(1))
            else
               read(22,*)
            endif

            do while(.true.)
               
               read(22,*) string
               if(string.eq.'END') exit
               backspace 22
               read(22,*) pstring,atname1,atname2,string,atname3,string,etype
               if(etype.eq.1.and.(.not.use_sumtype)) then
                  write(*,*) 'ERROR: MNET.DAT USES SUMTYPES, BUT SUMTYPE NOT SET IN CONTROL FILE'
                  stop
               endif
               
               call get_atom_id(atname1,id1)
               call get_atom_id(atname2,id2)         
               if(pstring.eq.'cpair:') id1 = 0 
               read(22,*) string,ncols
               do m1 = 1,ncols
                  if(atomid_central.ne.0) then 
                     cindex = col_index(etype,m1,id2,id1,atomid_central)
                  else
                     cindex = 0
                  endif
                  if(cindex.ne.0) then 
                     read(22,*) (mnet_wmat1(m2,cindex),m2=1,mnet_size(1))
                     do m2 = 1,mnet_size(1)
                        if(mnet_wmat1(m2,cindex).gt.1.d+3) mnet_wmat1(m2,cindex) = 1.0d0
                        if(mnet_wmat1(m2,cindex).lt.-1.d+3) mnet_wmat1(m2,cindex) = -1.0d0
                     end do
                  else
                     read(22,*) 
                  endif
               end do
            end do !do while(.true.)

            read(22,*) string
            if(atomid_central.ne.0) then 
               read(22,*) (mnet_w(1,m1,atomid_central,2), m1=1,mnet_size(1)) 
            else
               read(22,*)
            endif

            read(22,*) string

            if(atomid_central.ne.0) then 
               read(22,*) (mnet_b(m1,atomid_central,1),m1=1,mnet_size(1))
            else
               read(22,*)
            endif
!            read(22,*) string
!            read(22,*) mnet_const(atomid_central)
         end do
            
         !for putting on a 2nd layer of neurons


         scale = 1000.0d0 
         do id = 1,natom_list
            do m1 = 1,mnet_size(1)
               do m2 = 1,mnet_size(2)
                  call random_number(rand)
                  mnet_w(m2,m1,id,2) = mnet_w(1,m1,id,2) + (rand - 0.5d0)* 0.1d0 
               end do
            end do
         end do
         mnet_w(:,:,:,2) = mnet_w(:,:,:,2) / scale

         do itype = 1,natom_list
            do m2 = 1,mnet_size(2)
               call random_number(rand)
               mnet_w(1,m2,itype,3) = 4.0d0*scale/mnet_size(2) + (rand - 0.5d0) * 0.1d0 
            end do 
         end do



!!$         scale = 100.0d0 
!!$         mnet_w(:,:,:,2) = mnet_w(:,:,:,2) / scale
!!$
!!$         do itype = 1,natom_list
!!$            call random_number(rand)
!!$            mnet_w(1,1,itype,3) = 4.0d0 * scale
!!$         end do


      end subroutine read_mnet0


    subroutine read_mnet
      use molneural_data
      use parse_text
      implicit none
      character(len=32) :: string,pstring,string1,string2
      character(len=2048) :: buffer
      character(len=2048), dimension(64) :: textlist
      character(len=32), dimension(0:10) :: tlist
      integer itype, i, j, k, n1, n2, ios, nitems, ii, molmax0,n
      integer nc,ncols,etype,nfid,layer
      integer atomid_central,id1,id2,cindex,m1,m2,nstruct,frag_id
      real(8) :: uu
      character(len=32) :: name, net
      character(len = 3) :: atname,atname1,atname2,atname3
      real(8), allocatable, dimension(:) :: vec
      real(8) :: rand,p
      logical begin_pass
      
      begin_pass = .false.
      do while(.true.)
         read(22,*,iostat = ios) string
         if(ios.ne.0) exit
         if(string.eq.'BEGIN') then
            begin_pass = .true.
            exit
         endif
      end do
      if(.not.begin_pass) then
         write(*,*) "ERROR. 'BEGIN DATA' NOT FOUND IN MNET.DAT"
         stop
      endif
         
         net_type = "mnet"

         if(.not.restart) write(*,*) 'READ_MNET'
         do while(.true.)
            read(22,*,iostat = ios) string,atname
            if(ios.ne.0) exit

            if(string.eq.'MNET_BIAS2') then 
               backspace(22)
               read(22,*) string,nfid
               do frag_id = 1,nfid
                  read(22,*) i,mnet_const(frag_id)
               end do
               exit
            endif


            call get_atom_id(atname,atomid_central)
            read(22,*) string
            
!   charge term
            read(22,*) string
            if(atomid_central.ne.0) then 
               cindex = charge_col_index(atomid_central)
            else
               cindex = 0
            endif

            if(cindex.ne.0) then 
               read(22,*) (mnet_wmat1(m2,cindex),m2=1,mnet_size(1))
            else
               read(22,*)
            endif

            do while(.true.)
               read(22,*) string
               if(string.eq.'END') exit
               backspace 22
               read(22,*) pstring,atname1,atname2,string,atname3,string,etype
               if(etype.eq.1.and.(.not.use_sumtype)) then
                  write(*,*) 'ERROR: MNET.DAT USES SUMTYPES, BUT SUMTYPE NOT SET IN CONTROL FILE'
                  stop
               endif
               
               call get_atom_id(atname1,id1)
               call get_atom_id(atname2,id2)         
               if(pstring.eq.'cpair:') id1 = 0 
               read(22,*) string,ncols
               do m1 = 1,ncols
                  if(atomid_central.ne.0) then 
                     cindex = col_index(etype,m1,id2,id1,atomid_central)
                  else
                     cindex = 0
                  endif
                  if(cindex.ne.0) then 
                     read(22,*) (mnet_wmat1(m2,cindex),m2=1,mnet_size(1))
                     do m2 = 1,mnet_size(1)
                        if(mnet_wmat1(m2,cindex).gt.1.d+3) mnet_wmat1(m2,cindex) = 1.0d0
                        if(mnet_wmat1(m2,cindex).lt.-1.d+3) mnet_wmat1(m2,cindex) = -1.0d0
                     end do
                  else
                     read(22,*) 
                  endif
               end do
            end do !do while(.true.)

            do while(.true.) 
               read(22,*) string

               if(string.eq.'WEIGHT') then 
                  backspace(22) 
                  read(22,*) string,string1,string2,layer
                  do m2 = 1,mnet_size(layer)
                     if(atomid_central.ne.0) then 
                        read(22,*) (mnet_w(m2,m1,atomid_central,layer), m1=1,mnet_size(layer-1)) 
                     else
                        read(22,*)
                     endif
                  end do
               else if(string.eq.'BIAS') then 
                  backspace(22) 
                  read(22,*) string,string1,string2,layer                  
                  if(atomid_central.ne.0) then 
                     read(22,*) (mnet_b(m1,atomid_central,layer),m1=1,mnet_size(layer))
                  else
                     read(22,*)
                  endif
               else if(string.eq.'ATOM:'.or.string.eq.'END') then 
                  backspace(22)
                  exit
               endif
            end do

            if(string.eq.'END') then  
               read(22,*) string
               read(22,*) string,nfid
               mnet_const = 0.0d0 
               do n = 1,nfid
                  read(22,*) i,mnet_const(n)
               end do
               exit
            endif
         end do
      end subroutine read_mnet



      subroutine read_cnet
         use molneural_data
         use parse_text
         implicit none
         character(len=32) :: string,pstring
         character(len=2048) :: buffer
         character(len=2048), dimension(64) :: textlist
         character(len=32), dimension(0:10) :: tlist
         integer itype, i, j, k, n1, n2, ios, nitems, ii, molmax0
         integer nc,ncols,etype
         integer atomid_central,id1,id2,cindex,m1,m2,nstruct,frag_id,nfrag
         real(8) :: uu
         logical begin_pass
         character(len=32) :: name, net
         character(len = 3) :: atname,atname1,atname2,atname3
         
!!$         begin_pass = .false.
!!$         do while(.true.)
!!$            read(32,*,iostat = ios) string
!!$            if(ios.ne.0) exit
!!$            if(string.eq.'BEGIN') then
!!$               begin_pass = .true.
!!$               exit
!!$            endif
!!$         end do
!!$         if(.not.begin_pass) then
!!$            write(*,*) "ERROR. 'BEGIN DATA' NOT FOUND IN CNET.DAT"
!!$            stop
!!$         endif
!!$
!!$         net_type = "cnet"
!!$
!!$         write(*,*) 'READ_CNET'
!!$         do while(.true.)
!!$            read(32,*,iostat = ios) string,atname
!!$            if(ios.ne.0) exit
!!$
!!$            if(string.eq.'CNET_BIAS2') then 
!!$               backspace(32)
!!$               read(32,*) string,nfrag
!!$               do frag_id = 1,nfrag
!!$                  read(32,*) i,cnet_b2vec(frag_id)
!!$               end do
!!$               exit
!!$            endif
!!$
!!$
!!$            call get_atom_id(atname,atomid_central)
!!$            read(32,*) string
!!$            
!!$            do while(.true.)
!!$               
!!$               read(32,*) string
!!$               if(string.eq.'END') exit
!!$               backspace 32
!!$               read(32,*) pstring,atname1,atname2,string,atname3,string,etype
!!$               if(etype.eq.1.and.(.not.use_sumtype)) then
!!$                  write(*,*) 'ERROR: CNET.DAT USES SUMTYPES, BUT SUMTYPE NOT SET IN CONTROL FILE'
!!$                  stop
!!$               endif
!!$               
!!$               call get_atom_id(atname1,id1)
!!$               call get_atom_id(atname2,id2)         
!!$               if(pstring.eq.'cpair:') id1 = 0 
!!$               read(32,*) string,ncols
!!$               do m1 = 1,ncols
!!$                  cindex = col_index(etype,m1,id2,id1,atomid_central)
!!$                  if(cindex.ne.0) then 
!!$                     read(32,*) (cnet_wmat1(m2,cindex),m2=1,mnet_s1)
!!$                     do m2 = 1,mnet_s1
!!$                        if(cnet_wmat1(m2,cindex).gt.1.d+3) cnet_wmat1(m2,cindex) = 1.0d0
!!$                        if(cnet_wmat1(m2,cindex).lt.-1.d+3) cnet_wmat1(m2,cindex) = -1.0d0
!!$                     end do
!!$                  else
!!$                     read(32,*) 
!!$                  endif
!!$               end do
!!$            end do !do while(.true.)
!!$            read(32,*) string
!!$            read(32,*) (cnet_wmat2(m2,atomid_central),m2=1,mnet_s1)
!!$            read(32,*) string
!!$            read(32,*) (cnet_b1vec(m2,atomid_central),m2=1,mnet_s1)         
!!$            read(32,*) string
!!$            read(32,*) cnet_b2vec(atomid_central)
!!$         end do
!!$            
      end subroutine read_cnet

      
      subroutine print_mnet(fileno)
         use molneural_data
         implicit none
         integer itype, i, j, n, nc_used,id1,id2,atomid_central,cindex,m1,m2,id,etype
         integer fileno,frag_id,layer
         
         write(*,*) 'PRINT_MNET'
         
         REWIND(fileno)
         write(fileno, *) '**** MNETWORK ****', ' error: ftrain ', error_train
         write(fileno, *) 'NLAYERS ',mnet_nlayers
         write(fileno, *) 'SIZE',(mnet_size(m1),m1=1,mnet_nlayers)
         write(fileno, *) 'RADIAL_CUT  ', radial_cut
         write(fileno, *) 'ARC_CUT  ', arc_cut
         write(fileno, *) 'NCOLS_TOTAL ',ncols_wmat1
         write(fileno, *) 'MNET_ALPHA ',mnet_alpha
         write(fileno, *) 
         write(fileno, *) 'BEGIN DATA'

         do atomid_central = 1,natom_list
            write(fileno,*) 'ATOM: ',atom_list(atomid_central)
            write(fileno, *) 'WEIGHT | LAYER',1

            write(fileno,*) 'charge'
            cindex = charge_col_index(atomid_central)
            write(fileno,*) (mnet_wmat1(m2,cindex),m2 = 1,mnet_size(1))
            
            do id1 = 0,natom_list
               do id2 = id1,natom_list
                  do etype = 0,3
                     if(etype.eq.1.and.(.not.use_sumtype)) exit
                     if(col_used(etype,id2,id1,atomid_central)) then
                        
                        nc_used = ncol_used(etype,id2,id1,atomid_central)
                        id = id1
                        if(id1.eq.0) then
                           write(fileno,*) 'cpair: ',atom_list(atomid_central),atom_list(id2),'central: ',&
                                atom_list(atomid_central), 'type: ',etype
                           write(fileno,*) 'ncols',nc_used
                        else
                           write(fileno,*) 'pair:  ',atom_list(id1),atom_list(id2),'central: ',&
                                atom_list(atomid_central), 'type: ',etype
                           write(fileno,*) 'ncols',nc_used
                        endif
                        do m1 = 1,nc_used
                           cindex = col_index(etype,m1,id2,id1,atomid_central)
                           write(fileno,*) (mnet_wmat1(m2,cindex),m2 = 1,mnet_size(1))
                        end do
                     endif
                  end do
               end do
            end do
            write(fileno, *) 'END'

            do layer = 1,mnet_nlayers
               if(layer.gt.1) then 
                  write(fileno, *) 'WEIGHT | LAYER',layer
                  do m2 = 1,mnet_size(layer) 
                     write(fileno, *) (mnet_w(m2,m1,atomid_central,layer), m1=1,mnet_size(layer-1)) 
                  end do 
               endif
               if(layer.lt.mnet_nlayers) then
                  write(fileno, *) 'BIAS   | LAYER',layer
                  write(fileno, *) (mnet_b(m1,atomid_central,layer),m1=1,mnet_size(layer))
               endif
            end do

         end do

         write(fileno, *) 'END'
         write(fileno,*) 'nfrag',nfrag_id
         do frag_id = 1,nfrag_id
            write(fileno, *) frag_id,mnet_const(frag_id)
         end do
         call flush (fileno)

       end subroutine print_mnet

      subroutine print_cnet(fileno)
         use molneural_data
         implicit none
         integer itype, i, j, n, nc_used,id1,id2,atomid_central,cindex,m1,m2,id,etype
         integer fileno
         
!!$         write(*,*) 'PRINT_CNET'
!!$         
!!$         REWIND(fileno)
!!$         write(fileno, *) '**** CNETWORK ****', ' error: ', error2_charge + error2_c_reg
!!$         write(fileno, *) 'NHIDDEN ',mnet_size(1)
!!$         write(fileno, *) 'RADIAL_CUT  ', radial_cut
!!$         write(fileno, *) 'ARC_CUT  ', arc_cut
!!$         write(fileno, *) 'NCOLS_TOTAL ',ncols_wmat1
!!$         write(fileno, *) 'MNET_ALPHA ',mnet_alpha         
!!$         write(fileno, *)
!!$         write(fileno, *) 'BEGIN DATA'         
!!$
!!$         do atomid_central = 1,natom_list
!!$            write(fileno,*) 'ATOM: ',atom_list(atomid_central)
!!$            write(fileno, *) 'MNET_WMAT1'
!!$            do id1 = 0,natom_list
!!$               do id2 = id1,natom_list
!!$                  do etype = 0,3
!!$                     if(etype.eq.1.and.(.not.use_sumtype)) exit
!!$                     if(col_used(etype,id2,id1,atomid_central)) then
!!$                        
!!$                        nc_used = ncol_used(etype,id2,id1,atomid_central)
!!$                        id = id1
!!$                        if(id1.eq.0) then
!!$                           write(fileno,*) 'cpair: ',atom_list(atomid_central),atom_list(id2),'central: ',&
!!$                                atom_list(atomid_central), 'type: ',etype
!!$                           write(fileno,*) 'ncols',nc_used
!!$                        else
!!$                           write(fileno,*) 'pair:  ',atom_list(id1),atom_list(id2),'central: ',&
!!$                                atom_list(atomid_central), 'type: ',etype
!!$                           write(fileno,*) 'ncols',nc_used
!!$                        endif
!!$                        do m1 = 1,nc_used
!!$                           cindex = col_index(etype,m1,id2,id1,atomid_central)
!!$                           write(fileno,*) (cnet_wmat1(m2,cindex),m2 = 1,mnet_size(1))
!!$                        end do
!!$                     endif
!!$                  end do
!!$               end do
!!$            end do
!!$            write(fileno, *) 'END'
!!$            
!!$            write(fileno, *) 'CNET_WMAT2'
!!$            write(fileno, *) (cnet_wmat2(m2, atomid_central), m2 = 1,mnet_size(1))
!!$            write(fileno, *) 'CNET_BIAS1'
!!$            write(fileno, *) (cnet_b1vec(m2, atomid_central), m2 = 1,mnet_size(1))             
!!$            write(fileno, *) 'CNET_BIAS2'
!!$            write(fileno, *) cnet_b2vec(atomid_central)
!!$         end do
!!$
!!$         call flush(fileno)
      end subroutine print_cnet


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

      subroutine read_control(control_list, param_list, string_list, ncontrol)
         use parse_text
         use molneural_data
         implicit none
         character(len=32), dimension(1024) :: control_list
         real(8), dimension(1024, 8) :: param_list
         integer ncontrol, line, ios, nitems, m, n, nn
         character(len=2048) :: buffer
         character(2048), dimension(1024) :: string_list
         character(len=2048), dimension(64) :: textlist
         real(8) :: y
         character(len=16) :: string
         character(len=2048) :: fname
         character(len=32) :: type
         integer i, j, k
         integer, allocatable, dimension(:) :: seed
         integer :: seed_size
         character(32) :: logfile
         real(8) :: rand

         restart = .false.         
         !      write(*,*) 'READ CONTROL'

         line = 0
         ios = 0
         logfile = ""
         do while (ios .eq. 0)
            read(100, '(A)', iostat=ios) buffer
            if (ios .eq. 0) then

               line = line + 1
               call split_text_to_list(buffer, textlist, nitems)

               select case (upcase(textlist(1)))
               case ('STDOUT')
                  read(textlist(2), *, iostat=ios) logfile
               case ('NNSIZE')
                  write(*,*) 'nitems = ',nitems
                  mnet_nlayers = nitems - 1
                  allocate(mnet_size(nitems-1))
                  do m = 2,nitems
                     read(textlist(m),*) mnet_size(m-1)
                  end do

               case ('RANDOMIZE')
                  read(textlist(2), *, iostat=ios) randval
                  write(*,*) 'RANDOMIZE',randval
                  additive = .false.
               case ('RANDOMIZE+')
                  read(textlist(2), *, iostat=ios) randval
                  write(*,*) 'RANDOMIZE+',randval                  
                  additive = .true.
               case ('SUMTYPE')
                  write(*,*) 'USE SUMTYPE'
                  use_sumtype = .true.
               case ('READ')
                  read(textlist(2), *, iostat=ios) string
                  string = upcase(string)
                  if (string .eq. 'MNET') then
                     read_mnet_flag = .true.
                  else if(string.eq.'CNET') then 
                     read_cnet_flag = .true.
                  end if
               case ('RUN')
                  read(textlist(2), *, iostat=ios) string
                  string = upcase(string)
                  if (string .eq. 'MNET') then
                     run_mnet_flag = .true.
                     run_cnet_flag = .false.
                  else if(string.eq.'CNET') then 
                     run_cnet_flag = .true.
                     run_mnet_flag = .false.
                  endif
               case('SEED')
                  write(*,*) 'SEED'
                  call random_seed(size=seed_size)
                  allocate(seed(seed_size))
                 read(textlist(2),*,iostat = ios) i
                  seed = i
                  call random_seed( put=seed )
                  call random_number(rand)
               case ('ITMAX')
                  read(textlist(2),*,iostat = ios) itmax
                  write(*,*) 'ITMAX = ',itmax
               case('RESTART')
                  read(textlist(2),*,iostat = ios) logfile
                  restart = .true.
                  restart_calc_valid = .true.
                  call system('cp mnet.out mnet.dat')
                  write(*,*) '**************************'
                  write(*,*) '*** RESTART'
                  write(*,*) '**************************'
               case('RADIAL_CUT')
                  read(textlist(2),*,iostat = ios) radial_cut
                  write(*,*) 'RADIAL_CUT',radial_cut
               case('ARC_CUT')
                  read(textlist(2),*,iostat = ios) arc_cut
                  write(*,*) 'ARC_CUT',arc_cut
               case('MNET_ALPHA')
                  read(textlist(2),*,iostat = ios) mnet_alpha
                  write(*,*) 'MNET_ALPHA',mnet_alpha
               end select
               
            end if
         end do

         logfile = trim(logfile)
         if(logfile.ne."") then 
            if(restart) then
               open(5,file = logfile)
               restart_logfile = logfile
            else
               open(6,file = logfile)
            endif
         endif
      end subroutine read_control

      subroutine open_force_energy_files
         use molneural_data
         implicit none
         integer fileno
         character(64) :: filename
         logical file_exists

         fileno = 400
         inquire (fileno, EXIST=file_exists)
         if (file_exists) then
            close (fileno)
         end if
         filename = "force.dat"
         open(fileno, file=filename)
         
         fileno = 700
         inquire (fileno, EXIST=file_exists)
         if (file_exists) then
            close (fileno)
         end if
         filename = "energy.dat"
         open(fileno, file=filename)

         fileno = 800
         inquire (fileno, EXIST=file_exists)
         if (file_exists) then
            close (fileno)
         end if
         filename = "charge.dat"
         open(fileno, file=filename)


       end subroutine open_force_energy_files
    
       subroutine find_order(index1, index0, ndis, order)
         implicit none
         integer, dimension(256, 2) :: index1, index0
         integer, dimension(256) :: order
         integer ndis, j, k

!  do j = 1,ndis
!     write(*,*) 'aaa1',index1(j,1),index1(j,2),'aaa0',index0(j,1),index0(j,2)
!  end do

         do j = 1, ndis
            do k = 1, ndis
               if ((index1(j, 1) .eq. index0(k, 1) .and. index1(j, 2) .eq. index0(k, 2)) &
            & .or. (index1(j, 1) .eq. index0(k, 2) .and. index1(j, 2) .eq. index0(k, 1))) then
                  order(k) = j
               end if
            end do
         end do

      end subroutine find_order


      subroutine open_append_file(fileno, filename)
         implicit none
         integer fileno
         character(*) :: filename
         logical file_exists

         inquire (FILE=filename, EXIST=file_exists)
         if (file_exists) then
            open(fileno, file=filename, status='old', position='append')
         else
            open(fileno, file=filename)
         end if

      end subroutine open_append_file


subroutine allocate_initial
  use molneural_data
  use parse_text
  implicit none
  integer nat,ios,i,ns,itype,j,ntb_max,ntb,nn,nnmax,np_max,np,n,iatom,ndis,k
  character(32) :: string
  character(3) :: atname
  logical nb,match

  allocate(train_force(3,natoms_max,nsteps_total))
  allocate(train_force1(3,natoms_max,nsteps_total))
  allocate(train_energy(nsteps_total))
  allocate(train_energy1(nsteps_total))
  
  allocate(charge(natoms_max,nsteps_total))
  allocate(charge0(natoms_max,nsteps_total))
  allocate(atname_step(natoms_max,nsteps_total))
  allocate(natoms_step(nsteps_total))
  allocate(atomid_step(natoms_max,nsteps_total))

  allocate(coord(3,natoms_max,nsteps_total))
  
  call flush(6)
end subroutine allocate_initial

subroutine check_files_exist
  logical file_exists,found_error

  found_error = .false.

   INQUIRE(FILE="nncontrol.dat", EXIST=file_exists)
   if(.not.file_exists) then
      write(*,*) "ERROR: COULDN'T FIND FILE nncontrol.dat"
      found_error = .true.
   endif

     if(found_error) stop
end subroutine check_files_exist

subroutine calc_variances
  use molneural_data
  implicit none
  real(8) :: fx,fy,fz,x,y,z,f_rms,e_rms,fsum2,q
  real(8) :: energy,e,delta_e
  real(8), dimension(512) :: esum,esum2,eave
  integer nat,iatom,fcount
  integer, dimension(512) :: ecount
  character(3) :: atname
  integer ios,struct_id,frag_id,ns,itype,n,nn
  character(len=2048) :: buffer
  character(len=2048), dimension(64) :: textlist
  integer nitems,mult
  real(8), allocatable, dimension(:) :: qsum,qsum2
  integer, allocatable, dimension(:) :: qcount

  allocate(qsum(natom_list))
  allocate(qsum2(natom_list))
  allocate(qcount(natom_list))
  allocate(var_charge(natom_list))

  qsum = 0.0d0 
  qsum2 = 0.0d0 
  qcount = 0.0d0 


  evfac = 96.487d0
  evfac = 1.0d0 

  rewind(10)
  rewind(11)
  fcount = 0
  ecount = 0 
  qcount = 0 
  esum = 0.0d0
  esum2 = 0.0d0 
  fsum2 = 0.0d0 

  do ns = 1,nsteps_total
     read(10,*,iostat = ios) nat
     if(ios.ne.0) exit
     do iatom = 1,nat
        read(10,*) atname
     end do
     
     read(11,*) energy
     energy = energy * evfac
     energy = energy/dble(nat)
     
     if(nfrag_ns(ns).eq.1) then 
        frag_id = frag_id_ns(ns,1)
        mult = frag_mult_ns(ns,1)
        
        esum(frag_id) = esum(frag_id) + energy * dble(mult)
        ecount(frag_id) = ecount(frag_id) + dble(mult)
     endif
  end do
  do frag_id = 1,nfrag_id
     eave(frag_id) = esum(frag_id)/ecount(frag_id)
  end do
  
  esum = 0.0d0 
  ecount = 0
  rewind(10) 
  rewind(11) 

  !calculate the average energy for each fragment

  do ns = 1,nsteps_total
     read(10,*,iostat = ios) nat
     if(ios.ne.0) exit

     read(11,*) energy
 
     energy = energy * evfac
     energy = energy/dble(nat)

     if(nfrag_ns(ns).eq.1) then 
        frag_id = frag_id_ns(ns,1)
        mult = frag_mult_ns(ns,1)

        delta_e = energy - eave(frag_id)

        esum(frag_id) = esum(frag_id) + delta_e 
        esum2(frag_id) = esum2(frag_id) + (delta_e**2) 
        ecount(frag_id) = ecount(frag_id) + 1
     endif

     do iatom = 1,nat
        read(10,*) atname,x,y,z,fx,fy,fz,q
        call get_atom_id(atname,itype)
        
        fx = fx * evfac
        fy = fy * evfac
        fz = fz * evfac
        fsum2 = fsum2 + fx**2 + fy**2 + fz**2
        fcount = fcount + 3

        qsum(itype) = qsum(itype) + q
        qsum2(itype) = qsum2(itype) + q**2
        qcount(itype) = qcount(itype) + 1
     end do
  end do

  if(.not.restart) write(*,*) 'ecount',ecount(3),esum2(3)/ecount(3),(esum(3)/ecount(3))**2

  do frag_id = 1,nfrag_id
     var_energy(frag_id) = esum2(frag_id)/ecount(frag_id) - (esum(frag_id)/ecount(frag_id))**2
     var_energy(frag_id) = max(var_energy(frag_id),20.0d0)
  end do 

  do n = 1,natom_list
     var_charge(n) = qsum2(n)/qcount(n) - (qsum(n)/qcount(n))**2
     if(.not.restart) write(*,*) 'var_charge: ',atom_list(n),var_charge(n)
  end do

  var_force = fsum2 / fcount
  if(.not.restart) write(*,*) 'fcount = ',fcount

  do frag_id = 1,nfrag_id
     if(.not.restart) then 
        write(*,*) 'frag: ',frag_id_name(frag_id),'var_energy',var_energy(frag_id),'count',ecount(frag_id)
     endif
  end do 
  if(.not.restart) then 
     write(*,*) 'nstruct_id = ',nstruct_id,'nfrag_id',nfrag_id
     write(*,*) 'var_force',var_force
  endif
end subroutine calc_variances

subroutine print_error(step)
   use molneural_data
  implicit none
  integer step
  character(16) :: step_type

  if(minmode.eq.'energy') then 
     step_type = 'ESTEP'
  else if(minmode.eq.'all')  then 
     if(net_type.eq.'mnet') then 
        step_type = 'MSTEP'
     else if(net_type.eq.'cnet') then 
        step_type = 'CSTEP'
     endif
  else
     step_type = 'ASTEP'
  endif


  if(net_type.eq.'mnet') then 
     write(*,*) step_type,step,&
          'TOTAL_ERROR ',real(error_train),&
          'TRAIN_ERROR ',real(error_train-m_reg_error_train),&
          'ENERGY ',real(1.0d0 - error2_energy), &
          'FORCE ', real(1.0d0 - error2_force),&
          'M_REG ',real(m_reg_error_train), & 
          'VALID_ERROR ',real(error_val), &
          'ENERGY ',real(1.0d0 - error2_energy_val),'FORCE ', real(1.0d0 - error2_force_val)
  else if(net_type.eq.'cnet') then 
     write(*,*) step_type,step,&
          'TOTAL_ERROR ',real(error2_charge + error2_c_reg),&
          'CHARGE ',real(1.0d0 - error2_charge),&
          'C_REG ',real(error2_c_reg), & 
          'VALID ERROR ',real(error2_charge_val), & 
          'CHARGE ',real(1.0d0 - error2_charge_val)
  endif
  call flush(6)
end subroutine print_error


 subroutine sort_string_list(string_list,string_index_list,string_list_sorted,string_index_list_sorted,&
      atomid_list,range,nrange,ilength)
  implicit none
  character(3), dimension(512) :: string_list,string_list_sorted
  integer, dimension(512) :: string_index_list,string_index_list_sorted
  character(3), allocatable, dimension(:) :: slist
  integer, allocatable, dimension(:) :: nlist
  integer istart,ilength,ilength0,nrange
  integer, dimension(0:128) :: range
  integer i,j,atomid
  integer, dimension(128) :: atomid_list
  !sorts a string_list of length ilength starting from position istart  

  istart = 2
  ilength0 = ilength-istart+1
  allocate(slist(ilength0))
  allocate(nlist(ilength0))

  j = 0 
  do i = istart,ilength
     j = j + 1
     slist(j) = string_list(i)
     nlist(j) = string_index_list(i)
  end do 

  call piksr2(ilength0,slist,nlist)

  !need to also calculate the ranges

  nrange = 1
  range(1) = istart
  string_index_list_sorted(1) = string_index_list(1)
  string_list_sorted(1) = string_list(1)
  range(0) = 1
  atomid_list = 0 
  atomid_list(1) = 0
  if(ilength0.eq.0) return
  call get_atom_id(slist(1),atomid)
  atomid_list(2) = atomid

  j = 0 
  do i = istart,ilength
     j = j + 1
     string_list_sorted(i) = slist(j)
     string_index_list_sorted(i) = nlist(j)
     if(j.gt.1) then 
        if(slist(j).ne.slist(j-1)) then 
           nrange = nrange + 1
           range(nrange) = i
           call get_atom_id(slist(j),atomid)

           atomid_list(nrange+1) = atomid
        endif
     endif
  end do
  nrange = nrange + 1
  range(nrange) = ilength+1

  
!  write(*,*) '***************'

!  write(*,*) 'string:  ',(string_list_sorted(j),j=1,ilength)
!  write(*,*) 'indices',(string_index_list_sorted(j),j=1,ilength)
!  write(*,*) 'nrange',nrange
!  write(*,*) 'range',(range(j),j=1,nrange)

!  write(*,*) '***************'

end subroutine sort_string_list


SUBROUTINE piksr2(n,arr,brr)
  INTEGER n
  character(3) :: arr(n)
  integer :: brr(n)
!     Sorts an array arr(1:n) into ascending numerical order, by
!     straight insertion, while making the corresponding rearrangement
!     of the array brr(1:n).
  INTEGER i,j
  character(3) a
  integer b
  do 12 j=2,n               ! Pick out each element in turn.
     a=arr(j)
     b=brr(j)
     do 11 i=j-1,1,-1       ! Look for the place to insert it.
        if(arr(i).le.a)goto 10
        arr(i+1)=arr(i)
        brr(i+1)=brr(i)
11   enddo
     i=0
10   arr(i+1)=a             ! Insert it.
     brr(i+1)=b
12 enddo
  return
END SUBROUTINE piksr2

subroutine eigen0(xvec,couple,nmat,evalue,evector,sort_direction)
  use nr_mod
  implicit none
  real(8), dimension(nmat,nmat) :: mat,evector
  real(8), dimension(nmat) :: evalue,xvec
  integer nmat,nrot,n
  real(8) :: couple
  character(16) :: sort_direction

  mat = couple
  do n = 1,nmat
     mat(n,n) = xvec(n)
  end do

  call jacobi(mat,nmat,nmat,evalue,evector,nrot)
  if(sort_direction.eq."ascending") then 
     call eigsrt(evalue,evector,nmat,nmat)
  else
     call eigsrt_descending(evalue,evector,nmat,nmat)
  endif

  
end subroutine eigen0


subroutine eigen(xvec,couple,nmat,evalue,evector,sort_direction)
  use nr_mod
  implicit none
  real(8), dimension(nmat,nmat) :: mat,evector
  real(8), dimension(nmat) :: xvec,evalue
  real(8), allocatable, dimension(:) :: work
  integer nmat,nrot,n
  real(8) :: couple
  integer lwork, info
  character jobz, uplo 
  character(16) :: sort_direction

  mat = couple
  do n = 1,nmat
     mat(n,n) = xvec(n)
  end do

  allocate(work(3 * nmat - 1))

  lwork = 3 * nmat - 1

!  allocate(work(1024))
!  lwork = 1024

  jobz = 'V'
  uplo = 'U'

  call dsyev(jobz, uplo, nmat, mat, nmat, evalue, work, lwork, info)

  if (info .ne.0) then
      write(*, '(a)') ' '
      write(*, '(a,i8)') ' DSYEV returned nonzero INFO= ', info
  else
      if (jobz .eq. 'V') then
          evector = mat
      end if
  end if

!  call jacobi(mat1,nmat,nmat,evalue,evector,nrot)
  if(sort_direction.eq.'descending') then 
     call eigsrt_descending(evalue,evector,nmat,nmat)
  endif

end subroutine eigen


 subroutine allocate_cols
   use molneural_data
   implicit none        
   integer n,n1,n2,n3
   real(8) :: rand
   integer m1,m2

   allocate(col_used(0:3,0:natom_list,0:natom_list,natom_list))
   allocate(ncol_used(0:3,0:natom_list,0:natom_list,natom_list))   

   allocate(col_index(0:3,max_col,0:natom_list,0:natom_list,natom_list))
   allocate(col_atom_index(0:3,max_col,0:natom_list,0:natom_list,natom_list))
   allocate(ncols_atom_wmat1(natom_list))  

 end subroutine allocate_cols

 subroutine get_atom_id(atname,atomid)
   use molneural_data
   implicit none
   character(3) :: atname
   integer atomid,n

   atomid = 0 
   do n = 1,natom_list
      if(atname.eq.atom_list(n)) then
         atomid = n
         return
      endif
   end do
   
 end subroutine get_atom_id


 
 subroutine find_used_cols
   use molneural_data   
   use parse_text
   implicit none
   integer ns,iatom,j,itype,k,iat,jat,natoms0,nrange,atomid_central,n_types,n_types0
   integer n1_start,n1_end,n2_start,n2_end,n1,n2,ipair,npair,k1,k2,m1,m2
   integer iatomid,jatomid
   integer minid,maxid
   integer :: pair_atno(2,512) 
   integer, dimension(512) :: atindex
   integer, dimension(128) :: atomid_list   
   integer, dimension(0:128) :: range
   integer, dimension(512) :: atindex_sorted
   character(3), dimension(512) :: atname_sorted
   character(3), dimension(512) :: atname   
   integer, dimension(256) :: npairs_of_type
   integer, dimension(512,0:2) :: pairlist
   integer nc,etype,itype0,ndis_central,ndis_non_central,ndis,ndis2
   integer id,nn,cindex,ios
   integer, dimension(4) :: et   
   integer nplus,nminus
   real(8) :: rand
   real(8) :: xdif,ydif,zdif,rdis,rdis2,arc_cut2,qj
   character(32) :: string
   logical inside_cutoff   

   col_used = .false.
   ncol_used = 0
   col_index = 0 
   
   ncols_wmat1 = 0 
   arc_cut2 = arc_cut**2

   ! ncols_atom_wmat1 indexes the columns for each atom
   ncols_atom_wmat1 = 0 

   !charges
   do itype = 1,natom_list
      ncols_wmat1 = ncols_wmat1 + 1      
      charge_col_index(itype) = ncols_wmat1
   end do

   do ns = 1,nsteps_total

      do iatom = 1,natoms_step(ns)
            natoms0 = nbonded(iatom,ns) + 1
            j = 0 
            do k = 0,nbonded(iatom,ns)
               iat = atlist(k,iatom,ns)
               j = j + 1
               atname(j) = atname_step(iat,ns)
               atindex(j) = iat
            end do! do k = 0,ntypebonded(itype)

            call sort_string_list(atname,atindex,atname_sorted,atindex_sorted,&
                 atomid_list,range,nrange,natoms0)
            call get_atom_id(atname_sorted(1),atomid_central)

            npair = 0             
            n_types = 0
            n_types0 = 0 

            npairs_of_type = 0 
            ndis = 0 
            ndis_central = 0 
            ndis_non_central = 0 

            do n1 = 1,nrange
               n1_start = range(n1-1)
               n1_end = range(n1) - 1
               do n2 = n1,nrange
                  n2_start = range(n2-1)
                  n2_end = range(n2) - 1
                  
                  if(n1.eq.n2) then 
                     !same type
                     !check to see if there is more than zero pairs
                     if(n1_end-n1_start.gt.0) then 
                        n_types = n_types + 1
                        ipair = 0 
                        do k1 = 1,n1_end - n1_start + 1
                           do k2 = k1+1,n1_end - n1_start + 1
                              npair = npair + 1
                              ipair = ipair + 1

                           iat = atindex_sorted(n1_start + k1 - 1)
                           jat = atindex_sorted(n2_start + k2 - 1)

                           xdif = coord(1,jat,ns) - coord(1,iat,ns)
                           ydif = coord(2,jat,ns) - coord(2,iat,ns)
                           zdif = coord(3,jat,ns) - coord(3,iat,ns)
                           rdis2 = xdif**2 + ydif**2 + zdif**2 

                           inside_cutoff = .true.
                           if(rdis2.gt.arc_cut2) inside_cutoff = .false.
                           if(inside_cutoff) ndis_non_central = ndis_non_central + 1

                           if(inside_cutoff) then 
                              rdis = dsqrt(rdis2)
                              
                              ndis = ndis + 1
                              npairs_of_type(n_types) = npairs_of_type(n_types) + 1
                              pair_atno(1,ndis) = iat
                              pair_atno(2,ndis) = jat

                           endif


                           end do
                        end do
                     endif
                  else
                     !                 i.e. n1 .ne. n2
                     n_types = n_types + 1
                     if(n1.eq.1) n_types0 = n_types0 + 1
                     
                     ipair = 0 
                     do k1 = 1,n1_end - n1_start + 1
                        do k2 = 1,n2_end - n2_start + 1
                           npair = npair + 1
                           ipair = ipair + 1
                        iat = atindex_sorted(n1_start+k1 - 1)
                        jat = atindex_sorted(n2_start +k2 - 1)
                        
                        xdif = coord(1,jat,ns) - coord(1,iat,ns)
                        ydif = coord(2,jat,ns) - coord(2,iat,ns)
                        zdif = coord(3,jat,ns) - coord(3,iat,ns)
                        rdis2 = xdif**2 + ydif**2 + zdif**2

                        inside_cutoff = .true.
                        if(n1.eq.1) then 
                           !central atom
                           ndis_central = ndis_central + 1
                        else
                           if(rdis2.gt.arc_cut2) inside_cutoff = .false.
                           if(inside_cutoff) ndis_non_central = ndis_non_central + 1
                        endif
                        if(inside_cutoff) then 
                           rdis = dsqrt(rdis2)
                           ndis = ndis + 1
                           npairs_of_type(n_types) = npairs_of_type(n_types) + 1
                           pair_atno(1,ndis) = iat
                           pair_atno(2,ndis) = jat

                        endif

                        end do
                     end do
                  endif
                  
            end do
         end do


         ndis2 = ndis
         if(use_sumtype) ndis2 = ndis2 + ndis_non_central
         
         npair = 0 
         n_types = 0 
         do n1 = 1,nrange
            n1_start = range(n1-1)
            n1_end = range(n1) - 1
            iatomid = atomid_list(n1)
            do n2 = n1,nrange
               jatomid = atomid_list(n2)
               if(n1.eq.n2.and.n1_start.eq.n1_end) cycle
               n_types = n_types + 1

               maxid = max(iatomid,jatomid)
               minid = min(iatomid,jatomid)

               nplus = 0 
               nminus = 0 

               do k = 1,npairs_of_type(n_types)
                  npair = npair + 1

                  pairlist(npair,0) = k
                  pairlist(npair,1) = maxid
                  pairlist(npair,2) = minid

                  if(n_types.le.n_types0) then
                     !i.e. it's a central type
                     jat = pair_atno(2,npair)
                     qj = charge(jat,ns)
                     if(qj.eq.0) qj = 1.d-10

                     if(qj.ge.0) then 
                        nplus = nplus + 1
                     else
                        nminus = nminus + 1
                     endif

                     !central types = r,qrplus,qrminus
                     et = (/1,0,1,1/)
                  else
                     !non central types = r,sum
                     et = (/1,1,0,0/)
                  endif
                  
                  do etype = 0,3
                     if(et(etype+1).eq.0) cycle                     
                     col_used(etype,maxid,minid,atomid_central) = .true.
                     
                     if(etype.eq.0.or.etype.eq.1) then 
                        cindex = col_index(etype,k,maxid,minid,atomid_central)                     
                     else if(etype.eq.2.and.qj.ge.0) then 
                        cindex = col_index(etype,nplus,maxid,minid,atomid_central)                     
                     else if(etype.eq.3.and.qj.lt.0) then 
                        cindex = col_index(etype,nminus,maxid,minid,atomid_central)                     
                     endif

                     if(cindex.eq.0) then
                        ncols_wmat1 = ncols_wmat1 + 1
                        ncols_atom_wmat1(atomid_central) = ncols_atom_wmat1(atomid_central) + 1
                        if(k.gt.max_col) then
                           write(*,*) 'ERROR: max_col exceeded. Please set to larger than ',max_col
                           stop
                        endif
                        cindex = ncols_wmat1

                        if(etype.eq.0.or.etype.eq.1) then 
                           col_index(etype,k,maxid,minid,atomid_central) = cindex
                           col_atom_index(etype,k,maxid,minid,atomid_central) = ncols_atom_wmat1(atomid_central)   
                           if(k.gt.ncol_used(etype,maxid,minid,atomid_central)) then
                              ncol_used(etype,maxid,minid,atomid_central) = k
                           endif
                           
                        else if(etype.eq.2.and.qj.ge.0) then 
                           col_index(etype,nplus,maxid,minid,atomid_central) = cindex
                           col_atom_index(etype,nplus,maxid,minid,atomid_central) = ncols_atom_wmat1(atomid_central)   
                           if(nplus.gt.ncol_used(etype,maxid,minid,atomid_central)) then
                              ncol_used(etype,maxid,minid,atomid_central) = nplus
                           endif

                        else if(etype.eq.3.and.qj.lt.0) then 
                           col_index(etype,nminus,maxid,minid,atomid_central) = cindex
                           col_atom_index(etype,nminus,maxid,minid,atomid_central) = ncols_atom_wmat1(atomid_central)   
                           if(nminus.gt.ncol_used(etype,maxid,minid,atomid_central)) then
                              ncol_used(etype,maxid,minid,atomid_central) = nminus
                           endif
                        endif
                     endif
                  end do
               end do
               
            end do
         end do

         
      end do
      end do


      allocate(mnet_wmat1(mnet_size(1),ncols_wmat1))
      allocate(mnet_de2_dw1(mnet_size(1),ncols_wmat1))

      allocate(cnet_wmat1(mnet_size(1),ncols_wmat1))
      allocate(cnet_de2_dw1(mnet_size(1),ncols_wmat1))

      nn = 0 
      ncols_atom_max = 0 
      do id = 1,natom_list
         nn = nn + ncols_atom_wmat1(id)
         if(ncols_atom_wmat1(id).gt.ncols_atom_max) ncols_atom_max = ncols_atom_wmat1(id)
      end do
      
      if(.not.restart) then 
         write(*,*) 'ncols_wmat1 = ',ncols_wmat1
         !      write(*,*) 'ncols_atom_wmat1 = ',(ncols_atom_wmat1(id),id = 1,natom_list)
         !      write(*,*) 'ncols_atom_max = ',ncols_atom_max
      endif
    end subroutine find_used_cols

    subroutine find_atlist
      use molneural_data
      implicit none      
      integer iatom,jatom,ns
      real(8) :: xdif,ydif,zdif,rdis
      integer ni,nj,n,j

      allocate(nbonded(natoms_max,nsteps_total))
      allocate(atlist(0:256,natoms_max,nsteps_total))
      
      nbonded = 0 
      
      do ns = 1,nsteps_total
         do iatom = 1,natoms_step(ns)
            atlist(0,iatom,ns) = iatom
            do jatom = iatom+1,natoms_step(ns)
               xdif = coord(1,jatom,ns) - coord(1,iatom,ns)
               ydif = coord(2,jatom,ns) - coord(2,iatom,ns)
               zdif = coord(3,jatom,ns) - coord(3,iatom,ns)
               rdis = dsqrt(xdif**2 + ydif**2 + zdif**2)                

               if(rdis.le.radial_cut) then
                  nbonded(iatom,ns) = nbonded(iatom,ns) + 1
                  nbonded(jatom,ns) = nbonded(jatom,ns) + 1
                  ni = nbonded(iatom,ns)
                  nj = nbonded(jatom,ns)
                  atlist(ni,iatom,ns) = jatom
                  atlist(nj,jatom,ns) = iatom                  
               endif
            end do
         end do
         
!         write(*,*) 

      end do !do ns = 1,nsteps_total
!      stop
      
    end subroutine find_atlist

    subroutine read_natoms
      use molneural_data
      implicit none
      integer nat,ns,ios,iatom,n,id1,id2,nn,j
      character(3) :: atname
      logical match

      natoms_max = 0 

      if(.not.restart) write(*,*) 'READ_NATOMS'
      
      rewind(10)
      natom_list = 0 
      ns = 0 
      do while(.true.)
         read(10,*,iostat = ios) nat
         if(ios.ne.0) exit

         ns = ns + 1
            if(nat.gt.natoms_max) natoms_max = nat
         do iatom = 1,nat
            read(10,*) atname

            match = .true.
            do n = 1,natom_list
               if(atname.eq.atom_list(n)) then
                  match = .false.
                  exit
               endif
            end do
            if(match) then
               natom_list = natom_list + 1
               atom_list(natom_list) = atname
            endif
            

         end do
      end do
      
      rewind(10)      
      if(.not.restart) then 
         write(*,*) 'natom_list = ',natom_list
         write(*,*) 'atom_list = ',(atom_list(j),j=1,natom_list)
      endif

    end subroutine read_natoms

    subroutine run
      use molneural_data
      use networks_mod
      
      
      trainmode = 'ALTERNATE'
      
      call print_results
      
      !      call test_coulomb
!      call test_mnet
      !stop
!        call test_mnet_forces
!        stop
      
      if(.not.restart) then 
         nminstep = 0

         minmode = 'energy'
         call minimize
         nminsteps = 0 
      endif

         minmode = 'all'
         call minimize

      if(net_type.eq.'mnet') then 
         call run_train_mnet
      else if(net_type.eq.'cnet') then 
         call run_train_cnet
      endif
         
         nminstep = nminstep + 1
         call print_results
      
    end subroutine run
    
subroutine GETBFUNCS(eps,eps_sqrtpii,rr,rri,bfunc,nmax)
!---------------------------------------------------
!     calculate the B functions from W. Smith's 
  !     Ewald Sum revisited paper
  !---------------------------------------------------
  implicit none
  real(8),intent(in) :: rr,rri
  real(8) :: eps,efac,rr2i,pfac
  real(8) :: pow,rr_eps,eps_sqrtpii
  real(8), dimension(0:10),intent(out) :: bfunc
  integer :: i,nmax
  
  rr2i = rri * rri
  rr_eps = rr * eps
  
  bfunc(0) = derfc(rr_eps) * rri
  efac = dexp(-rr_eps**2) * eps_sqrtpii
  
  pow = 1.d0
  pfac = 2.d0 * eps*eps
  do i = 1,nmax
     pow = pow * pfac
     bfunc(i) = rr2i*(dble(2*i-1)*bfunc(i-1)+pow * efac)
  end do
  
end subroutine getbfuncs

subroutine GETPFUNCS(rri,rr2i,pfunc,nmax)
!---------------------------------------------------
!     calculate the P functions - the point charge analgoue to the B functions
!---------------------------------------------------
  implicit none
  real(8), intent(in) :: rri,rr2i
  integer, intent(in) :: nmax
  integer  n
  real(8), dimension(0:10),intent(out) :: pfunc
  
  pfunc(0) = rri 
  
  do n = 0,nmax-1
     pfunc(n+1) = (2*n+1) * pfunc(n) * rr2i 
  end do
  
end subroutine getpfuncs


subroutine calc_coulomb
  use molneural_data
  implicit none

  integer ns,iatom,jatom
  real(8) :: xdif,ydif,zdif,rdis,ri,r2i
  real(8), dimension(0:10) :: erfc_func,gfunc,tfunc,pfunc
  real(8) :: eps_damp,eps_damp_sqrtpii,damp_width
  real(8) :: ffix,ffiy,ffiz
  real(8) :: ggg(0:4)
  real(8) :: cij,ci,cj
  real(8) :: u
  integer nmax

  damp_width = 0.5d0 
  pi = 4.0d0 * datan(1.0d0)
  
  eps_damp = 1.0d0/(dsqrt(2.0d0) * damp_width)
  eps_damp_sqrtpii = 1.0d0/(eps_damp*dsqrt(pi))
  
  nmax = 1

  do ns = 1,nsteps_total
     u = 0.0d0 
     do iatom = 1,natoms_step(ns)
        ci = charge(iatom,ns)
        do jatom = iatom+1,natoms_step(ns)
           cj = charge(jatom,ns)

           cij = ci*cj
           
           xdif = coord(1,jatom,ns) - coord(1,iatom,ns)
           ydif = coord(2,jatom,ns) - coord(2,iatom,ns)
           zdif = coord(3,jatom,ns) - coord(3,iatom,ns)
           rdis = dsqrt(xdif**2 + ydif**2 + zdif**2) 
           ri = 1.0d0/rdis
           r2i = ri*ri

           call getpfuncs(ri,r2i,pfunc,nmax)
           call getbfuncs(eps_damp,eps_damp_sqrtpii,rdis,ri,tfunc,nmax)
           gfunc = pfunc - tfunc

           ffix = cij * gfunc(1) * xdif*fac(2)*fac(3)
           ffiy = cij * gfunc(1) * ydif*fac(2)*fac(3)
           ffiz = cij * gfunc(1) * zdif*fac(2)*fac(3)

           u = u + cij * gfunc(0)*fac(2)*fac(3)
           
           train_force(1,iatom,ns) = train_force(1,iatom,ns) + ffix
           train_force(2,iatom,ns) = train_force(2,iatom,ns) + ffiy
           train_force(3,iatom,ns) = train_force(3,iatom,ns) + ffiz

           train_force(1,jatom,ns) = train_force(1,jatom,ns) - ffix
           train_force(2,jatom,ns) = train_force(2,jatom,ns) - ffiy
           train_force(3,jatom,ns) = train_force(3,jatom,ns) - ffiz
           
        end do
     end do

     train_energy(ns) = train_energy(ns) - u
     
  end do! do ns = 1,nsteps_total
  

end subroutine calc_coulomb


subroutine set_units()
  use molneural_data
  implicit none

  !     sets all the units used in the simulation
  
  pi = 4.0d0 * atan(1.0d0)
  sqrtpi = dsqrt(pi)
  eps0 = 8.854187817d-12
  avsno = 6.022d23
  
  pifac = 180.0d0/(4.0d0*atan(1.0d0))
  
  !     UNITS
  !----------------------------------------------
  !     (1)  TIME
  !     (2)  DISTANCE
  !     (3)  MASS
  !     (4)  CHARGE OF ONE ELECTRON
  !     (5)  TEMPERATURE
  !----------------------------------------------
  
  unit(1) = 1.d-15          !1 fs
  unit(2) = 1.d-10          !1 Angstrom
  unit(3) = 1.67262178d-27  !proton mass
  unit(4) = 1.60217646d-19  !charge of one electron
  unit(5) = 1.0             !Kelvin
  
  !     speed of light in internal units
  lightspeed = 299792458.0d0 * unit(1) / unit(2)
  
  !----------------------------------------------
  !     derived units/factors
  !----------------------------------------------
  !     (1) conversion of ENERGY from internal to SI
  !     (2) e^2 / 4 pi * eps0 in internal units
  !     (3) conversion of internal ENERGY to KJ/MOL
  !     (4) conversion of internal ENERGY to KCAL/MOL
  !     (5) conversion of Bohr to internal units (Angstroms)
  !     (6) conversion of FIELD from internal to Volts/Angstrom
  !     (7) conversion of PRESSURE from internal to bar
  !     (8) conversion of density from internal to g/cm^3
  !     (9) conversion of cm-1 (frequency) to internal units
  !     (10) conversion of internal units to Debye
  !----------------------------------------------
  
  fac(1) = unit(3) * unit(2)**2 / unit(1)**2
  fac(2) = (unit(4) * unit(4) / (fac(1) * unit(2))) * 1.d0/(4.0d0 * pi * eps0)
  fac(3) = fac(1) * (avsno / 1000.0d0)
  fac(4) = (fac(3) / 4.184d0)
  fac(5) = 0.529177d0
  fac(6) = fac(1) / unit(4)
  fac(7) = (fac(1) / unit(2)**3) / (1.d+5)
  fac(8) = (unit(3) / unit(2)**3)/1000.0d0 
  fac(9) = lightspeed * 100.0d0 * unit(2) * 2.0d0 * pi
  !      fac(9) = 29.9792458 * 1.d+9 * unit(1) * 2.0d0 * pi
  fac(10) = 1.0d0 / 0.20819434d0
  
  boltz = (1.38064852d-23) * unit(5)/fac(1)
  hbar = 1.0545718d-34 /(unit(1)*fac(1))
  
end subroutine set_units

subroutine shift_energies
  use molneural_data
  implicit none
  integer ns
  real(8) :: emin

  emin = 1.d+23
  do ns = 1,nsteps_total
     if(train_energy(ns).le.emin) emin = train_energy(ns)
  end do
  do ns = 1,nsteps_total
!     train_energy(ns) = train_energy(ns) - emin
  end do
     
end subroutine shift_energies

subroutine test_coulomb
  use molneural_data
  implicit none
  real(8) :: u0,u1
  integer ns,iatom,ivec
  real(8) :: tiny
  real(8) :: numerical,analytical
  
  tiny = 1.0d-6
  ns = 1
  

  do ns = 1,nsteps_total

     train_energy = 0.0d0
     train_force = 0.0d0 
     call calc_coulomb
     u0 = train_energy(ns)
     
     do iatom = 1,natoms_step(ns)
        do ivec = 1,3
           coord(ivec,iatom,ns) = coord(ivec,iatom,ns) + tiny
           train_energy = 0.0d0
           train_force = 0.0d0 
           call calc_coulomb
           u1 = train_energy(ns)
           coord(ivec,iatom,ns) = coord(ivec,iatom,ns) - tiny        
           numerical = (u0-u1)/tiny
           analytical = train_force(ivec,iatom,ns)
           write(*,*) ns,iatom,ivec,numerical,analytical,numerical/analytical
        end do
     end do
  end do
     
end subroutine test_coulomb

subroutine print_energies_and_forces
   use molneural_data
   implicit none
   integer iatom,ivec,ns
   real(8) :: nfac
   
   write(700, *) "## Train energy GROUND TRUTH --- PREDICT "
   write(400, *) "## Train force GROUND TRUTH --- PREDICT "
   do ns = 1, nsteps_total
      if (train_step(ns)) then
         nfac = 1.0d0 / dble(natoms_step(ns))
         write(700, *) train_energy(ns)*nfac, train_energy1(ns)*nfac, &
              struct_id_name(struct_id_ns(ns)),ns              
         do iatom = 1, natoms_step(ns)
            do ivec = 1, 3
               write(400, *) train_force(ivec, iatom, ns), & 
                    train_force1(ivec, iatom, ns), &
                    ns, iatom, ivec
            end do
         end do
      end if
   end do
   write(700, *) ""
   write(700, *) "## Valid energy GROUND TRUTH --- PREDICTION "
   write(400, *) ""
   write(400, *) "## Valid force GROUND TRUTH --- PREDICT "
   do ns = 1, nsteps_total
      if (.not.train_step(ns)) then
         nfac = 1.0d0 / dble(natoms_step(ns))
         write(700, *) train_energy(ns)*nfac, train_energy1(ns)*nfac, &
              struct_id_name(struct_id_ns(ns)),ns              
         do iatom = 1, natoms_step(ns)
            do ivec = 1, 3
               write(400, *) train_force(ivec, iatom, ns), & 
                    train_force1(ivec, iatom, ns), &
                    ns, iatom, ivec
            end do
         end do
      end if
   end do
   
   call flush (700)
   call flush (400)
   REWIND(700)
   REWIND(400)
   
end subroutine print_energies_and_forces

subroutine print_charges
   use molneural_data
   implicit none
   integer iatom,ns

   write(*,*) 'PRINT CHARGES'
   write(800, *) "## Train charge GROUND TRUTH --- PREDICT "
   do ns = 1,nsteps_total
      if (train_step(ns)) then
         do iatom = 1,natoms_step(ns)
            write(800,*) charge(iatom,ns),charge0(iatom,ns),atname_step(iatom,ns) 
         end do
      endif
   end do

   write(800,*) 
   write(800, *) "## Valid charge GROUND TRUTH --- PREDICT "
   do ns = 1,nsteps_total
      if (.not.train_step(ns)) then
         do iatom = 1,natoms_step(ns)
            write(800,*) charge(iatom,ns),charge0(iatom,ns),atname_step(iatom,ns) 
         end do
      endif
   end do


   call flush(800)
   rewind(800)
end subroutine print_charges

subroutine get_frag_id_from_string(frag_id_string,frag_id)
  use molneural_data
  implicit none
  character(32) :: frag_id_string
  integer mm,frag_id
  logical match

  match = .true. 
  do mm = 1,nfrag_id
     if(frag_id_string.eq.frag_id_name(mm)) then
        match = .false.
        frag_id = mm
        exit
     endif
  end do
  if(match) then
     nfrag_id = nfrag_id + 1
     frag_id = nfrag_id
     frag_id_name(frag_id) = frag_id_string
  endif
  
end subroutine get_frag_id_from_string


subroutine get_struct_id_from_string(struct_id_string,struct_id)
  use molneural_data
  implicit none
  character(32) :: struct_id_string
  integer mm,struct_id
  logical match

  match = .true. 
  do mm = 1,nstruct_id
     if(struct_id_string.eq.struct_id_name(mm)) then
        match = .false.
        struct_id = mm
        exit
     endif
  end do
  if(match) then
     nstruct_id = nstruct_id + 1
     struct_id = nstruct_id
     struct_id_name(struct_id) = struct_id_string
  endif
  
end subroutine get_struct_id_from_string

subroutine find_nsteps_total
  use molneural_data
  implicit none
  character(3) :: aname
  integer nat,ns,iatom,ios
  real(8) :: e
  
  rewind(11)
  rewind(10)
  ns = 0 
  do while(.true.)
     read(11,*,iostat = ios) e
     if(ios.ne.0) exit
     read(10,*,iostat = ios) nat
     if(ios.ne.0) exit

     do iatom = 1,nat
        read(10,*) aname
     end do
     ns = ns + 1
  end do
  nsteps_total = ns

  if(.not.restart) write(*,*) 'nsteps_total = ',nsteps_total
  rewind(11)
  rewind(10)
end subroutine find_nsteps_total

subroutine get_mnet_reg
  use molneural_data
  implicit none
  integer id1,id2,etype,m1,m2,cindex,atomid_central,nc_used
  integer k,itype,layer,count
  real(8) :: lambda_red

  lambda_m = 0.020d0 
  count = 0 


  count = ncols_wmat1*mnet_size(1)

!!$  do layer = 2,mnet_nlayers
!!$     count = count + mnet_size(layer-1)*mnet_size(layer)*natom_list
!!$  end do

  lambda_red = lambda_m/dble(count)

  do cindex = 1,ncols_wmat1
     do m2 = 1,mnet_size(1)
        error2_m_reg = error2_m_reg + 0.5d0 * lambda_red * mnet_wmat1(m2,cindex)**2
        mnet_de2_dw1(m2,cindex) = mnet_de2_dw1(m2,cindex) + lambda_red * mnet_wmat1(m2,cindex)
     end do
  end do

!!$  do layer = 2,mnet_nlayers
!!$     do itype = 1,natom_list
!!$        do m1 = 1,mnet_size(layer-1)
!!$           do m2 = 1,mnet_size(layer)
!!$              error2_m_reg = error2_m_reg + 0.5d0 * lambda_red * mnet_w(m2,m1,itype,layer)**2
!!$              mnet_de2_dw(m2,m1,itype,layer) = mnet_de2_dw(m2,m1,itype,layer) &
!!$                   + lambda_red * mnet_w(m2,m1,itype,layer)              
!!$           end do
!!$        end do
!!$     end do 
!!$  end do

end subroutine get_mnet_reg

subroutine get_cnet_reg
  use molneural_data
  implicit none
  integer id1,id2,etype,m1,m2,cindex,atomid_central,nc_used
  integer k,ncols_wmat1c
  real(8) :: lambda_red

  !the first natom_list indices are reserved for the charge columns, which aren't used for cindex

  ncols_wmat1c = ncols_wmat1 - natom_list

  lambda_m = 0.1d0  
  lambda_red = lambda_m/dble(ncols_wmat1c*mnet_size(1))

  do cindex = natom_list+1,ncols_wmat1
     do m2 = 1,mnet_size(1)
        error2_c_reg = error2_c_reg + 0.5d0 * lambda_red * cnet_wmat1(m2,cindex)**2
        cnet_de2_dw1(m2,cindex) = cnet_de2_dw1(m2,cindex) + lambda_red * cnet_wmat1(m2,cindex)
     end do
  end do

end subroutine get_cnet_reg

   subroutine count_atoms_and_molecules
     use parse_text
     use molneural_data
     implicit none
     integer iatom,ns,itype
     character(3) :: atname
     character(len=2048) :: buffer
     character(len=2048), dimension(64) :: textlist
     character(32) :: train_string
     integer nitems,ios
     real(8) :: e

     charge_count = 0.0d0 
     ncharges_total = 0 
     ncharges_train = 0
     ncharges_valid = 0 

     do ns = 1,nsteps_total

        do iatom = 1,natoms_step(ns)
           atname = atname_step(iatom,ns)
           call get_atom_id(atname,itype)        
           charge_count(itype) = charge_count(itype) + 1
           ncharges_total = ncharges_total + 1
           if(train_step(ns)) then 
              ncharges_train(itype) = ncharges_train(itype) + 1
           else
              ncharges_valid(itype) = ncharges_valid(itype) + 1
           endif
        end do
     end do
     
   end subroutine count_atoms_and_molecules

   subroutine read_efile
     use molneural_data
     use parse_text
     implicit none
     character(len=2048) :: buffer
     character(len=2048), dimension(64) :: textlist
     integer ns,ios,nitems,frag_id,struct_id,nground,n
     real e
     character(32) :: frag_id_string,train_string
     character(32), dimension(32) :: stringlist
     character(32) :: string
     integer, dimension(32) :: multlist
     character(4) :: aaa

     nsteps_train = 0 
     nsteps_val = 0 
     nsteps_train_struct_id = 0
     nsteps_val_struct_id = 0          
     struct_id_count = 0 

     allocate(train_step(nsteps_total))

     nstruct_id = 0
     nfrag_id = 0
     rewind(11)
     do ns = 1,nsteps_total
        read(11, '(A)', iostat=ios) buffer   
        call split_text_to_list(buffer, textlist, nitems)
        call parse_ground(buffer,stringlist,multlist,nground)

        do n = 1,nground
           
           call get_frag_id_from_string(stringlist(n),frag_id)           
           
           !frag_id_ns contains the id numbers of the ground state energies for each fragment.
           !A particular structure can include more than one fragment, n fragments in total. 
           !Similarly, frag_mult_ns contains the multipliers for each fragment. e.g. if a structure 
           !contains 3 of the same fragment, the the multiplier would be 3.
           
           frag_id_ns(ns,n) = frag_id
           frag_mult_ns(ns,n) = multlist(n)
           nfrag_ns(ns) = n

        end do

        read(textlist(1),*) string


        if(nground.eq.1) then 
           string = trim(stringlist(1))
        else
           string = ""
           do n = 1,nground
                 
              if(multlist(n).ne.1) then 
                 if(n.gt.1) string = trim(string)//"+"
                 write(aaa,'(I4)') multlist(n)
                 string = trim(string)//adjustl(trim(aaa))
              endif
              string = trim(string)//adjustl(trim(stringlist(n)))
           end do
        endif



        read(textlist(2),*) train_string
        if(ios.ne.0) exit
        
        if(upcase(train_string).eq.'TRAIN') then
           train_step(ns) = .true.
        else
           train_step(ns) = .false.
        endif


        call get_struct_id_from_string(string,struct_id)           

        struct_id_ns(ns) = struct_id
        struct_id_count(struct_id) = struct_id_count(struct_id) + 1

        if(train_step(ns)) then 
           nsteps_train = nsteps_train + 1
           nsteps_train_struct_id(struct_id) = nsteps_train_struct_id(struct_id) + 1
        else
           nsteps_val = nsteps_val + 1
           nsteps_val_struct_id(struct_id) = nsteps_val_struct_id(struct_id) + 1
        endif
     end do

   end subroutine read_efile

  subroutine parse_ground(buffer,stringlist,multlist,nground)
    use parse_text
    implicit none
    integer ios,nitems,n,nitems2,n2,mult
    character(len=2048) :: buffer,text1,text2
    character(len = 2048), dimension(64) :: textlist,textlist2
    character(32) :: ggg
    character(32), dimension(32) :: stringlist
    integer, dimension(32) :: multlist
    integer nground
    
    call split_text_to_list2(buffer,'|',textlist,nitems)
    buffer = textlist(2)
    
    call split_text_to_list2(buffer,'+',textlist,nitems)

    do n = 1,nitems
       call split_text_to_list2(textlist(n),'*',textlist2,nitems2)
       
       mult = 0 
       if(nitems2.eq.2) then 
          read(textlist2(1),*) mult
          read(textlist2(2),*) ggg
       else if(nitems2.eq.1) then 
          mult = 1
          read(textlist2(1),*) ggg
       endif
       
       stringlist(n) = ggg
       multlist(n) = mult

    end do

    nground = nitems
    
  end subroutine parse_ground

  subroutine restart_log
    use molneural_data
    implicit none
    integer ios
    character(16) :: string

    !    write(*,*) 'GET STEP NO'

    nminsteps = 0
     do while(.true.)
       read(*,*,iostat = ios) string
       if(ios.ne.0) exit
       if(string.eq.'MSTEP') nminsteps = nminsteps + 1
    end do
    write(*,*) 'STEPS: ',nminsteps

    open(6,file = restart_logfile,status='old',access='append')
  end subroutine restart_log
