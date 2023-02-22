module molneural_data
  use netparam
  logical :: readed_atomlist = .false.

!******* MNET *********

  real(8), allocatable, dimension(:,:) :: mnet_wmat2,mnet_wmat3
  real(8), allocatable, dimension(:,:) :: mnet_b1vec,mnet_b2vec,mnet_b3vec
  real(8), allocatable, dimension(:,:) :: mnet_de2_db1,mnet_de2_db2,mnet_de2_db3
  
  real(8), allocatable, dimension(:,:,:) :: mnet_de2_db
  real(8), allocatable, dimension(:,:,:,:) :: mnet_de2_dw

  real(8), allocatable, dimension(:) :: mnet_de2_dconst
  real(8), allocatable, dimension(:) :: mnet_const

  real(8), allocatable, dimension(:,:,:) :: mnet_de2_dw2,mnet_de2_dw3

  real(8), allocatable, dimension(:,:,:) :: cnet_wmat2,cnet_wmat3
  real(8), allocatable, dimension(:,:,:) :: cnet_de2_dw2,cnet_de2_dw3
  real(8), allocatable, dimension(:,:) :: cnet_b1vec
  real(8), allocatable, dimension(:) :: cnet_b2vec
  real(8), allocatable, dimension(:,:) :: cnet_de2_db1,cnet_de2_db2

  real(8) mnet_denergy_dmol_energy
  real(8), allocatable, dimension(:) :: de2_dmol_energy

  integer mnet_n1max,mnet_s1max
  integer :: mnet_s1

  integer, allocatable, dimension(:,:,:) :: type_perm_list
  integer, allocatable, dimension(:) :: type_nperm

!******* RNET *********

  real(8), allocatable, dimension(:,:) :: rnet_wmat1,rnet_de2_dw1
  real(8), allocatable, dimension(:,:) :: rnet_wmat2,rnet_de2_dw2
  real(8), allocatable, dimension(:,:) :: rnet_b1vec,rnet_de2_db1
  real(8), allocatable, dimension(:) :: rnet_b2vec,rnet_de2_db2
  real(8), allocatable, dimension(:) :: rnet_x1vec

  integer rnet_n2max,ndis_max
  integer :: rnet_n1,rnet_n2

  integer, allocatable, dimension(:,:) :: rnet_ndis

  integer mtypes

!****** END RNET ********


  integer nsteps_train,nsteps_val,nsteps_train0,nsteps_total0,nsteps_total
  integer ntypes
  character(32), allocatable, dimension(:,:) :: typelist
  integer, allocatable, dimension(:) :: ntypebonded,ntypebonded1,ntypebonded2

  real(8), allocatable, dimension(:,:,:) :: train_force,train_force0,train_force1
  real(8), allocatable, dimension(:) :: train_energy,train_energy0,train_energy1

  integer, allocatable, dimension(:) :: natoms
  integer, allocatable, dimension(:,:) :: atomtype

  real(8) :: error2_force,error2,error2_energy,error2_m_reg,error2_c_reg, error2_force_val, error2_energy_val
  real(8) :: error2_charge,error2_charge_val
  
  real(8) :: aufac,aufac_f
  real(8) :: error_train,m_reg_error_train,r_reg_error_train
  real(8) :: error_val,  force_error_val,  energy_error_val, min_error_val

  logical print_forces,restart_calc_valid

  integer nminstep
  character(len = 16) :: minmode

  real(8), dimension(8192) :: dismin_b,dismax_b,dismin_nb

  integer rnet_npairtypes,mnet_npairtypes

  logical, allocatable, dimension(:) :: rnet_used,mnet_used
  logical, allocatable, dimension(:,:) :: mnet_used_mol

  integer, allocatable, dimension(:) :: mol_train_count
  integer nmol_types
  character(len = 32), allocatable, dimension(:) :: molname
  integer, allocatable, dimension(:) :: train_mol

  integer natoms_train,natoms_val
  integer :: imol_max

  logical first_read
  real(8), allocatable, dimension(:) :: mol_force_error_train,mol_energy_error_train
  integer, allocatable, dimension(:,:,:) :: atomlist      

  integer, allocatable, dimension(:,:) :: rhist_b,rhist_nb
  integer nbinmax
  real(8) :: rmin_nb,rmax_nb,rmin_b,rmax_b

  integer, allocatable, dimension(:,:) :: rnet_train_npairtype

  real(8), allocatable, dimension(:) :: mol_energy

  integer molmax
  integer neighbor
  integer imol0,ii0

  real(8) :: evfac,evfac_f

  integer nprint,neighbor0,natoms_max
  integer itmax,nminsteps
  logical mnet_allocated,rnet_allocated

  character(32) :: mnet_activation_type,rnet_activation_type

  real(8) :: lambda_m,lambda_r      
  character(16) :: trainmode
  
  real(8) :: var_force
  real(8), dimension(512) :: var_energy

  integer, dimension(128,128) :: bondtable
  real(8), allocatable, dimension(:,:,:) :: coord
  character(32), dimension(128) :: mol_list

  integer nmol_list

  character(3), dimension(128,128) :: atomic_name

  integer natom_list
  character(3), dimension(128) :: atom_list

  logical, allocatable, dimension(:,:,:,:) :: col_used
  integer, allocatable, dimension(:,:,:,:) :: ncol_used  

  integer, allocatable, dimension(:,:,:,:,:) :: col_index,col_atom_index
  real(8), allocatable, dimension(:,:) :: mnet_wmat1,mnet_de2_dw1
  real(8), allocatable, dimension(:,:) :: cnet_wmat1,cnet_de2_dw1
  real(8), allocatable, dimension(:,:,:,:) :: mnet_dforce_dw1  

  !dimension by the number of atom-types, H, C, O, etc.
  type(net_param), allocatable, dimension(:) :: mnet
  type(net_param), allocatable, dimension(:) :: cnet
  
  real(8), allocatable, dimension(:,:,:,:) :: mnet_w
  real(8), allocatable, dimension(:,:,:) :: mnet_b

  integer ncols_wmat1
  integer, allocatable, dimension(:) :: ncols_atom_wmat1
  integer ncols_atom_max
  real(8) :: radial_cut,rnet_rcut,arc_cut
  integer, allocatable, dimension(:,:) :: nbonded
  integer, allocatable, dimension(:,:,:) :: atlist
  real(8) :: rdis_store
  integer, allocatable, dimension(:,:) :: pairtype_lookup
  integer, allocatable, dimension(:,:) :: atomid_step
  integer, allocatable, dimension(:,:) :: id_lookup
  integer, allocatable, dimension(:) :: natoms_step
  character(3), allocatable, dimension(:,:) :: atname_step
  real(8), dimension(10) :: fac,unit
  real(8) :: avsno,boltz,eps0,pi,sqrtpi,hbar,lightspeed,pifac
  real(8), allocatable, dimension(:,:) :: charge,charge0
  integer max_col
  logical read_mnet_flag,read_cnet_flag
  real(8) :: randval
  logical additive
  logical use_sumtype
  integer nstruct_id,nfrag_id
  integer, dimension(2048) :: struct_id_ns
  integer, dimension(2048) :: struct_id_count
  integer, dimension(2048) :: nsteps_train_struct_id,nsteps_val_struct_id
  integer, dimension(128) :: charge_col_index
  character(16) :: net_type

  logical run_mnet_flag,run_cnet_flag  
  
  character(32), dimension(2048) :: struct_id_name,frag_id_name
  logical, allocatable, dimension(:) :: train_step
  
  logical calc_valid
  real(8) :: mnet_alpha,cnet_alpha
  real(8), allocatable, dimension(:) :: var_charge
  integer ncharges_total
  integer, dimension(1024) :: charge_count,ncharges_train,ncharges_valid
  
  integer, dimension(8192) :: nfrag_ns
  integer, dimension(8192,16) :: frag_id_ns,frag_mult_ns
  logical restart
  character(32) :: restart_logfile

   integer, allocatable, dimension(:) :: mnet_size,cnet_size
   integer mnet_nlayers,cnet_nlayers

   integer mnet_size_max
   logical, dimension(8192) :: use_type
   logical calc_forces

end module molneural_data
