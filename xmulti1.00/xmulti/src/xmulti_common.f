      module common_data
      real(8), dimension(:,:,:,:), allocatable :: site_coord0
      real(8), dimension(:,:), allocatable :: lab_charge
      real(8), dimension(:,:,:), allocatable :: lab_dipole
      real(8), dimension(:,:,:,:), allocatable :: lab_quad
      real(8), dimension(:,:,:,:,:), allocatable :: lab_oct
      real(8), dimension(:,:,:), allocatable :: lab_coord,lab_relative_coord
      character(len = 3), dimension(:,:), allocatable :: site_name
      integer, dimension(:), allocatable :: mol_natoms
      integer, dimension(:), allocatable :: mol_nmassatoms
      integer, dimension(:), allocatable :: mol_type
      real(8), dimension(:,:), allocatable :: mol_euler,mol_com
      integer nmol,natoms,nmassatoms
      real(8) :: utotal,upair,ucharge,udipole,uquad,uoct,uconformer
      real(8), dimension(:,:,:), allocatable :: force
      real(8), dimension(:,:), allocatable :: mass
      real(8), dimension(:,:), allocatable :: forcemol,fracforcemol,torquemol
      real(8), dimension(:,:), allocatable :: dudeulermol
      character(len = 16), dimension(:), allocatable :: molecule_name_list

      real(8), dimension(:,:,:), allocatable :: site_charge0
      real(8), dimension(:,:,:,:), allocatable :: site_dipole0
      real(8), dimension(:,:,:,:,:), allocatable :: site_quad0
      real(8), dimension(:,:,:,:,:,:), allocatable :: site_oct0

      real(8), dimension(:,:), allocatable :: cfield
      real(8), dimension(:,:,:), allocatable :: dfield
      real(8), dimension(:,:,:,:), allocatable :: qfield
      real(8), dimension(:,:,:,:,:), allocatable :: ofield

      real(8), dimension(:,:,:), allocatable :: lab_spherical2,lab_spherical3
      real(8), dimension(:,:,:,:), allocatable :: lab_spherical2_1

      real(8), dimension(:,:,:), allocatable :: fieldtens1,fieldtens2,fieldtens3
      real(8), dimension(:,:), allocatable :: fieldtens0
      real(8), dimension(:,:,:,:), allocatable :: fieldtens1_1,fieldtens2_1

      integer nmoltypes
      real(8), dimension(10) :: unit,fac
      real(8) :: lightspeed,hbar,boltz,pi
      real(8) :: eps0,avsno
      real(8), dimension(3,3) :: cvec,cveci
      real(8) :: ewald_eps,kcut,rcut,rcut2
      real(8) :: eps_sqrtpii
      real(8), dimension(3,3) :: dudcvec,dudcvec_coulomb,dudcvec_pair,dudcvec_intra,dudcvec_pv
      logical :: rigid,periodic
      real(8), dimension(:,:), allocatable :: rfrac
      real(8) :: cmax
      integer, dimension(:), allocatable :: dimensionality
      integer :: n2bonds,n3bonds,n2bondtypes,n3bondtypes
      real(8), dimension(:,:), allocatable :: bond_splinecoeff,bond_xspline,bond_yspline
      real(8) :: ubond
      real(8) :: pressure,pressure_internal,density,cell_mass,upv
      integer, dimension(:), allocatable :: nmoloftype
      real(8) :: pefac
      character(len=10) :: input_energy_unit,output_energy_unit
      real(8) :: r_closest
      character(len=16) :: read_mode,read_cvec,write_cvec,rigid_output_mode,rigid_input_mode
      integer, dimension(:,:), allocatable :: pair_index
      real(8), dimension(:,:,:,:), allocatable :: fm_force
      real(8), dimension(:,:,:), allocatable :: fm_forcemol,fm_torquemol
      integer, dimension(:,:), allocatable :: model_list
      real(8), dimension(:,:), allocatable :: model_coeff
      integer, dimension(:,:), allocatable :: model_power
      integer, dimension(:,:), allocatable :: model_param
      integer, dimension(:), allocatable :: model_nterms
      integer :: nmodel_interaction
      logical :: calc_fm
      real(4) :: pad
      real(8), dimension(:,:), allocatable :: rhist
      real(8) :: rhmax,rhmin
      real(8), dimension(:), allocatable :: rsplinemax
      integer :: nhbins
      integer :: nspline
      logical :: use_spline
      real(8) :: damp_width,rdamp_cutoff,eps_damp,eps_damp_sqrtpii
      character(len = 8), dimension(:), allocatable :: bond_name
      integer(8) nbondlist,npair,nmodel_param
      character(len = 32) :: batch_file
      real(8) :: swap_chance,swap_conformer_chance
      real(8) :: mintol
      integer max_rank
      integer, dimension(:), allocatable :: mol_conformer,mol_nconformers
      real(8), dimension(:,:), allocatable :: mol_conformer_energy
      integer model_index_max
      real(8), dimension(:,:,:,:), allocatable :: trans_mat
      integer, dimension(:,:), allocatable :: site_rank

      integer nmoltype_max,nconformer_max,natom_max,nmol_max
      integer nspline_max,nbondlist_max,nmodel_param_max,nterms_max,nbin_max

      end module common_data

