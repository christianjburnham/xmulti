      module allocate
      contains
      
      subroutine allocate_arrays
      use common_data
      use parse_text
      implicit none
      character(len=2048) :: buffer
      character(len = 2048), dimension(64) :: textlist
      integer ios,nitems,line
      logical file_exists

!    defaults

      natom_max = 100
      nconformer_max = 64 
      nmoltype_max = 500 
      nmol_max = 100
      nspline_max = 256
      nbondlist_max = 500
      nmodel_param_max = 200
      nterms_max = 200
      nbin_max = 500

!     read in the size.dat file

      inquire(file = 'size.dat',exist = file_exists)
      if(.not.file_exists) call file_doesnt_exist('size.dat')
      open(45,file = 'size.dat')

      line = 0

      do while(.true.)
         read(45, '(A)', iostat=ios) buffer
         if(ios.ne.0) exit
         line = line + 1
         call split_text_to_list(buffer,textlist,nitems)
         select case(upcase(textlist(1)))
      case('NATOM_MAX')
         read(textlist(2),*,iostat = ios) natom_max
      case('NCONFORMER_MAX')
         read(textlist(2),*,iostat = ios) nconformer_max
      case('NMOLTYPE_MAX')
         read(textlist(2),*,iostat = ios) nmoltype_max
      case('NMOL_MAX')
         read(textlist(2),*,iostat = ios) nmol_max
      case('NSPLINE_MAX')
         read(textlist(2),*,iostat = ios) nspline_max
      case('NBONDLIST_MAX')
         read(textlist(2),*,iostat = ios) nbondlist_max
      case('NMODEL_PARAM_MAX')
         read(textlist(2),*,iostat = ios) nmodel_param_max
      case('NTERMS_MAX') 
         read(textlist(2),*,iostat = ios) nterms_max
      case('NBIN_MAX')
         read(textlist(2),*,iostat = ios) nbin_max
      case default
         print *, 'SIZE FILE: SKIPPING INVALID LABEL AT LINE', line               
      end select
      end do 
      
      allocate(site_coord0(3,natom_max,nconformer_max,nmoltype_max))
      allocate(lab_charge(natom_max,nmol_max))
      allocate(lab_dipole(3,natom_max,nmol_max))
      allocate(lab_quad(3,3,natom_max,nmol_max))
      allocate(lab_oct(3,3,3,natom_max,nmol_max))
      allocate(lab_coord(3,natom_max,nmol_max))
      allocate(lab_relative_coord(3,natom_max,nmol_max))
      allocate(site_name(nmoltype_max,natom_max))
      allocate(mol_natoms(nmoltype_max))
      allocate(mol_nmassatoms(nmoltype_max))
      allocate(mol_type(nmol_max))
      allocate(mol_euler(3,nmol_max))
      allocate(mol_com(3,nmol_max))
      allocate(force(3,natom_max,nmol_max))
      allocate(mass(natom_max,nmoltype_max))
      allocate(forcemol(3,nmol_max))
      allocate(fracforcemol(3,nmol_max))
      allocate(torquemol(3,nmol_max))
      allocate(dudeulermol(3,nmol_max))
      allocate(molecule_name_list(nmoltype_max))
      
      allocate(site_charge0(natom_max,nconformer_max,nmoltype_max))
      allocate(site_dipole0(3,natom_max,nconformer_max,nmoltype_max))
      allocate(site_quad0(3,3,natom_max,nconformer_max,nmoltype_max))
      allocate(site_oct0(3,3,3,natom_max,nconformer_max,nmoltype_max))

      allocate(cfield(natom_max,nmol_max))
      allocate(dfield(3,natom_max,nmol_max))
      allocate(qfield(3,3,natom_max,nmol_max))
      allocate(ofield(3,3,3,natom_max,nmol_max))

      allocate(lab_spherical2(5,natom_max,nmol_max))
      allocate(lab_spherical3(7,natom_max,nmol_max))
      allocate(lab_spherical2_1(5,3,natom_max,nmol_max))
      
      allocate(fieldtens0(natom_max,nmol_max))
      allocate(fieldtens1(3,natom_max,nmol_max))
      allocate(fieldtens2(5,natom_max,nmol_max))
      allocate(fieldtens3(7,natom_max,nmol_max))
      allocate(fieldtens1_1(3,3,natom_max,nmol_max))
      allocate(fieldtens2_1(5,3,natom_max,nmol_max))

      allocate(rfrac(3,nmol_max))
      allocate(dimensionality(nmoltype_max))

      allocate(bond_splinecoeff(nspline_max,nbondlist_max))
      allocate(bond_xspline(nspline_max,nbondlist_max))
      allocate(bond_yspline(nspline_max,nbondlist_max))

      allocate(nmoloftype(nmoltype_max))

      allocate(pair_index(natom_max,nmoltype_max))

      allocate(fm_force(nmodel_param_max,3,natom_max,nmol_max))
      allocate(fm_forcemol(nmodel_param_max,3,nmol_max))
      allocate(fm_torquemol(nmodel_param_max,3,nmol_max))

      allocate(model_list(nbondlist_max,3))
      allocate(model_coeff(nterms_max,nbondlist_max))
      allocate(model_power(nterms_max,nbondlist_max))
      allocate(model_param(nterms_max,nbondlist_max))
      allocate(model_nterms(nbondlist_max))
      allocate(rhist(nbondlist_max,0:nbin_max))
      allocate(rsplinemax(nbondlist_max))
      allocate(bond_name(nbondlist_max))

      allocate(mol_conformer(nmol_max))
      allocate(mol_nconformers(nmoltype_max))
      allocate(mol_conformer_energy(nconformer_max,nmoltype_max))

      allocate(trans_mat(0:3,nconformer_max,nconformer_max,nmoltype_max))
      allocate(site_rank(natom_max,nmoltype_max))
      end subroutine allocate_arrays
      end module allocate
