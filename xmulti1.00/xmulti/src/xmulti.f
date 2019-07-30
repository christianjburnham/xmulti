!*****************************************************************
!***   X-MULTI
!***  
!***   A code for finding crystal structures in triclinic periodic boundary conditions.
!***
!***   Christian J Burnham, University College Dublin, 2019
!***
!*****************************************************************

      call start()
      end 

      subroutine start()
      use common_data
      use allocate
      implicit none
      character(len = 10) :: mintype
      character(len = 8) :: antype
      real(8) :: basin_temp,press,radius,d_radius,ewald_error,schance,tolerance
      real(8) :: delta_frac,delta_xyz,delta_cvec,delta_angle,delta
      character(len = 32), dimension(1024) :: control_list
      real(8), dimension(1024,8) :: param_list
      real(8) :: rmin,rdamp
      integer ncontrol
      logical :: restart,file_exists
      integer i,n,imol
      real(8) :: dismin
      logical :: too_close
      character(len=6) :: version
      character(len=32) :: arg

      version = '1.03'
      
      do i = 1,command_argument_count()
         call get_command_argument(i, arg)
         select case(arg)
            case('-v')
               write(*,*) 'VERSION ',version
               stop
         end select
      end do 

      inquire(file = 'control.dat',exist = file_exists)
      if(.not.file_exists) call file_doesnt_exist('control.dat')

      open(40,file = 'control.dat')
      open(70,file = 'output.xyz')
      open(71,file = 'output.rigid')

      write(*,*)   "       -------------------------     "
      write(*,*)   "      /                       /      " 
      write(*,*)   "     /       X-MULTI         /       "
      write(*,*)   "    /                       /        "
      write(*,*)   "   -------------------------         "

      write(*,*) '*************************************&
     &*************************************************'
      write(*,*) '***   A code for finding crystal structures&
      & in triclinic periodic boundary conditions.'
      write(*,*) '***   Christian J Burnham, UCD, 2019'
      write(*,*) '*************************************&
     &*************************************************'

      call allocate_arrays

      call set_units()
!     default values

      mintol = 0.0d0 
      swap_chance = 0.0d0 
      use_spline = .true.
      rdamp_userset = .false.
      calc_fm = .false.
      pefac = fac(3)
      output_energy_unit = 'kJ/mol'
      write_cvec = 'XYZ'
      rigid_output_mode = 'FRAC'
      rdamp = 0.0d0
      call set_damp(rdamp)
      press = 0.0d0 
      call set_pressure(press)
      ewald_error = 1.d-7
      rcut = 8.0d0 
      radius = rcut
      call ewald_setup(rcut,ewald_error)
      r_closest = 0.6d0

      basin_temp = 250.0d0 
      delta_frac = 0.08d0
      delta_xyz = 1.0d0 
      delta_cvec = 0.05d0 
      delta_angle = 1.0d0 
      mol_conformer = 1

      call read_control(control_list,param_list,ncontrol)
      call read_molecular_geometry
      call get_cell_mass()

      if(periodic) then
         call generate_cartesians_from_frac()
      else
         call generate_cartesians()
      endif
      do n = 1,ncontrol
         select case(control_list(n))
         case('OUTPUT KJ')
            write(*,*)
            write(*,*) 'SETTING OUTPUT UNITS TO KJ/MOL'
            pefac = fac(3)
            output_energy_unit = 'kJ/mol'
         case('OUTPUT KCAL')
            write(*,*)
            write(*,*) 'SETTING OUTPUT UNITS TO KCAL/MOL'
            pefac = fac(4)
            output_energy_unit = 'kcal/mol'
         case('CVEC OUTPUT XYZ')
            write(*,*)
            write(*,*) 'SETTING OUTPUT CVEC MODE TO XYZ'
            write_cvec = 'XYZ'
         case('CVEC OUTPUT ANGLE')
            write(*,*)
            write(*,*) 'SETTING OUTPUT CVEC MODE TO ANGLE'
            write_cvec = 'ANGLE'            
         case('RIGID OUTPUT FRAC')
            write(*,*) 
            write(*,*) 'SETTING RIGID BODY OUTPUT TO FRACTIONAL + EULER'
            rigid_output_mode = 'FRAC'
         case('RIGID OUTPUT XYZ')
            write(*,*) 
            write(*,*) 'SETTING RIGID BODY OUTPUT TO XYZ + EULER'
            rigid_output_mode = 'XYZ'            
         case('PRESSURE')
            press = param_list(n,1)
            write(*,*) 
            write(*,23) press
 23         format('SETTING PRESSURE TO',f8.2,' bar')
            call set_pressure(press)            
         case('TOLERANCE')
            tolerance = param_list(n,1)
            write(*,*)
            write(*,29) tolerance,output_energy_unit
 29         format('SETTING LOCAL MINIMIZATION TOLERANCE TO',f14.6,a8)
            call set_tolerance(tolerance)
         case('ENERGY')
            write(*,*)
            write(*,*) 'CALLING ENERGY'
            call energy()
            write(*,*) 'ENERGIES IN ',output_energy_unit
            write(*,17) 'utotal = ',utotal * pefac
            write(*,17) 'upair = ',upair * pefac
            write(*,17) 'ucharge = ',ucharge * pefac
            write(*,17) 'udipole = ',udipole * pefac
            write(*,17) 'uquad = ',uquad * pefac
            write(*,17) 'uoct = ',uoct * pefac
            write(*,17) 'upv = ',upv * pefac
            write(*,17) 'uconformer = ',uconformer * pefac
 17        format(a12,f14.6,' ',a8)


         case('BATCH')
            write(*,*) 'BATCH READ'
            call batch()
         case('BATCH_XYZ')
            write(*,*) 'BATCH_XYZ READ'
            call batch_xyz()
         case('FORCES')
            write(*,*) 
            write(*,*) 'TESTING FORCES'
            delta = param_list(n,1)
            call forcetest(delta)     
         case('FFORCES')
            write(*,*) 
            write(*,*) 'TESTING FRAC FORCES'
            delta = param_list(n,1)
            call fracforcetest(delta)     
         case('RFORCES')
            write(*,*) 
            write(*,*) 'TESTING RIGID BODY FORCES'
            delta = param_list(n,1)
            call forcetest_rigid(delta)
         case('EULER')
            write(*,*) 
            write(*,*) 'TESTING EULER DERIVATIVES'
            delta = param_list(n,1)
            call eulertest(delta)
         case('CVEC')
            write(*,*) 
            write(*,*) 'TESTING CELL-VECTOR DERIVATIVES'
            delta = param_list(n,1)
            call test_cvec_derivatives(delta)     
         case('TORQUES')
            write(*,*) 
            write(*,*) 'TESTING TORQUES'
            delta = param_list(n,1)
            call test_torques(delta)     
         case('STOP')
            write(*,*)
            write(*,*) 'RECEIVED STOP INSTRUCTION IN CONTROL.'
            stop
         case('RCUT')
            radius = param_list(n,1)
            write(*,*)
            write(*,30) radius
 30         format('SETTING CUT-OFF RADIUS TO',f10.6,' ANGSTROMS')
            if(periodic) call ewald_setup(radius,ewald_error)
         case('RDAMP_CUT')
            d_radius = param_list(n,1)
            write(*,*)
            write(*,31) d_radius
 31         format('SETTING CHARGE-DAMPING &
     &CUT-OFF RADIUS TO',f10.6,' ANGSTROMS')
            call set_rdampcut(d_radius)
         case('BASIN RMIN')
            rmin = param_list(n,1)
            write(*,*)
            write(*,24) rmin
 24         format('SETTING BASIN RMIN TO',f8.4,' ANGSTROMS')
            call set_r_closest(rmin)
         case('PRINT_MULTIPOLES')
            if(periodic) then 
               call generate_cartesians_from_frac()
            else
               call generate_cartesians()
            endif
            write(*,*) 
            write(*,*) 'PRINTING MULTIPOLES'
            call print_multipoles
         case('PRINT_XYZ')
            if(periodic) then 
               call generate_cartesians_from_frac()
            else
               call generate_cartesians()
            endif
            call print_coordinates(70)
            write(*,*)
            write(*,*) 'PRINTING CARTESIAN COORDINATES TO OUTPUT.XYZ'
         case('PRINT_RIGID')
            call print_rigid_coordinates(71)
            write(*,*)
            write(*,*) 'PRINTING RIGID-BODY COORDINATES TO OUTPUT.RIGID'
         case('EWALD_ERROR')
            ewald_error = param_list(n,1)
            write(*,*)
            write(*,43) ewald_error
 43         format('SETTING EWALD ERROR TERM TO ',e14.6)
            if(periodic) call ewald_setup(radius,ewald_error)
         case('FINDMIN')
            write(*,*)
            write(*,*) 'MINIMIZING WITH TOLERANCE ',mintol*pefac
            call minimize()
         end select
      enddo

      end subroutine start

      
      subroutine ewald_setup(radius,ewald_error)
!     sets the parameters for the Ewald sum
      use common_data
      use nr_mod
      implicit none
      logical okflag
      real(8) :: ewald_error,radius

      call mat3inverse(cvec,cveci,okflag)

      rcut = radius
      rcut2 = rcut*rcut

!     calculate the ewald parameter from the ewald error term
      ewald_eps = dsqrt(-dlog(ewald_error)) / rcut

!     calculate the reciprocal-space cut-off. 
      kcut = 2.0d0 * ewald_eps**2 * rcut

!     reciprocal-space cut-off

      eps_sqrtpii = 1.0d0/(ewald_eps*dsqrt(pi))

      end subroutine ewald_setup
      

      subroutine copyvec(a,b)
      use common_data
      implicit none
      integer ivec,imol
      real(8), dimension(3,*) :: a,b
      
      do imol = 1,nmol
         do ivec = 1,3
            b(ivec,imol) = a(ivec,imol)
         end do 
      end do 
      
      end subroutine copyvec


      subroutine set_units()
      use common_data
      use nr_common_data
      implicit none

!     sets all the units used in the simulation

      pi = 4.0d0 * atan(1.0d0)
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
      fac(10) = 1.0d0 / 0.20819434
      
      boltz = (1.38064852d-23) * unit(5)/fac(1)
      hbar = 1.0545718d-34 /(unit(1)*fac(1))

      end subroutine set_units


      subroutine minimize()
!     finds the local minimum
      use common_data
      use nr_common_data
      use nr_mod
      implicit none
      integer i,j
      real(8) :: xvec(60000)
      real(8) :: ftol,fret
      integer :: imol,n,ieuler,ivec,jvec,iter,moltype,imode
      character(len = 10) :: mintype
      
      n = 0
      do imol = 1,nmol
         moltype = mol_type(imol)
         do ivec = 1,3
            n = n + 1
            if(periodic) then
               xvec(n) = rfrac(ivec,imol)
            else
               xvec(n) = mol_com(ivec,imol)
            endif
         end do 
         if(dimensionality(moltype).gt.0) then 
            do ieuler = 1,3
               n = n + 1
               xvec(n) = mol_euler(ieuler,imol)
            end do 
         endif
      end do 
      if(periodic) then
         do ivec = 1,3
            do jvec = 1,ivec
               n = n + 1
               xvec(n) = cvec(jvec,ivec)
            end do 
         end do 
      endif
      ftol = mintol
      linmin_param = 1.d-9
      nprintdata = 2000

      call frprmn(xvec,n,ftol,iter,fret)

      if(periodic) then 
           call generate_cartesians_from_frac()
        else
           call generate_cartesians()
        endif
      call energy()
      end subroutine minimize

      real(8) function func(xvec)
      use common_data
      implicit none
      real(8) :: xvec(60000)
      integer :: imol,n,ieuler,ivec,jvec,iter,moltype,imode

      func = 0.0d0 

      n = 0 
      do imol = 1,nmol
         moltype = mol_type(imol)
         do ivec = 1,3
            n = n + 1
            if(periodic) then
               rfrac(ivec,imol) = xvec(n)
            else
               mol_com(ivec,imol) = xvec(n)
            endif
         end do 
         if(dimensionality(moltype).gt.0) then
            do ieuler = 1,3
               n = n + 1
               mol_euler(ieuler,imol) = xvec(n)
            end do 
         endif
      end do 
      if(periodic) then
         do ivec = 1,3
            do jvec = 1,ivec
               n = n + 1
               if(ivec.ne.jvec) then 
                  cvec(jvec,ivec) = xvec(n)
               else
                  cvec(jvec,ivec) = dabs(xvec(n))
               endif
            end do 
         end do 
      endif

      if(periodic) then
         call generate_cartesians_from_frac()
      else
         call generate_cartesians()
      endif
      call energy()

      call calculate_rigid_body_forces()
      if(periodic) call get_rigid_frac_forces()

!      call print_coordinates(88)
!      write(*,*) 'uuu',utotal * pefac

      func = utotal * pefac

      end function func

      subroutine dfunc(xvec,gvec)
      use common_data
      implicit none
      real(8) :: xvec(60000),gvec(60000)
      integer :: imol,n,ieuler,ivec,jvec,iter,moltype,imode

      n = 0 
      do imol = 1,nmol
         moltype = mol_type(imol)
         do ivec = 1,3
            n = n + 1
            if(periodic) then
               gvec(n) = -fracforcemol(ivec,imol) * pefac
            else
               gvec(n) = -forcemol(ivec,imol) * pefac
            endif
         end do 
         if(dimensionality(moltype).gt.0) then 
            do ieuler = 1,3
               n = n + 1
               gvec(n) = -dudeulermol(ieuler,imol) * pefac
            end do 
         endif
      end do 
      if(periodic) then
         do ivec = 1,3
            do jvec = 1,ivec
               n = n + 1
               gvec(n) = dudcvec(jvec,ivec) * pefac
            end do 
         end do 
      endif

      end subroutine dfunc

      subroutine TEST_CVEC_DERIVATIVES(delta)
!     tests the derivatives with respect to the cell vectors.
      use common_data
      implicit none
      real(8) :: utot0,utot1
      integer:: ivec,jvec
      real(8) :: delta,force_numeric
      real(8), dimension(3,3) :: dudcvec0
      integer natoms_in_mol,iatom,imol,imoltype
      character(len = 2), dimension(3,3) :: cvec_label

      cvec_label(:,1) = (/'ax','ay','az'/)
      cvec_label(:,2) = (/'bx','by','bz'/)
      cvec_label(:,3) = (/'cx','cy','cz'/)

      if(.not.periodic) then 
         write(*,*) 'ERROR: CANNOT TEST CELL VECTOR DERIVATIVES. &
     &NO PERIODIC BOUNDARY CONDITIONS.'
         return
      endif

      write(*,*) 'derivatives with respect to cell vector components'
      call get_frac_coords()
      call generate_cartesians_from_frac()
      call energy()

      utot0 = utotal
      dudcvec0 = dudcvec

      do ivec = 1,3
         do jvec = ivec,3
            cvec(ivec,jvec) = cvec(ivec,jvec) + delta
            
            call generate_cartesians_from_frac()
            call energy()
            
            utot1 = utotal
            force_numeric = (utot0-utot1)/delta

            write(*,17) cvec_label(ivec,jvec),force_numeric * pefac,-dudcvec0(ivec,jvec)*pefac,&
     &-dudcvec0(ivec,jvec)*delta/(utot0-utot1)
 
            cvec(ivec,jvec) = cvec(ivec,jvec) - delta
         end do 
      end do 
      
 17   format(a2,' numerical ',f14.6,' analytical ',f14.6,' ratio ',f14.6)
      end subroutine test_cvec_derivatives

      subroutine test_torques(delta)
      use common_data
      implicit none
      integer imol,moltype
      real(8) :: phi0,theta0,psi0,phi,theta,psi,stheta
      real(8), dimension(3,3) :: rotmat,rotmat2,axismat
      real(8), dimension(3) :: uvec
      real(8) :: delta,sum,utot0,s,g
      real(8), dimension(:,:), allocatable :: torquemol0
      integer ivec,jvec,kvec,lvec,mvec
      character(len = 1),dimension(3) :: direction

      direction = (/'x','y','z'/)

      allocate(torquemol0(3,nmol_max))
      
      if(periodic) then
         call generate_cartesians_from_frac()
      else
         call generate_cartesians()
      endif
      call energy
      call calculate_rigid_body_forces()
      utot0 = utotal 
      torquemol0 = torquemol

      do imol = 1,nmol
         moltype = mol_type(imol)
         
         phi0 = mol_euler(1,imol)
         theta0 = mol_euler(2,imol)
         psi0 = mol_euler(3,imol)
         
         call get_rotmat(phi0,theta0,psi0,rotmat)

!     loop over rotation axes, x,y,z
         do ivec = 1,3

!     uvec is the rotation axes about which we're going to calculate the torque
            uvec = 0.0d0 
            uvec(ivec) = 1.0d0 

!     calculate the rotation matrix for a rotation about this axis
            call get_rotmat_for_axis(uvec,delta,axismat)

!     now multiply matrices to generate a combined rotation matrix
            
            do jvec = 1,3
               do kvec = 1,3
                  sum = 0.0d0 
                  do lvec = 1,3
                     sum = sum + axismat(jvec,lvec) * rotmat(lvec,kvec) 
                  end do 
                  rotmat2(jvec,kvec) = sum 
               end do 
            end do 

!     find the new values of phi,theta,psi

            call get_euler_from_rotmat(rotmat2,phi,theta,psi)
            
            mol_euler(1,imol) = phi
            mol_euler(2,imol) = theta
            mol_euler(3,imol) = psi 

            if(periodic) then 
               call generate_cartesians_from_frac()
            else 
               call generate_cartesians()
            endif
            call energy

            write(*,17) imol,direction(ivec),-pefac * (utotal - utot0) / delta,pefac * torquemol0(ivec,imol),& 
     &           -(utotal - utot0) / (delta * torquemol0(ivec,imol))

 17            format('mol ',i5,a3,' numerical ',f14.6,&
     &'    analytical ',f14.6,' ratio ',f14.6)
            
            mol_euler(1,imol) = phi0
            mol_euler(2,imol) = theta0 
            mol_euler(3,imol) = psi0 
         end do 
         
      end do 
         
      end subroutine test_torques

      subroutine generate_cartesians_from_frac()
!     calculates the cartesian coordinates from the fractional coordinates
      use common_data
      use nr_mod
      implicit none
      integer natoms_in_mol,imol,imoltype,iatom
      logical okflag
      real(8), dimension(3) :: vec1,vec2

      call mat3inverse(cvec,cveci,okflag)

      do imol = 1,nmol
         vec1(1) = rfrac(1,imol) - anint(rfrac(1,imol))
         vec1(2) = rfrac(2,imol) - anint(rfrac(2,imol))
         vec1(3) = rfrac(3,imol) - anint(rfrac(3,imol))
         
         call mat3multvec(cvec,vec1,vec2)
         
         mol_com(1,imol) = vec2(1)
         mol_com(2,imol) = vec2(2) 
         mol_com(3,imol) = vec2(3)
         
      end do 

      call generate_cartesians()

      end subroutine generate_cartesians_from_frac


      subroutine get_frac_coords()
!     calculates the fractional coordinates
      use common_data
      use nr_mod
      implicit none
      integer imol,imoltype,natoms_in_mol,iatom
      logical okflag
      real(8), dimension(3) :: vec1,vec2
      
      call mat3inverse(cvec,cveci,okflag)

      do imol = 1,nmol
         vec1(1) = mol_com(1,imol)
         vec1(2) = mol_com(2,imol)
         vec1(3) = mol_com(3,imol)
         
         call mat3multvec(cveci,vec1,vec2)
         
         rfrac(1,imol) = vec2(1)
         rfrac(2,imol) = vec2(2)
         rfrac(3,imol) = vec2(3)
      end do 

      end subroutine get_frac_coords

      subroutine forcetest(delta)
!     tests the XYZ cartesian forces
      use common_data
      implicit none
      integer moltype,natoms_in_mol,i,j,iatom
      real(8) :: utot0,utot1
      integer ivec,imol
      real(8) :: delta,force_numeric
      real(8), dimension(:,:,:), allocatable :: force0
      character(len = 1),dimension(3) :: direction
      
      direction = (/'x','y','z'/)

      allocate(force0(3,natom_max,nmol_max))

      call energy()
      utot0 = utotal

      do imol = 1,nmol
         moltype = mol_type(imol)
         natoms_in_mol = mol_natoms(moltype)
         do iatom = 1,natoms_in_mol
            do ivec = 1,3
            force0(ivec,iatom,imol) = force(ivec,iatom,imol)
            end do 
         end do 
         
      end do 

      do imol = 1,nmol
         moltype = mol_type(imol)
         natoms_in_mol = mol_natoms(moltype)

         do iatom = 1,natoms_in_mol
            do ivec = 1,3
               lab_coord(ivec,iatom,imol) = lab_coord(ivec,iatom,imol) + delta

               call energy()
               utot1 = utotal

               force_numeric = (utot0 - utot1) / delta

               write(*,17) imol,iatom,direction(ivec),force_numeric * pefac,force0(ivec,iatom,imol) * pefac,&
     & force_numeric / force0(ivec,iatom,imol)
 17            format('mol ',i5,' atom ',i5,a3,' numerical ',f14.6,&
     &'    analytical ',f14.6,' ratio ',f14.6)

               lab_coord(ivec,iatom,imol) = lab_coord(ivec,iatom,imol) - delta

            end do 
         end do 
      end do 

      end subroutine forcetest

      subroutine eulertest(delta)
!     tests the euler angle derivatives
      use common_data
      implicit none
      integer imol,natoms_in_mol,moltype,ieuler
      real(8) :: utot0,utot1,delta,force_numeric
      real(8), dimension(:,:), allocatable :: dudeulermol0

      allocate(dudeulermol0(3,nmol_max))

      call generate_cartesians()
      call energy()
      call calculate_rigid_body_forces()

      utot0 = utotal

      do imol = 1,nmol
         moltype = mol_type(imol)
         natoms_in_mol = mol_natoms(moltype)

         dudeulermol0(:,imol) = dudeulermol(:,imol)
      end do 

      write(*,*) 'Derivatives with respect to Euler angles'
      do imol = 1,nmol

!     loop over Euler angles
         do ieuler = 1,3
            mol_euler(ieuler,imol) = mol_euler(ieuler,imol) + delta

            call generate_cartesians()
            call energy()

            utot1 = utotal
            force_numeric = (utot0 - utot1) / delta

            write(*,17) imol,ieuler,force_numeric * pefac,dudeulermol0(ieuler,imol) * pefac,&
     & force_numeric/dudeulermol0(ieuler,imol)
 17            format('mol ',i5,' angle ',i3,' numerical ',f14.6,&
     &'    analytical ',f14.6,' ratio ',f14.6)

            mol_euler(ieuler,imol) = mol_euler(ieuler,imol) - delta
         enddo 

      end do 

      end subroutine eulertest

      subroutine fracforcetest(delta)
!     tests the fractional derivatives
      use common_data
      implicit none
      integer moltype,natoms_in_mol,i,j,iatom
      real(8) :: utot0,utot1
      integer ivec,imol,ieuler
      real(8) :: delta,force_numeric
      real(8), dimension(:,:),allocatable :: forcemol0,dudeulermol0
      character(len = 1),dimension(3) :: direction
      
      direction = (/'1','2','3'/)

      allocate(forcemol0(3,nmol_max))
      allocate(dudeulermol0(3,nmol_max))

      call get_frac_coords()
      call energy()
      call calculate_rigid_body_forces()
      call get_rigid_frac_forces()

      utot0 = utotal

      do imol = 1,nmol
         do ivec = 1,3
            forcemol0(ivec,imol) = fracforcemol(ivec,imol)
         end do 
         do ieuler = 1,3
            dudeulermol0(ieuler,imol) = dudeulermol(ieuler,imol)
         end do 
      end do 

      write(*,*) 'DERIVATIVES WRT FRAC COM'

      do imol = 1,nmol
         do ivec = 1,3
            rfrac(ivec,imol) = rfrac(ivec,imol) + delta
            
            call generate_cartesians_from_frac()
            call energy()
            call calculate_rigid_body_forces()
            call get_rigid_frac_forces()

            utot1 = utotal
            force_numeric = (utot0 - utot1) / delta
            
            write(*,17) imol,direction(ivec),force_numeric * pefac,forcemol0(ivec,imol) * pefac,&
     & force_numeric / forcemol0(ivec,imol)

            rfrac(ivec,imol) = rfrac(ivec,imol) - delta
!     need one last call to make sure coordinates are reset to initial
            call generate_cartesians_from_frac()
            
         end do 
      end do 

 17         format('mol ',i5,a3,' numerical ',f14.6,&
     &'    analytical ',f14.6,' ratio ',f14.6)

      end subroutine fracforcetest

      subroutine forcetest_rigid(delta)
!     tests the rigid-body forces
      use common_data
      implicit none
      integer moltype,natoms_in_mol,i,j,iatom
      real(8) :: utot0,utot1
      integer ivec,imol,ieuler
      real(8) :: delta,force_numeric
      character(len = 1),dimension(3) :: direction
      real(8), dimension(:,:), allocatable :: forcemol0, dudeulermol0

      direction = (/'x','y','z'/)

      allocate(forcemol0(3,nmol_max))
      allocate(dudeulermol0(3,nmol_max))

      call generate_cartesians()
      call energy()
      call calculate_rigid_body_forces()

      utot0 = utotal

      do imol = 1,nmol
         moltype = mol_type(imol)
         natoms_in_mol = mol_natoms(moltype)

         forcemol0(:,imol) = forcemol(:,imol)
         dudeulermol0(:,imol) = dudeulermol(:,imol)
      end do 

      write(*,*) 
      write(*,*) 'Derivatives with respect to COM coords'
      do imol = 1,nmol
         moltype = mol_type(imol)
         natoms_in_mol = mol_natoms(moltype)

         do ivec = 1,3
            mol_com(ivec,imol) = mol_com(ivec,imol) + delta
            
            call generate_cartesians()
            call energy()

            utot1 = utotal
            force_numeric = (utot0 - utot1) / delta

            write(*,17) imol,direction(ivec),force_numeric * pefac,forcemol0(ivec,imol) * pefac,&
     &  force_numeric / forcemol0(ivec,imol)

            mol_com(ivec,imol) = mol_com(ivec,imol) - delta
         end do 
      end do 

 17         format('mol ',i5,a3,' numerical ',f14.6,&
     &'    analytical ',f14.6,' ratio ',f14.6)

      end subroutine forcetest_rigid


      subroutine initialize()
!     initializes the energies / forces
      use common_data
      implicit none
      integer imol,imoltype,iatom,imode
      integer m,n,j,ivec,jvec,kvec,ii,jj

      utotal = 0.0d0 

      do imol = 1,nmol
         imoltype = mol_type(imol)
      end do 

      do imol = 1,nmol
         forcemol(1,imol) = 0.0d0 
         forcemol(2,imol) = 0.0d0 
         forcemol(3,imol) = 0.0d0 

         torquemol(1,imol) = 0.0d0 
         torquemol(2,imol) = 0.0d0 
         torquemol(3,imol) = 0.0d0 

         dudeulermol(1,imol) = 0.0d0 
         dudeulermol(2,imol) = 0.0d0 
         dudeulermol(3,imol) = 0.0d0 

         imoltype = mol_type(imol)
         do iatom = 1,mol_natoms(imoltype)
            cfield(iatom,imol) = 0.0d0 
            if(max_rank.ge.1) then
               do ivec = 1,3
                  dfield(ivec,iatom,imol) = 0.0d0 
                  if(max_rank.ge.2) then 
                     do jvec = 1,3
                        qfield(ivec,jvec,iatom,imol) = 0.0d0 
                        if(max_rank.ge.3) then
                           do kvec = 1,3
                              ofield(ivec,jvec,kvec,iatom,imol) = 0.0d0 
                           end do 
                        endif
                     end do 
                  endif
               end do 
            endif

            fieldtens0(iatom,imol) = 0.0d0 
            if(max_rank.ge.1) then 
               do ii = 1,3
                  fieldtens1(ii,iatom,imol) = 0.0d0 
               end do 
            endif
            if(max_rank.ge.2) then 
               do ii = 1,5
                  fieldtens2(ii,iatom,imol) = 0.0d0
               end do 
               do ii = 1,3
                  do jj = 1,3
                     fieldtens1_1(ii,jj,iatom,imol) = 0.0d0 
                  end do 
               end do 
            endif
            if(max_rank.ge.3) then 
               do ii = 1,7
                  fieldtens3(ii,iatom,imol) = 0.0d0 
               end do 
               do ii = 1,5
                  do jj = 1,3
                     fieldtens2_1(ii,jj,iatom,imol) = 0.0d0 
                  end do 
               end do 
            endif

            force(1,iatom,imol) = 0.0d0
            force(2,iatom,imol) = 0.0d0 
            force(3,iatom,imol) = 0.0d0 

         end do  
      end do 

      end subroutine initialize


      subroutine RECIPROCAL()
!----------------------------------------------
!     CALCULATES THE RECIPROCAL SPACE ENERGY + FORCES
!----------------------------------------------
      use common_data
      use nr_common_data
      use nr_mod
      implicit none

      integer :: ik1,ik2,ik3
      integer :: i,j,k,ivec,jvec,kvec,lvec,mvec,iatom,imoltype
      real(8) :: a,b,c,d,e,f,g,h,o,p,q,r,s,t,u,v,w,x,y,z
      real(8) :: sinkri,coskri,sinkrj,coskrj
      real(8) :: kcut2,kk,kx,ky,kz
      real(8) :: k1x,k1y,k1z,k2x,k2y,k2z,k3x,k3y,k3z
      integer :: k1max,k2max,k3max
      real(8) :: cos_charge_sum,sin_charge_sum,k_d_r,sum2,kk2,efac
      real(8) :: vol0,voli,prefac,multfac
      real(8) :: coskrij_charge_sum,sinkrij_charge_sum
      real(8) :: cos_dipi_d_k_sum,sin_dipi_d_k_sum
      real(8) :: coskrij,sinkrij
      real(8) :: chargei
      real(8) :: du
      real(8), dimension(:,:), allocatable :: cosk_d_r,sink_d_r
      real(8) :: vol_rec
      real(8) :: xdif,ydif,zdif
      real(8) :: ff,gg,fprefac,dfprefac
      real(8) :: sum_real,sum_imag,sum_ewald
      real(8) :: ax,ay,az,bx,by,bz,cx,cy,cz
      integer :: iox,ih1,ih2,im,imol,jmol,natoms_in_mol
      integer :: ii
      real(8), dimension(3) :: vec1,vec2
      real(8) :: rcx,rcy,rcz
      real(8) :: k_d_rc
      real(8) :: twopi,twopii,twopivoli
      real(8) :: a_rec_x,a_rec_y,a_rec_z,b_rec_x,b_rec_y,b_rec_z,c_rec_x,c_rec_y,c_rec_z
      real(8) :: axbmag_rec,bxcmag_rec,cxamag_rec
      real(8) :: axbx_rec,axby_rec,axbz_rec,bxcx_rec,bxcy_rec,bxcz_rec,cxax_rec,cxay_rec,cxaz_rec
      real(8) :: amag,bmag,cmag
      real(8) :: dipi_d_k
      real(8), dimension(3,3) :: cos_ddip_d_k_dcvec_sum,sin_ddip_d_k_dcvec_sum
      real(8) :: sinkrij_dip_sum,coskrij_dip_sum
      real(8) :: sinkrij_quad_sum,coskrij_quad_sum
      real(8) :: sinkrij_oct_sum,coskrij_oct_sum
      real(8), dimension(:,:), allocatable :: dip_d_k,quad_dd_kk,oct_ddd_kkk
      real(8) :: quadi_dd_kk,octi_ddd_kkk
      real(8) :: cos_quadi_dd_kk_sum,sin_quadi_dd_kk_sum
      real(8) :: cos_octi_ddd_kkk_sum,sin_octi_ddd_kkk_sum
      real(8), dimension(3,3) :: quad
      real(8), dimension(3,3,3) :: oct
      real(8), dimension(3) :: k_vec
      real(8), dimension(3,3,3) :: dkdcvec
      real(8) :: ddip_d_k_dcvec,dquad_dd_kk_dcvec,doct_ddd_kkk_dcvec
      real(8), dimension(3,3) :: cos_dquad_dd_kk_dcvec_sum,sin_dquad_dd_kk_dcvec_sum
      real(8), dimension(3,3) :: cos_doct_ddd_kkk_dcvec_sum,sin_doct_ddd_kkk_dcvec_sum
      real(8), dimension(3) :: dip
      real(8) :: kdkdcvec,dk_d_rdcvec
      real(8) :: xx,xy,xz,yy,yz,zz
      real(8) :: xxx,xxy,xyy,yyy,xxz,xyz,yyz,xzz,yzz,zzz
      real(8), dimension(5) :: ksolid2
      real(8), dimension(7) :: ksolid3
      real(8), dimension(5) :: spherical2
      real(8), dimension(7) :: spherical3

      allocate(cosk_d_r(natom_max,nmol_max))
      allocate(sink_d_r(natom_max,nmol_max))
      allocate(dip_d_k(natom_max,nmol_max))
      allocate(quad_dd_kk(natom_max,nmol_max))
      allocate(oct_ddd_kkk(natom_max,nmol_max))

      kcut2 = kcut*kcut

      ax = cvec(1,1) 
      ay = cvec(2,1) 
      az = cvec(3,1)

      bx = cvec(1,2)
      by = cvec(2,2)
      bz = cvec(3,2)

      cx = cvec(1,3)
      cy = cvec(2,3)
      cz = cvec(3,3)

      vol0 = ax * by * cz
      voli = 1.d0 / vol0

      ff = 8.0d0 * pi * voli 
      gg = 1.d0/ (4.0d0 * ewald_eps*ewald_eps)

!      write(*,*) 'gg',dexp(-kcut2 *gg) 
!      stop

      twopi = 2.0d0 * pi
      twopii = 1.0d0 / twopi
      twopivoli = 2.0d0 * pi * voli

      a_rec_x = 2.0d0 * pi * cveci(1,1)
      a_rec_y = 2.0d0 * pi * cveci(1,2) 
      a_rec_z = 2.0d0 * pi * cveci(1,3)

      b_rec_x = 2.0d0 * pi * cveci(2,1)
      b_rec_y = 2.0d0 * pi * cveci(2,2)
      b_rec_z = 2.0d0 * pi * cveci(2,3)

      c_rec_x = 2.0d0 * pi * cveci(3,1)
      c_rec_y = 2.0d0 * pi * cveci(3,2)
      c_rec_z = 2.0d0 * pi * cveci(3,3)

!     get the reciprocal volume a_rec.(b_rec X c_rec)
      vol_rec = (2.0d0 * pi)**3 / vol0

!     finds the size of the supercell which contains the cut-off sphere

      amag = dsqrt(cvec(1,1)**2 + cvec(2,1)**2 + cvec(3,1)**2)
      bmag = dsqrt(cvec(1,2)**2 + cvec(2,2)**2 + cvec(3,2)**2)
      cmag = dsqrt(cvec(1,3)**2 + cvec(2,3)**2 + cvec(3,3)**2)
      
      k1max = int(kcut * amag / (2.0d0 * pi) ) 
      k2max = int(kcut * bmag / (2.0d0 * pi) ) 
      k3max = int(kcut * cmag / (2.0d0 * pi) ) 

      sum_ewald = 0.0d0 

!     loop over k-vectors

      do ik1 = 0,k1max
         k1x = dble(ik1) * a_rec_x
         k1y = dble(ik1) * a_rec_y
         k1z = dble(ik1) * a_rec_z
         do ik2 = -k2max,k2max
            k2y = dble(ik2) * b_rec_y
            k2z = dble(ik2) * b_rec_z
            do ik3 = -k3max,k3max
!     only sum over half the lattice
               if((ik1.gt.0).or.(((ik2.gt.0.or.ik2.eq.0.and.ik3.gt.0)))) then 
               k3z = dble(ik3) * c_rec_z

!     works for upper triangular
               kx = k1x 
               ky = k1y + k2y 
               kz = k1z + k2z + k3z

               kk2 = kx**2 + ky**2 + kz**2
               if(kk2 < kcut2) then 
                  
                  k_vec(1) = kx
                  k_vec(2) = ky
                  k_vec(3) = kz
                  
!     derivatives of the k vector wrt the cell-vector components
                  dkdcvec(1,1,1) = -kx * voli * by * cz
                  
                  dkdcvec(1,1,2) = 0.0d0 
                  dkdcvec(1,2,2) = voli * (-kx * ax * cz + twopi * (ik1 * cz - ik3 * az))
                  
                  dkdcvec(1,1,3) = 0.0d0 
                  dkdcvec(1,2,3) = twopivoli * (ik2 * az - ik1 * bz)
                  dkdcvec(1,3,3) = voli * (-kx * ax * by + twopi * (ik1 * by - ik2 * ay))
                  
                  dkdcvec(2,1,1) = voli * (-ky * by * cz + twopi * (ik2 * cz - ik3 * bz))
                  
                  dkdcvec(2,1,2) = twopivoli * (ik3 * az - ik1 * cz)
                  dkdcvec(2,2,2) = -ky * voli * ax * cz
                  
                  dkdcvec(2,1,3) = twopivoli * (ik1 * bz - ik2 * az)
                  dkdcvec(2,2,3) = 0.0d0 
                  dkdcvec(2,3,3) = voli * (-ky * ax * by + twopi * (ik2 * ax - ik1 * bx))
                  
                  dkdcvec(3,1,1) = voli * (-kz * by * cz + twopi * (ik3 * by - ik2 * cy))
                  
                  dkdcvec(3,1,2) = twopivoli * (ik1 * cy - ik3 * ay)
                  dkdcvec(3,2,2) = voli * (-kz * ax * cz + twopi * (ik3 * ax - ik1 * cx))
                  
                  dkdcvec(3,1,3) = twopivoli * (ik2 * ay - ik1 * by)
                  dkdcvec(3,2,3) = twopivoli * (ik1 * bx - ik2 * ax)
                  dkdcvec(3,3,3) = -kz * voli * ax * by
                  

!     fprefac is the exp(-k^2/4 eps^2)/(V eps0 k^2) term
                  
                     fprefac = ff * dexp(-kk2 *gg) / kk2
                  
!     dfprefac is 1/k * the derivative of fprefac wrt k
                  
                     dfprefac = - fprefac * (4.0d0 * ewald_eps**2 + kk2) / (2.0d0 * ewald_eps**2 * kk2) 
                     
                     cos_charge_sum = 0.0d0
                     sin_charge_sum = 0.0d0

                     cos_dipi_d_k_sum = 0.0d0
                     sin_dipi_d_k_sum = 0.0d0
                     
                     cos_dquad_dd_kk_dcvec_sum = 0.0d0 
                     sin_dquad_dd_kk_dcvec_sum = 0.0d0 

                     cos_doct_ddd_kkk_dcvec_sum = 0.0d0 
                     sin_doct_ddd_kkk_dcvec_sum = 0.0d0 

                     cos_ddip_d_k_dcvec_sum = 0.0d0 
                     sin_ddip_d_k_dcvec_sum = 0.0d0 

                     cos_quadi_dd_kk_sum = 0.0d0 
                     sin_quadi_dd_kk_sum = 0.0d0 
                     
                     cos_octi_ddd_kkk_sum = 0.0d0 
                     sin_octi_ddd_kkk_sum = 0.0d0 

!     pre-calculate sums for later use
                     do imol = 1,nmol
                        imoltype = mol_type(imol)
                        natoms_in_mol = mol_natoms(imoltype)
                        do iatom = 1,natoms_in_mol
                           k_d_r = kx*lab_coord(1,iatom,imol) + ky*lab_coord(2,iatom,imol) + kz*lab_coord(3,iatom,imol)
                           
                           call sincos(k_d_r*pifac,s,c)
                           
                           cosk_d_r(iatom,imol) = c
                           sink_d_r(iatom,imol) = s
                           
                           chargei = lab_charge(iatom,imol)
                           
                           cos_charge_sum = cos_charge_sum + chargei * c
                           sin_charge_sum = sin_charge_sum + chargei * s
                           
                           if(max_rank.ge.1) then
                              
                              dip(1) = lab_dipole(1,iatom,imol)
                              dip(2) = lab_dipole(2,iatom,imol)
                              dip(3) = lab_dipole(3,iatom,imol)
                              
                              dipi_d_k = dip(1) * kx + dip(2) * ky + dip(3) * kz
                              
                              cos_dipi_d_k_sum = cos_dipi_d_k_sum + dipi_d_k * c
                              sin_dipi_d_k_sum = sin_dipi_d_k_sum + dipi_d_k * s
                              
                              dip_d_k(iatom,imol) = dipi_d_k
                              
!     calculate the derivatives of kx,ky,kz wrt the cell-vectors
                              
                              do jvec = 1,3
                                 do ivec = 1,jvec
                                    ddip_d_k_dcvec = dip(1) * dkdcvec(1,ivec,jvec) &
     &                                             + dip(2) * dkdcvec(2,ivec,jvec) &
     &                                             + dip(3) * dkdcvec(3,ivec,jvec)
                                    cos_ddip_d_k_dcvec_sum(ivec,jvec) = cos_ddip_d_k_dcvec_sum(ivec,jvec) + ddip_d_k_dcvec * c
                                    sin_ddip_d_k_dcvec_sum(ivec,jvec) = sin_ddip_d_k_dcvec_sum(ivec,jvec) + ddip_d_k_dcvec * s
                                 end do 
                              end do 
                           endif
                              
                           if(max_rank.ge.2) then 
!     quadrupoles
                              do ivec = 1,3
                                 do jvec = 1,3
                                    quad(ivec,jvec) = lab_quad(ivec,jvec,iatom,imol)
                                 end do 
                              end do 
                           
                              xx = kx * kx
                              xy = kx * ky
                              xz = kx * kz
                              yy = ky * ky
                              yz = ky * kz
                              zz = kz * kz
                              call convert_quad_to_spherical(xx,xy,xz,yy,yz,zz,ksolid2)

                              do ii = 1,5
                                 spherical2(ii) = lab_spherical2(ii,iatom,imol)
                              end do 

                              quadi_dd_kk = 0.0d0 
                              do ii = 1,5
                                 quadi_dd_kk = quadi_dd_kk + spherical2(ii) * ksolid2(ii)
                              end do 
                              
                              cos_quadi_dd_kk_sum = cos_quadi_dd_kk_sum + quadi_dd_kk * c
                              sin_quadi_dd_kk_sum = sin_quadi_dd_kk_sum + quadi_dd_kk * s
                              
                              quad_dd_kk(iatom,imol) = quadi_dd_kk
                              
                              do jvec = 1,3
                                 do ivec = 1,jvec
                                    
                                    dquad_dd_kk_dcvec = 0.0 
                                    do kvec = 1,3
                                       do lvec = 1,3
                                          dquad_dd_kk_dcvec = dquad_dd_kk_dcvec & 
     &                                  + 2.0d0 * quad(kvec,lvec) * k_vec(kvec) * dkdcvec(lvec,ivec,jvec)
                                       end do 
                                    end do 

                                    cos_dquad_dd_kk_dcvec_sum(ivec,jvec) = cos_dquad_dd_kk_dcvec_sum(ivec,jvec) + & 
     &                                   dquad_dd_kk_dcvec * c
                                    sin_dquad_dd_kk_dcvec_sum(ivec,jvec) = sin_dquad_dd_kk_dcvec_sum(ivec,jvec) + & 
     &                                   dquad_dd_kk_dcvec * s
                                 end do 
                              end do 
                           endif
                           if(max_rank.ge.3) then 
!     octopoles
                              do ivec = 1,3
                                 do jvec = 1,3
                                    do kvec = 1,3
                                       oct(ivec,jvec,kvec) = lab_oct(ivec,jvec,kvec,iatom,imol)
                                    end do 
                                 end do 
                              end do 

                              xxx = xx * kx
                              xxy = xx * ky
                              xyy = yy * kx
                              yyy = yy * ky
                              xxz = xx * kz
                              xyz = xy * kz
                              yyz = yy * kz
                              xzz = zz * kx
                              yzz = zz * ky
                              zzz = zz * kz
                              
                              call  convert_oct_to_spherical(xxx,xxy,xyy,yyy,xxz,xyz,yyz,xzz,yzz,zzz,ksolid3)

                              do ii = 1,7
                                 spherical3(ii) = lab_spherical3(ii,iatom,imol)
                              end do 
                              
                              octi_ddd_kkk = 0.0d0 
                              do ii = 1,7
                                 octi_ddd_kkk = octi_ddd_kkk + spherical3(ii) * ksolid3(ii)
                              end do 

                              cos_octi_ddd_kkk_sum = cos_octi_ddd_kkk_sum + octi_ddd_kkk * c
                              sin_octi_ddd_kkk_sum = sin_octi_ddd_kkk_sum + octi_ddd_kkk * s

                              oct_ddd_kkk(iatom,imol) = octi_ddd_kkk


                              do jvec = 1,3
                                 do ivec = 1,jvec
                                    
                                    doct_ddd_kkk_dcvec = 0.0 
                                    do kvec = 1,3
                                       do lvec = 1,3
                                          do mvec = 1,3
                                          doct_ddd_kkk_dcvec = doct_ddd_kkk_dcvec & 
     &          + 3.0d0 * oct(kvec,lvec,mvec) * k_vec(kvec) * k_vec(lvec) * dkdcvec(mvec,ivec,jvec)
                                          end do 
                                       end do 
                                    end do 
                                    
                                    cos_doct_ddd_kkk_dcvec_sum(ivec,jvec) = cos_doct_ddd_kkk_dcvec_sum(ivec,jvec) + & 
     &                                   doct_ddd_kkk_dcvec * c
                                    sin_doct_ddd_kkk_dcvec_sum(ivec,jvec) = sin_doct_ddd_kkk_dcvec_sum(ivec,jvec) + & 
     &                                   doct_ddd_kkk_dcvec * s
                                 end do 
                              end do 

                           endif

                        end do 
                     end do 
                     
!     main particle loop

                  iloop: do imol = 1,nmol
                  imoltype = mol_type(imol)
                  natoms_in_mol = mol_natoms(imoltype)
                  do iatom = 1,natoms_in_mol
                     
                     chargei = lab_charge(iatom,imol)
                     coskri = cosk_d_r(iatom,imol)
                     sinkri = sink_d_r(iatom,imol)
                     
                     sinkrij_charge_sum = sinkri*cos_charge_sum - coskri * sin_charge_sum
                     coskrij_charge_sum = coskri*cos_charge_sum + sinkri * sin_charge_sum

                     f = fprefac * sinkrij_charge_sum * chargei
                     fieldtens0(iatom,imol) = fieldtens0(iatom,imol) + fprefac * coskrij_charge_sum

                     dipi_d_k = 0.0d0 
                     if(max_rank.ge.1) then 
                        sinkrij_dip_sum = sinkri*cos_dipi_d_k_sum - coskri * sin_dipi_d_k_sum
                        coskrij_dip_sum = coskri*cos_dipi_d_k_sum + sinkri * sin_dipi_d_k_sum
                        dipi_d_k = dip_d_k(iatom,imol)

                        f = f + fprefac * dipi_d_k * sinkrij_dip_sum
                        f = f - fprefac * (chargei * coskrij_dip_sum - dipi_d_k * coskrij_charge_sum)

                        fieldtens0(iatom,imol) = fieldtens0(iatom,imol) + fprefac * sinkrij_dip_sum

                        d = fprefac * (coskrij_dip_sum - sinkrij_charge_sum)
                        
                        do ii = 1,3
                           fieldtens1(ii,iatom,imol) = fieldtens1(ii,iatom,imol) + k_vec(ii) * d
                        end do 
                     endif

                     if(max_rank.ge.2) then 
!     quadrupoles
                        quadi_dd_kk = quad_dd_kk(iatom,imol)

                        sinkrij_quad_sum = sinkri*cos_quadi_dd_kk_sum - coskri*sin_quadi_dd_kk_sum
                        coskrij_quad_sum = coskri*cos_quadi_dd_kk_sum + sinkri*sin_quadi_dd_kk_sum

                        f = f - fprefac * (chargei * sinkrij_quad_sum + quadi_dd_kk * sinkrij_charge_sum)
                        f = f + fprefac * (quadi_dd_kk * coskrij_dip_sum - dipi_d_k * coskrij_quad_sum &
     &                                   + quadi_dd_kk * sinkrij_quad_sum)

                        fieldtens0(iatom,imol) = fieldtens0(iatom,imol) - fprefac * coskrij_quad_sum

                        a = fprefac * sinkrij_quad_sum                        
                        do ii = 1,3
                           fieldtens1(ii,iatom,imol) = fieldtens1(ii,iatom,imol) + k_vec(ii) * a
                        end do 

                        a = fprefac * (coskrij_charge_sum + sinkrij_dip_sum - coskrij_quad_sum)
                        do ii = 1,5
                           fieldtens2(ii,iatom,imol) = fieldtens2(ii,iatom,imol) - ksolid2(ii) * a
                        end do 

                     endif
                  
                     if(max_rank.ge.3) then 
!     octopoles
                        octi_ddd_kkk = oct_ddd_kkk(iatom,imol)
                        sinkrij_oct_sum = sinkri*cos_octi_ddd_kkk_sum - coskri*sin_octi_ddd_kkk_sum
                        coskrij_oct_sum = coskri*cos_octi_ddd_kkk_sum + sinkri*sin_octi_ddd_kkk_sum

                        fieldtens0(iatom,imol) = fieldtens0(iatom,imol) - fprefac * sinkrij_oct_sum

                        a = - fprefac * coskrij_oct_sum
                        do ii = 1,3
                              fieldtens1(ii,iatom,imol) = fieldtens1(ii,iatom,imol) + a * k_vec(ii) 
                        end do 

                        a = fprefac * sinkrij_oct_sum
                        do ii = 1,5
                           fieldtens2(ii,iatom,imol) = fieldtens2(ii,iatom,imol) + ksolid2(ii) * a
                        end do 

                        a = (sinkrij_charge_sum - coskrij_dip_sum - sinkrij_quad_sum + coskrij_oct_sum) * fprefac
                        do ii = 1,7
                           fieldtens3(ii,iatom,imol) = fieldtens3(ii,iatom,imol) + ksolid3(ii) * a
                        end do 

                        f = f + fprefac * (- octi_ddd_kkk * coskrij_charge_sum &
     &                             + chargei * coskrij_oct_sum & 
     &                       - octi_ddd_kkk * sinkrij_dip_sum  & 
     &                            - dipi_d_k * sinkrij_oct_sum & 
     &                       + octi_ddd_kkk * coskrij_quad_sum &
     &                        - quadi_dd_kk * coskrij_oct_sum & 
     &                       + octi_ddd_kkk * sinkrij_oct_sum)
                     endif

                     force(1,iatom,imol) = force(1,iatom,imol) + f * kx * fac(2) 
                     force(2,iatom,imol) = force(2,iatom,imol) + f * ky * fac(2) 
                     force(3,iatom,imol) = force(3,iatom,imol) + f * kz * fac(2) 
                     
                     rcx = lab_relative_coord(1,iatom,imol)
                     rcy = lab_relative_coord(2,iatom,imol)
                     rcz = lab_relative_coord(3,iatom,imol)
                     
                     k_d_rc = (kx * rcx + ky * rcy + kz * rcz)*twopii
                     
                     do jvec = 1,3
                        do ivec = 1,jvec
                           dk_d_rdcvec = rcx * dkdcvec(1,ivec,jvec) + rcy * dkdcvec(2,ivec,jvec) + rcz * dkdcvec(3,ivec,jvec)
                           dudcvec_intra(ivec,jvec) = dudcvec_intra(ivec,jvec) - f * dk_d_rdcvec * fac(2) 
                        end do 
                     end do 

               end do 
            end do iloop

            sum_real = cos_charge_sum 
            sum_imag = -sin_charge_sum 
            
            if(max_rank.ge.1) then
               sum_real = sum_real - sin_dipi_d_k_sum 
               sum_imag = sum_imag - cos_dipi_d_k_sum
            endif
            if(max_rank.ge.2) then 
               sum_real = sum_real - cos_quadi_dd_kk_sum
               sum_imag = sum_imag + sin_quadi_dd_kk_sum
            endif
            if(max_rank.ge.3) then 
               sum_real = sum_real + sin_octi_ddd_kkk_sum
               sum_imag = sum_imag + cos_octi_ddd_kkk_sum
            endif
            sum2 = sum_real**2 + sum_imag**2
               
!     reciprocal Ewald Sum contribution for this k vector
            sum_ewald = sum_ewald + 0.5d0 * fprefac * fac(2) * sum2
            
            do jvec = 1,3
               do ivec = 1,jvec
                  kdkdcvec = k_vec(1) * dkdcvec(1,ivec,jvec) &
     &                     + k_vec(2) * dkdcvec(2,ivec,jvec) &
     &                     + k_vec(3) * dkdcvec(3,ivec,jvec)
                  dudcvec_coulomb(ivec,jvec) = dudcvec_coulomb(ivec,jvec) + 0.5d0 * kdkdcvec * dfprefac * sum2 * fac(2) 
               end do 
            end do 

            do jvec = 1,3
               do ivec = 1,jvec
                  du = 0.0d0 
                  if(max_rank.ge.1) then 
                     du = du &
     &                    - (sum_real * sin_ddip_d_k_dcvec_sum(ivec,jvec) &
     &                    + sum_imag * cos_ddip_d_k_dcvec_sum(ivec,jvec))  
                  endif
                  if(max_rank.ge.2) then 
                     du = du & 
     &                    - (sum_real * cos_dquad_dd_kk_dcvec_sum(ivec,jvec) &
     &                    - sum_imag * sin_dquad_dd_kk_dcvec_sum(ivec,jvec)) 
                  endif
                  if(max_rank.ge.3) then 
                     du = du & 
     &                    + (sum_real * sin_doct_ddd_kkk_dcvec_sum(ivec,jvec) &
     &                    + sum_imag * cos_doct_ddd_kkk_dcvec_sum(ivec,jvec)) 
                  endif
                     dudcvec_coulomb(ivec,jvec) = dudcvec_coulomb(ivec,jvec) + du * fprefac * fac(2) 
               end do 
            end do 

         endif
      endif
      end do 
      end do 
      end do 
      
!     write(*,*) 'Reciprocal space energy = ',sum_ewald*pefac
      
      dudcvec_coulomb(1,1) = dudcvec_coulomb(1,1) - (by*cz - bz*cy) * voli * sum_ewald 

      dudcvec_coulomb(1,2) = dudcvec_coulomb(1,2) - (cy*az - cz*ay) * voli * sum_ewald 
      dudcvec_coulomb(2,2) = dudcvec_coulomb(2,2) - (cz*ax - cx*az) * voli * sum_ewald 

      dudcvec_coulomb(1,3) = dudcvec_coulomb(1,3) - (ay*bz - az*by) * voli * sum_ewald 
      dudcvec_coulomb(2,3) = dudcvec_coulomb(2,3) - (az*bx - ax*bz) * voli * sum_ewald 
      dudcvec_coulomb(3,3) = dudcvec_coulomb(3,3) - (ax*by - ay*bx) * voli * sum_ewald 

      end subroutine reciprocal
      
      subroutine energy()
!     calculates the energy of the system
      use common_data
      implicit none
      integer natoms_in_mol,imoltype,imol,iatom

      upair = 0.0d0 
      upv = 0.0d0 
      uconformer = 0.0d0 
      
      dudcvec_coulomb = 0.0d0 
      dudcvec_intra = 0.0d0 
      dudcvec_pair = 0.0d0 
      dudcvec_pv = 0.0d0 
      dudcvec = 0.0d0 

      call initialize

      if(periodic) then
         call reciprocal
         call self
      endif
      call real
      call calc_elec
      call calc_conformer_energy

      if(periodic) then 
         call pair_longrange
         call press
      endif

      dudcvec = dudcvec_coulomb + dudcvec_pair + dudcvec_intra + dudcvec_pv
      utotal = utotal + ucharge + udipole + uquad + uoct + upair + upv + uconformer

      end subroutine energy

      subroutine real()
!     calculates the real part of the energy of the system
      use common_data
      use nr_mod
      implicit none
      integer imol,jmol,n
      integer natoms_in_imol,natoms_in_jmol
      integer iatom,jatom
      integer imoltype,jmoltype
      integer :: jcell1,jcell2,jcell3
      integer :: jcell1max,jcell2max,jcell3max
      integer nmax,irank,jrank
      real(8) :: xi,yi,zi,xj,yj,zj
      real(8), dimension(3)  :: dipi,dipj
      real(8) :: dipi_d_r,dipj_d_r
      real(8) :: xdif,ydif,zdif,rdis
      real(8) :: sigma,sigma6,sigma12,epsilon,dupair
      real(8) :: ri,r2i,r3i,r4i,r6i,r12i
      real(8) :: charge_i,charge_j
      real(8) :: fx,fy,fz
      real(8) :: chargei_chargej,diff
      real(8) :: xforce_pair,yforce_pair,zforce_pair
      real(8) :: xforce,yforce,zforce
      real(8), dimension(3) :: force_elec
      real(8), dimension(0:10) :: pfunc,bfunc,gfunc,tfunc
      real(8) :: rdis2
      real(8) :: xxdif,yydif,zzdif
      real(8) :: upair0,switch,rc0,rc1,dupair0,dswitchdss,dswitchdr,ss,dssdr
      real(8) :: diff0,cfield0i,cfield0j
      real(8) :: sfac
      real(8) :: imolx,imoly,imolz,jmolx,jmoly,jmolz
      logical :: same_mol,same_cell,same_atom,loop,okflag,match
      real(8), dimension(3,3) :: scvec,scveci
      real(8), dimension(3) :: vec1,vec2
      integer, dimension(3) :: vec3
      real(8) :: f1dif,f2dif,f3dif
      real(8) :: f1,f2,f3
      real(8) :: ax,ay,az,bx,by,bz,cx,cy,cz
      real(8) :: vol
      real(8) :: xj0,yj0,zj0,xdif0,ydif0,zdif0
      real(8) :: dcellx,dcelly,dcellz,dcellx0,dcelly0,dcellz0
      real(8) :: xxdif0,yydif0,zzdif0
      real(8) :: ffix,ffiy,ffiz,ffjx,ffjy,ffjz
      real(8) :: dcyz,dcxz,dcxy
      real(8) :: astar_mag,bstar_mag,cstar_mag
      real(8) :: pair6,pair12
      real(8) :: dupair0dp1,dupair0dp2
      real(8) :: dupairdp1,dupairdp2
      real(8) :: rspline_max
      real(8) :: dipi_d_dipj
      real(8) :: x,y,z,xf,yf,zf
      real(8) :: quadi_dd_rr,quadj_dd_rr,quadi_d_dipj_d_r,quadj_d_dipi_d_r
      real(8) :: quadi_dd_quadj,quadi_d_r_d_quadj_r_r
      real(8), dimension(3) :: quadi_d_r,quadj_d_r,quadi_d_dipj,quadj_d_dipi
      real(8), dimension(3,3) :: quadi,quadj
      real(8), dimension(3,3) :: quadi_d_quadj,quadj_d_quadi
      real(8), dimension(3) :: quadi_d_quadj_d_r,quadj_d_quadi_d_r
      real(8), dimension(3,3,3) :: octi,octj
      real(8), dimension(3) :: dif
      real(8) :: rrrr,rrr,rr
      real(8) :: octi_ddd_rrr,octj_ddd_rrr
      real(8) :: di_d_pj_dd_r,dj_d_pi_dd_r
      real(8), dimension(3) :: octi_dd_rr,octj_dd_rr
      real(8), dimension(3) :: octi_d_r_d_dipj,octj_d_r_d_dipj
      real(8) :: dioctj_dddd_rrrr,djocti_dddd_rrrr
      real(8), dimension(3) :: quadi_d_octj_dd_rr,quadj_d_octi_dd_rr
      real(8) :: dipi_d_octj_dd_rr,dipj_d_octi_dd_rr
      real(8), dimension(3) :: quadi_d_r_d_octj_d_r,quadj_d_r_d_octi_d_r
      real(8) :: quadi_d_r_d_octj_dd_rr,quadj_d_r_d_octi_dd_rr
      real(8), dimension(3) :: quadj_dd_octi,quadi_dd_octj
      real(8) :: quadj_dd_octi_d_r,quadi_dd_octj_d_r
      real(8), dimension(3) :: octi_d_r_d_octj_dd_rr,octj_d_r_d_octi_dd_rr
      real(8) :: octi_dd_rr_d_octj_dd_rr
      real(8), dimension(3) :: octi_dd_octj_d_r,octj_dd_octi_d_r
      real(8) :: octi_d_r_dd_octj_d_r
      real(8) :: octi_ddd_octj
      real(8) :: octi_ijk,octj_ijk
      real(8), dimension(0:12) :: ggg
      real(8), dimension(5) :: spherical2i,spherical2j
      real(8), dimension(7) :: spherical3i,spherical3j
      real(8), dimension(5) :: octi_d_r,octj_d_r
      real(8), dimension(3,3) :: octi_d_r1_1,octj_d_r1_1
      real(8), dimension(5,3) :: spherical2_1i,spherical2_1j
      integer :: ivec,jvec,kvec,lvec,ii,jj
      integer :: i_index,j_index
      integer :: imoltype_list,jmoltype_list
      integer :: i_index_list,j_index_list
      integer :: m,iparam,nterms,j,nb,mm
      real(8) :: xx,xy,xz,yy,yz,zz
      real(8) :: xxx,xxy,xyy,yyy,xxz,xyz,yyz,xzz,yzz,zzz
      real(8), dimension(5) :: rsolid2
      real(8), dimension(7) :: rsolid3
      real(8), dimension(5) :: solid2
      real(8), dimension(3,3) :: solid1_1

      integer mmm

      rc0 = rcut * 0.9d0 
      rc1 =  rcut

      sfac = 1.0d0 / (rc1 - rc0)

      if(periodic) then 

!     finds the size of the supercell which contains the cut-off sphere
         vol = cvec(1,1) * cvec(2,2) * cvec(3,3)

!     magnitude of the reciprocal lattice vectors
         astar_mag = dsqrt(cveci(1,1)**2 + cveci(1,2)**2 + cveci(1,3)**2)
         bstar_mag = dsqrt(cveci(2,1)**2 + cveci(2,2)**2 + cveci(2,3)**2)
         cstar_mag = dsqrt(cveci(3,1)**2 + cveci(3,2)**2 + cveci(3,3)**2)

         jcell1max = int(2.0d0 * rcut * astar_mag) + 1
         jcell2max = int(2.0d0 * rcut * bstar_mag) + 1
         jcell3max = int(2.0d0 * rcut * cstar_mag) + 1

!     supercell cell-vectors

         scvec(:,1) = (jcell1max) * cvec(:,1)
         scvec(:,2) = (jcell2max) * cvec(:,2)
         scvec(:,3) = (jcell3max) * cvec(:,3)
         
         call mat3inverse(scvec,scveci,okflag)      
      else
         jcell1max = 1
         jcell2max = 1
         jcell3max = 1
      endif

      imol_loop: do imol = 1,nmol
         imoltype = mol_type(imol)
         natoms_in_imol = mol_natoms(imoltype)
         iatom_loop: do iatom = 1,natoms_in_imol

         xi = lab_coord(1,iatom,imol)
         yi = lab_coord(2,iatom,imol)
         zi = lab_coord(3,iatom,imol)

         irank = site_rank(iatom,imoltype)

         i_index = pair_index(iatom,imoltype)

         charge_i = lab_charge(iatom,imol)

         if(irank.ge.1) then 
            dipi(:) = lab_dipole(:,iatom,imol) 
         else
            dipi = 0.0d0 
         endif
         if(irank.ge.2) then 
            quadi(:,:) = lab_quad(:,:,iatom,imol)
            do ii = 1,5
               spherical2i(ii) = lab_spherical2(ii,iatom,imol)
            end do 
         else
            quadi = 0.0d0 
            spherical2i = 0.0d0 
         endif
         if(irank.ge.3) then 
            do ii = 1,7
               spherical3i(ii) = lab_spherical3(ii,iatom,imol)
            end do 
            do ii = 1,5
               do jj = 1,3
                  spherical2_1i(ii,jj) = lab_spherical2_1(ii,jj,iatom,imol)
               end do 
            end do 
         else
            do ii = 1,7
               spherical3i(ii) = 0.0d0 
            end do 
            do ii = 1,5
               do jj = 1,3
                  spherical2_1i(ii,jj) = 0.0d0 
               end do 
            end do 
         endif

         imolx = mol_com(1,imol)
         imoly = mol_com(2,imol)
         imolz = mol_com(3,imol)

         jmol_loop: do jmol = imol,nmol
         jmoltype = mol_type(jmol)
         natoms_in_jmol = mol_natoms(jmoltype)
         jatom_loop: do jatom = 1,natoms_in_jmol

         jrank = site_rank(jatom,jmoltype)
         
         j_index = pair_index(jatom,jmoltype)
         match = .false.
         mm = 0 
         do m = 1,nbondlist
            i_index_list = model_list(m,1)
            j_index_list = model_list(m,2)
            if(i_index.eq.i_index_list.and.j_index.eq.j_index_list) match = .true.
            if(j_index.eq.i_index_list.and.i_index.eq.j_index_list) match = .true.
            if(match) then
               nb = m
               exit
            endif
         end do 

         if(nb.ne.0) then 
            nterms = model_nterms(nb)
            rspline_max = rsplinemax(nb)
         else
            nterms = 0
            rspline_max = 0.0d0
         endif

         charge_j = lab_charge(jatom,jmol)

         if(jrank.ge.1) then 
            dipj(:) = lab_dipole(:,jatom,jmol) 
         else
            dipj = 0.0d0 
         endif
         if(jrank.ge.2) then 
            quadj(:,:) = lab_quad(:,:,jatom,jmol)
            do ii = 1,5
               spherical2j(ii) = lab_spherical2(ii,jatom,jmol)
            end do 
         else
            quadj = 0.0d0 
            spherical2j = 0.0d0 
         endif
         if(jrank.ge.3) then 
            do ii = 1,7
               spherical3j(ii) = lab_spherical3(ii,jatom,jmol)
            end do 
            do ii = 1,5
               do jj = 1,3
                  spherical2_1j(ii,jj) = lab_spherical2_1(ii,jj,jatom,jmol)
               end do 
            end do 
         else
            do ii = 1,7
               spherical3j(ii) = 0.0d0 
            end do 
            do ii = 1,5
               do jj = 1,3
                  spherical2_1j(ii,jj) = 0.0d0 
               end do 
            end do 
         endif

         jmolx = mol_com(1,jmol)
         jmoly = mol_com(2,jmol)
         jmolz = mol_com(3,jmol)

         if(periodic) then
            xj0 = lab_coord(1,jatom,jmol) 
            yj0 = lab_coord(2,jatom,jmol) 
            zj0 = lab_coord(3,jatom,jmol) 

            xdif0 = xj0 - xi
            ydif0 = yj0 - yi
            zdif0 = zj0 - zi

            xxdif0 = jmolx - imolx
            yydif0 = jmoly - imoly
            zzdif0 = jmolz - imolz
         else
            xj = lab_coord(1,jatom,jmol) 
            yj = lab_coord(2,jatom,jmol) 
            zj = lab_coord(3,jatom,jmol) 

            xdif = xj - xi
            ydif = yj - yi
            zdif = zj - zi
         endif


         same_atom = .false.
         if(imol.eq.jmol.and.iatom.eq.jatom) same_atom = .true.

         ffix = 0.0d0
         ffiy = 0.0d0 
         ffiz = 0.0d0 

         ffjx = 0.0d0 
         ffjy = 0.0d0 
         ffjz = 0.0d0 

!     PRECALCULATE MULTIPOLE DOT PRODUCTS BEFORE WE ENTER THE LOOP OVER BOXES

         chargei_chargej = charge_i * charge_j

         if(irank.ge.1.or.jrank.ge.1) then 
!     dipole precalculations

            dipi_d_dipj = 0.0d0 
            do ivec = 1,3
               dipi_d_dipj = dipi_d_dipj + dipi(ivec) * dipj(ivec)
            end do 
         endif
         if(irank.ge.2.or.jrank.ge.2) then 
!     quadrupole precalculations

            quadi_d_quadj = 0.0d0 
            quadj_d_quadi = 0.0d0 

            do ivec = 1,3
               do jvec = 1,3
                  do kvec = 1,3
                     quadi_d_quadj(ivec,jvec) = quadi_d_quadj(ivec,jvec) + quadi(ivec,kvec) * quadj(kvec,jvec)
                     quadj_d_quadi(ivec,jvec) = quadj_d_quadi(ivec,jvec) + quadj(ivec,kvec) * quadi(kvec,jvec)
                  end do 
               end do 
            end do 

            quadi_d_dipj = 0.0d0 
            quadj_d_dipi = 0.0d0 
            do ivec = 1,3
               do jvec = 1,3
                  quadi_d_dipj(ivec) = quadi_d_dipj(ivec) + quadi(ivec,jvec) * dipj(jvec)
                  quadj_d_dipi(ivec) = quadj_d_dipi(ivec) + quadj(ivec,jvec) * dipi(jvec)
               end do 
            end do 
         endif

         quadi_dd_quadj = 0.0d0 

         do ii = 1,5
            quadi_dd_quadj = quadi_dd_quadj + spherical2i(ii) * spherical2j(ii)
         end do 

         if(irank.ge.3.or.jrank.ge.3) then 
!     octopole precalculations

            quadj_dd_octi = 0.0d0 
            quadi_dd_octj = 0.0d0 
            
            do ii = 1,5
               do jj = 1,3
                  quadj_dd_octi(jj) = quadj_dd_octi(jj) + spherical2_1i(ii,jj)*spherical2j(ii)
                  quadi_dd_octj(jj) = quadi_dd_octj(jj) + spherical2_1j(ii,jj)*spherical2i(ii)
               end do 
            end do 


            octi_ddd_octj = 0.0d0 
            do ii = 1,7
               octi_ddd_octj = octi_ddd_octj + spherical3i(ii) * spherical3j(ii)
            end do 


         endif

         gfunc = 0.0d0 
         do jcell3 = 0,jcell3max-1
            dcellz = jcell3 * cvec(3,3)
            dcyz = jcell3 * cvec(2,3)
            dcxz = jcell3 * cvec(1,3)
            do jcell2 = 0,jcell2max-1
               dcelly = jcell2 * cvec(2,2) + dcyz
               dcxy = jcell2 * cvec(1,2)
               do jcell1 = 0,jcell1max-1
                  dcellx = jcell1 * cvec(1,1) + dcxy + dcxz

!     decide whether to loop
                  loop = .false.
                  same_mol = .false.
                  
                  if(periodic) then
                     same_cell = .false.
                     if(jcell1.eq.0.and.jcell2.eq.0.and.jcell3.eq.0) same_cell = .true.
                     if(same_cell.and.(imol.eq.jmol)) same_mol = .true.
                     
                     if((jmol.gt.imol).or.(imol.eq.jmol.and.jatom.gt.iatom)) loop = .true.

                     if(same_atom) then
                        loop = .false.
                        if(.not. same_cell) loop = .true.
                     endif
                  else
                     if(jmol.gt.imol) loop = .true.
                  endif
                  if(loop) then


!     find minimum image in supercell

                     if(periodic) then
                        xdif = xdif0 + dcellx
                        ydif = ydif0 + dcelly
                        zdif = zdif0 + dcellz
                        
                        call nearest(xdif,ydif,zdif,scvec,scveci)

!     If same-atom, only count nearest neighbors in half-sphere, to avoid double-counting. 
!     Don't follow through with this loop unless xdif > 0
                        if(same_atom) then
                           if((xdif.gt.0).or.(xdif.eq.0.and.ydif.gt.0)&
     &                          .or.(xdif.eq.0.and.ydif.eq.0.and.zdif.gt.0)) then 
                        else
                           cycle
                        endif
                     endif
                     
                  endif
                  
                  rdis2 = xdif**2 + ydif**2 + zdif**2
                  if(periodic.and.(rdis2.le.rcut2).or.(.not.periodic)) then

                  dif(1) = xdif
                  dif(2) = ydif
                  dif(3) = zdif

                  rdis = dsqrt(rdis2)                  
!                  if(iatom.eq.1.and.jatom.eq.1) write(38,*) rdis
                  
                  if(rdis.eq.0) then 
                     write(*,*) 'RDIS = 0, STOPPING',imol,jmol,iatom,jatom&
     &                    ,jcell1,jcell2,jcell3,same_cell,same_mol
                     write(*,*) xi,yi,zi
                     write(*,*) xj0+dcellx,yj0+dcelly,zj0+dcellz
                     stop
                  endif

                  ri = 1.0d0 / rdis
                  r2i = ri * ri

                  nmax = irank + jrank + 1 

                  if(periodic) call getbfuncs(ewald_eps,eps_sqrtpii,rdis,ri,gfunc,nmax)

                  if(periodic) then
                     if(same_mol) then
                        call getpfuncs(ri,r2i,pfunc,nmax)
                        gfunc = gfunc - pfunc
                     endif
                  else
                     call getpfuncs(ri,r2i,pfunc,nmax)
                     gfunc = pfunc
                  endif
                  
!     short-range electrostatic damping 
                  if(.not.same_mol.and.rdis.lt.rdamp_cutoff) then
                     call getbfuncs(eps_damp,eps_damp_sqrtpii,rdis,ri,tfunc,nmax)                  
                     gfunc = gfunc - tfunc
                  endif

                  if(periodic) then
!     calculate switch function
                     ss = (rc1 - rdis) * sfac
                     dssdr = -sfac
                     
                     switch = ss**2 * (3.0d0 - 2.0d0 * ss)
                     dswitchdss = 6.0d0 * ss * (1.0d0  - ss)
                     
!     switch = 6.0d0 * ss**5 - 15.0d0  * ss**4 + 10.0d0 * ss**3
!     dswitchdss = 30.0d0 * ss**4 - 60.0d0  * ss**3 + 30.0d0 * ss**2
                     
                     dswitchdr = dswitchdss * dssdr
                  endif
                  
                  dupair = 0.0d0 
                  xforce_pair = 0.0d0
                  yforce_pair = 0.0d0 
                  zforce_pair = 0.0d0 

                  if(.not.same_mol) then
                     if(use_spline.and.rdis.le.rspline_max) then
                        call calc_pair_from_spline(rdis,ri,dupair,rc0,rc1,switch,dswitchdr,nb)
                     else
                        call calc_pair_from_polynomial(rdis,ri,dupair,rc0,rc1,switch,dswitchdr,nb,nterms)
                     endif

                     xforce_pair = - dupair * xdif 
                     yforce_pair = - dupair * ydif 
                     zforce_pair = - dupair * zdif                      
                     
                  endif
                  
!     ELECTROSTATICS

                  force_elec = 0.0d0 

                  do n = 0,irank + jrank
                     ggg(n) = 0.0d0 
                  end do 

                  fieldtens0(iatom,imol) = fieldtens0(iatom,imol) + charge_j * gfunc(0) 
                  fieldtens0(jatom,jmol) = fieldtens0(jatom,jmol) + charge_i * gfunc(0) 

!     ggg are the Smith Gji(r) functions 

                  ggg(0) = ggg(0) + chargei_chargej

                  if(irank.ge.1.or.jrank.ge.1) then 
!     DIPOLE interactions

                     dipi_d_r = 0.0d0 
                     dipj_d_r = 0.0d0 
                     do ivec = 1,3
                        dipi_d_r = dipi_d_r + dipi(ivec) * dif(ivec)
                        dipj_d_r = dipj_d_r + dipj(ivec) * dif(ivec)
                     end do 

                     fieldtens0(iatom,imol) = fieldtens0(iatom,imol) - dipj_d_r * gfunc(1)
                     fieldtens0(jatom,jmol) = fieldtens0(jatom,jmol) + dipi_d_r * gfunc(1)

                     do ii = 1,3
                        fieldtens1(ii,iatom,imol) = fieldtens1(ii,iatom,imol) & 
     &                       + (dif(ii) * charge_j + dipj(ii)) * gfunc(1) &
     &                       - dif(ii) * dipj_d_r * gfunc(2)  

                        fieldtens1(ii,jatom,jmol) = fieldtens1(ii,jatom,jmol) & 
     &                       + (- dif(ii) * charge_i + dipi(ii)) * gfunc(1) &
     &                       - dif(ii) * dipi_d_r * gfunc(2)  
                     end do 


                     ggg(1) = ggg(1) + dipi_d_r * charge_j - dipj_d_r * charge_i  & 
     &                    + dipi_d_dipj
                     ggg(2) = ggg(2) - dipi_d_r * dipj_d_r 
                     do ivec = 1,3
                        force_elec(ivec) = force_elec(ivec) &
     &                       + (dipi(ivec) * charge_j - dipj(ivec) * charge_i) * gfunc(1) & 
     &                       - (dipi(ivec) * dipj_d_r + dipj(ivec) * dipi_d_r) * gfunc(2) 
                     end do 
                  endif


                  if(irank.ge.2.or.jrank.ge.2) then 
!     QUADRUPOLE INTERACTIONS

                  xx = xdif * xdif
                  xy = xdif * ydif
                  xz = xdif * zdif
                  yy = ydif * ydif
                  yz = ydif * zdif
                  zz = zdif * zdif
                  
                  call convert_quad_to_spherical(xx,xy,xz,yy,yz,zz,rsolid2)


                  quadi_d_r = 0.0d0 
                  quadj_d_r = 0.0d0 
                  
                  do ivec = 1,3
                     do jvec = 1,3
                        quadi_d_r(ivec) = quadi_d_r(ivec) + quadi(ivec,jvec)*dif(jvec)
                        quadj_d_r(ivec) = quadj_d_r(ivec) + quadj(ivec,jvec)*dif(jvec)
                     end do 
                  end do 
                  
                  quadi_d_dipj_d_r = 0.0d0 
                  quadj_d_dipi_d_r = 0.0d0 
                  quadi_d_r_d_quadj_r_r = 0.0d0 

                  do ivec = 1,3
                     quadi_d_dipj_d_r = quadi_d_dipj_d_r - quadi_d_dipj(ivec) * dif(ivec)
                     quadj_d_dipi_d_r = quadj_d_dipi_d_r + quadj_d_dipi(ivec) * dif(ivec)
                     quadi_d_r_d_quadj_r_r = quadi_d_r_d_quadj_r_r + quadi_d_r(ivec) * quadj_d_r(ivec) 
                  end do 
                  
                  quadi_dd_rr = 0.0d0 
                  quadj_dd_rr = 0.0d0 
                  do ii = 1,5
                     quadj_dd_rr = quadj_dd_rr + spherical2j(ii) * rsolid2(ii)
                     quadi_dd_rr = quadi_dd_rr + spherical2i(ii) * rsolid2(ii)
                  end do 
                  

                  quadi_d_quadj_d_r = 0.0d0 
                  quadj_d_quadi_d_r = 0.0d0 
                  
                  do ivec = 1,3
                     do jvec = 1,3
                        quadi_d_quadj_d_r(ivec) = quadi_d_quadj_d_r(ivec) + quadi_d_quadj(ivec,jvec) * dif(jvec)
                        quadj_d_quadi_d_r(ivec) = quadj_d_quadi_d_r(ivec) + quadj_d_quadi(ivec,jvec) * dif(jvec)
                     end do 
                  end do 
                  
                  fieldtens0(iatom,imol) = fieldtens0(iatom,imol) + quadj_dd_rr * gfunc(2)
                  fieldtens0(jatom,jmol) = fieldtens0(jatom,jmol) + quadi_dd_rr * gfunc(2)
                  
                  do ii = 1,3
                     fieldtens1(ii,iatom,imol) = fieldtens1(ii,iatom,imol) & 
     &                    - 2.0d0 * quadj_d_r(ii) * gfunc(2) &
     &                    + dif(ii) * quadj_dd_rr * gfunc(3) 
                     fieldtens1(ii,jatom,jmol) = fieldtens1(ii,jatom,jmol) &
     &                    + 2.0d0 * quadi_d_r(ii) * gfunc(2) &
     &                    - dif(ii) * quadi_dd_rr * gfunc(3) 
                  end do 
                  
                  do ii = 1,3
                     do jj = 1,3
                        fieldtens1_1(ii,jj,iatom,imol) = fieldtens1_1(ii,jj,iatom,imol)  &
     &                       + 2.0d0 * dipj(ii) * dif(jj) * gfunc(2) & 
     &                       - 4.0d0 * quadj_d_r(ii) * dif(jj) * gfunc(3) 
                        fieldtens1_1(ii,jj,jatom,jmol) = fieldtens1_1(ii,jj,jatom,jmol)  &
     &                       - 2.0d0 * dipi(ii) * dif(jj) * gfunc(2) & 
     &                       - 4.0d0 * quadi_d_r(ii) * dif(jj) * gfunc(3) 
                     end do 
                  end do 
                  
                  do ii = 1,5
                     fieldtens2(ii,iatom,imol) = fieldtens2(ii,iatom,imol) &
     &                    + (rsolid2(ii) * charge_j & 
     &                    +2.0d0 * spherical2j(ii)) * gfunc(2) & 
     &                    -rsolid2(ii) * dipj_d_r * gfunc(3) & 
     &                    +rsolid2(ii) * quadj_dd_rr * gfunc(4)
                     fieldtens2(ii,jatom,jmol) = fieldtens2(ii,jatom,jmol) & 
     &                    + (rsolid2(ii) * charge_i & 
     &                    +2.0d0 * spherical2i(ii)) * gfunc(2) & 
     &                    +rsolid2(ii) * dipi_d_r * gfunc(3) &
     &                    +rsolid2(ii) * quadi_dd_rr * gfunc(4)
                  end do 
                  
                  
                  ggg(2) = ggg(2) + charge_j * quadi_dd_rr + charge_i * quadj_dd_rr &
     &                 - 2.0d0 * (quadi_d_dipj_d_r + quadj_d_dipi_d_r) &
     &                 +2.0d0 * quadi_dd_quadj 
                  
                  ggg(3) = ggg(3) +  dipi_d_r * quadj_dd_rr - dipj_d_r * quadi_dd_rr & 
     &                               - 4.0d0 * quadi_d_r_d_quadj_r_r
                  
                  ggg(4) = ggg(4) + quadi_dd_rr * quadj_dd_rr 

                  do ivec = 1,3
                     force_elec(ivec) = force_elec(ivec) &
     &                    + (2.0d0 * (charge_j * quadi_d_r(ivec) + charge_i * quadj_d_r(ivec)) &
     &                    - 2.0d0 * (quadj_d_dipi(ivec) - quadi_d_dipj(ivec))) * gfunc(2)  & 
                     
     &                    + (- 2.0d0 * (-dipi_d_r * quadj_d_r(ivec) + dipj_d_r * quadi_d_r(ivec)) &
     &                    - (quadi_dd_rr * dipj(ivec) - quadj_dd_rr * dipi(ivec)) &
     &                    - 4.0d0 * (quadi_d_quadj_d_r(ivec) + quadj_d_quadi_d_r(ivec))) * gfunc(3) & 
                     
     &                    + 2.0d0 * (quadi_dd_rr * quadj_d_r(ivec) + quadj_dd_rr * quadi_d_r(ivec)) * gfunc(4) 
                  end do 
               endif
               
               if(irank.ge.3.or.jrank.ge.3) then 
!     OCTOPOLE INTERACTIONS
                  
                  xxx = xx * xdif
                  xxy = xx * ydif
                  xyy = yy * xdif
                  yyy = yy * ydif
                  xxz = xx * zdif
                  xyz = xy * zdif
                  yyz = yy * zdif
                  xzz = zz * xdif
                  yzz = zz * ydif
                  zzz = zz * zdif

                  call  convert_oct_to_spherical(xxx,xxy,xyy,yyy,xxz,xyz,yyz,xzz,yzz,zzz,rsolid3)

                  do ii = 1,3
                     octi_dd_rr(ii) = 0.0d0 
                     octj_dd_rr(ii) = 0.0d0 
                     do jj = 1,5
                        octi_dd_rr(ii) = octi_dd_rr(ii) + spherical2_1i(jj,ii) * rsolid2(jj)
                        octj_dd_rr(ii) = octj_dd_rr(ii) + spherical2_1j(jj,ii) * rsolid2(jj)
                     end do 
                  end do 
                  
                  
                  do ii = 1,5
                     octi_d_r(ii) = 0.0d0 
                     octj_d_r(ii) = 0.0d0 
                     do jj = 1,3
                        octi_d_r(ii) = octi_d_r(ii) + spherical2_1i(ii,jj) * dif(jj)
                        octj_d_r(ii) = octj_d_r(ii) + spherical2_1j(ii,jj) * dif(jj)
                     end do 
                  end do 
                  
                  
!     split octi_d_r
                  call convert_quad_to_cartesian(octi_d_r,xx,xy,xz,yy,yz,zz)
                  call make_cart2(xx,xy,xz,yy,yz,zz,octi_d_r1_1)
                  call convert_quad_to_cartesian(octj_d_r,xx,xy,xz,yy,yz,zz)
                  call make_cart2(xx,xy,xz,yy,yz,zz,octj_d_r1_1)
                  
                  octi_d_r_d_dipj = 0.0d0 
                  octj_d_r_d_dipj = 0.0d0 
                  
                  do ivec = 1,3
                     do jvec = 1,3
                        octi_d_r_d_dipj(ivec) = octi_d_r_d_dipj(ivec) + octi_d_r1_1(ivec,jvec) * dipj(jvec)
                        octj_d_r_d_dipj(ivec) = octj_d_r_d_dipj(ivec) + octj_d_r1_1(ivec,jvec) * dipi(jvec)
                     end do 
                  end do 
                  
                  octi_ddd_rrr = 0.0d0 
                  octj_ddd_rrr = 0.0d0 
                  do ivec = 1,3
                     octi_ddd_rrr = octi_ddd_rrr + octi_dd_rr(ivec) * dif(ivec)
                     octj_ddd_rrr = octj_ddd_rrr + octj_dd_rr(ivec) * dif(ivec)
                  end do 
                  
                  quadi_d_octj_dd_rr = 0.0d0 
                  quadj_d_octi_dd_rr = 0.0d0 
                  quadi_d_r_d_octj_d_r = 0.0d0 
                  quadj_d_r_d_octi_d_r = 0.0d0 
                  octi_d_r_dd_octj_d_r = 0.0d0 
                  octi_d_r_d_octj_dd_rr = 0.0d0 
                  octj_d_r_d_octi_dd_rr = 0.0d0 
                  
                  do ivec = 1,3
                     do jvec = 1,3
                        quadi_d_r_d_octj_d_r(ivec) = quadi_d_r_d_octj_d_r(ivec) + quadi_d_r(jvec) * octj_d_r1_1(ivec,jvec)
                        quadj_d_r_d_octi_d_r(ivec) = quadj_d_r_d_octi_d_r(ivec) + quadj_d_r(jvec) * octi_d_r1_1(ivec,jvec)
                        octi_d_r_dd_octj_d_r = octi_d_r_dd_octj_d_r + octi_d_r1_1(ivec,jvec) * octj_d_r1_1(ivec,jvec)
                        quadi_d_octj_dd_rr(ivec) = quadi_d_octj_dd_rr(ivec) + quadi(ivec,jvec) * octj_dd_rr(jvec)
                        quadj_d_octi_dd_rr(ivec) = quadj_d_octi_dd_rr(ivec) + quadj(ivec,jvec) * octi_dd_rr(jvec)
                        octi_d_r_d_octj_dd_rr(ivec) = octi_d_r_d_octj_dd_rr(ivec) + octi_d_r1_1(ivec,jvec) * octj_dd_rr(jvec) 
                        octj_d_r_d_octi_dd_rr(ivec) = octj_d_r_d_octi_dd_rr(ivec) + octj_d_r1_1(ivec,jvec) * octi_dd_rr(jvec) 
                     end do 
                  end do 
                  
                  
                  quadi_d_r_d_octj_dd_rr = 0.0d0 
                  quadj_d_r_d_octi_dd_rr = 0.0d0 
                  octi_dd_rr_d_octj_dd_rr = 0.0d0
                  dioctj_dddd_rrrr = 0.0d0 
                  djocti_dddd_rrrr = 0.0d0 
                  dipi_d_octj_dd_rr = 0.0d0 
                  dipj_d_octi_dd_rr = 0.0d0 
                  quadj_dd_octi_d_r = 0.0d0 
                  quadi_dd_octj_d_r = 0.0d0 
                  
                  do ivec = 1,3
                     quadi_d_r_d_octj_dd_rr = quadi_d_r_d_octj_dd_rr + quadi_d_r_d_octj_d_r(ivec) * dif(ivec)
                     quadj_d_r_d_octi_dd_rr = quadj_d_r_d_octi_dd_rr + quadj_d_r_d_octi_d_r(ivec) * dif(ivec)
                     octi_dd_rr_d_octj_dd_rr = octi_dd_rr_d_octj_dd_rr + &
     &                    octi_dd_rr(ivec) * octj_dd_rr(ivec)
                     dioctj_dddd_rrrr = dioctj_dddd_rrrr + dipi(ivec) * octj_ddd_rrr * dif(ivec)
                     djocti_dddd_rrrr = djocti_dddd_rrrr + dipj(ivec) * octi_ddd_rrr * dif(ivec)
                     dipi_d_octj_dd_rr = dipi_d_octj_dd_rr + dipi(ivec) * octj_dd_rr(ivec)
                     dipj_d_octi_dd_rr = dipj_d_octi_dd_rr + dipj(ivec) * octi_dd_rr(ivec)
                     quadj_dd_octi_d_r = quadj_dd_octi_d_r + quadj_dd_octi(ivec) * dif(ivec)
                     quadi_dd_octj_d_r = quadi_dd_octj_d_r + quadi_dd_octj(ivec) * dif(ivec)
                  end do 
                  
                  
                  do ii = 1,3
                     octi_dd_octj_d_r(ii) = 0.0d0 
                     octj_dd_octi_d_r(ii) = 0.0d0 
                     do jj = 1,5
                        octi_dd_octj_d_r(ii) = octi_dd_octj_d_r(ii) + spherical2_1i(jj,ii) * octj_d_r(jj)
                        octj_dd_octi_d_r(ii) = octj_dd_octi_d_r(ii) + spherical2_1j(jj,ii) * octi_d_r(jj)
                     end do 
                  end do 
                  
                  
                  fieldtens0(iatom,imol) = fieldtens0(iatom,imol) - gfunc(3) * octj_ddd_rrr
                  fieldtens0(jatom,jmol) = fieldtens0(jatom,jmol) + gfunc(3) * octi_ddd_rrr
                  
                  do ii = 1,3
                     fieldtens1(ii,iatom,imol) = fieldtens1(ii,iatom,imol) + 3.0d0 * octj_dd_rr(ii) * gfunc(3) & 
     &                    - dif(ii) * octj_ddd_rrr * gfunc(4) 
                     fieldtens1(ii,jatom,jmol) = fieldtens1(ii,jatom,jmol) + 3.0d0 * octi_dd_rr(ii) * gfunc(3) & 
     &                    - dif(ii) * octi_ddd_rrr * gfunc(4) 
                  end do 
                  
                  do ii = 1,5
                     fieldtens2(ii,iatom,imol) = fieldtens2(ii,iatom,imol) - octj_ddd_rrr * rsolid2(ii) * gfunc(5) & 
     &                    - 6.0d0 * octj_d_r(ii) * gfunc(3)
                     fieldtens2(ii,jatom,jmol) = fieldtens2(ii,jatom,jmol) + octi_ddd_rrr * rsolid2(ii) * gfunc(5) & 
     &                    + 6.0d0 * octi_d_r(ii) * gfunc(3)
                  end do 
                  
                  do ii = 1,3
                     do jj = 1,3
                        fieldtens1_1(ii,jj,iatom,imol) = fieldtens1_1(ii,jj,iatom,imol) &
     &                       + 6.0d0 * octj_dd_rr(ii) * dif(jj) * gfunc(4)
                        fieldtens1_1(ii,jj,jatom,jmol) = fieldtens1_1(ii,jj,jatom,jmol) & 
     &                       - 6.0d0 * octi_dd_rr(ii) * dif(jj) * gfunc(4)
                     end do 
                  end do 
                  
                  do ii = 1,7
                     fieldtens3(ii,iatom,imol) = fieldtens3(ii,iatom,imol) &
     &                    + (charge_j * rsolid3(ii) & 
     &                    + 6.0d0 * spherical3j(ii)) * gfunc(3) & 
     &                    - rsolid3(ii) * dipj_d_r  * gfunc(4) & 
     &                    + rsolid3(ii) * quadj_dd_rr * gfunc(5) & 
     &                    - rsolid3(ii) * octj_ddd_rrr * gfunc(6)
                     fieldtens3(ii,jatom,jmol) = fieldtens3(ii,jatom,jmol) & 
     &                    +(- charge_i * rsolid3(ii) & 
     &                    + 6.0d0 * spherical3i(ii)) * gfunc(3) & 
     &                    - rsolid3(ii) * dipi_d_r  * gfunc(4) & 
     &                    - rsolid3(ii) * quadi_dd_rr * gfunc(5) & 
     &                    - rsolid3(ii) * octi_ddd_rrr * gfunc(6)
                  end do 
                  
                  do ii = 1,5
                     do jj = 1,3
                        fieldtens2_1(ii,jj,iatom,imol) = fieldtens2_1(ii,jj,iatom,imol) & 
     &                       + (3.0d0 * rsolid2(ii) * dipj(jj) &
     &                       + 6.0d0 * dif(jj) * spherical2j(ii)) * gfunc(3) & 
     &                       - (6.0d0 * rsolid2(ii) * quadj_d_r(jj) & 
     &                       + 18.0d0 * dif(jj) * octj_d_r(ii)) * gfunc(4) & 
     &                       + 9.0d0 * rsolid2(ii) * octj_dd_rr(jj) * gfunc(5)
                        
                        fieldtens2_1(ii,jj,jatom,jmol) = fieldtens2_1(ii,jj,jatom,jmol) & 
     &                       + (3.0d0 * rsolid2(ii) * dipi(jj) & 
     &                       - 6.0d0 * dif(jj) * spherical2i(ii)) * gfunc(3) & 
     &                       + (6.0d0 * rsolid2(ii) * quadi_d_r(jj) & 
     &                       - 18.0d0 * dif(jj) * octi_d_r(ii)) * gfunc(4) & 
     &                       + 9.0d0 * rsolid2(ii) * octi_dd_rr(jj) * gfunc(5)
                     end do 
                  end do 
                  
                  ggg(3) = ggg(3) +  charge_j * octi_ddd_rrr - charge_i * octj_ddd_rrr  & 
     &                 + 3.0d0 * (dipi_d_octj_dd_rr + dipj_d_octi_dd_rr) & 
     &                 + 6.0d0 * (quadj_dd_octi_d_r - quadi_dd_octj_d_r  & 
     &                 +  octi_ddd_octj)
                  ggg(4) = ggg(4) - dioctj_dddd_rrrr - djocti_dddd_rrrr & 
     &                 + 6.0d0 * (quadi_d_r_d_octj_dd_rr - quadj_d_r_d_octi_dd_rr) & 
     &                 - 18.0d0 * (octi_d_r_dd_octj_d_r)
                  ggg(5) = ggg(5) + octi_ddd_rrr * quadj_dd_rr - octj_ddd_rrr * quadi_dd_rr & 
     &                 + 9.0d0 * (octi_dd_rr_d_octj_dd_rr)
                  ggg(6) = ggg(6) - octi_ddd_rrr * octj_ddd_rrr
                  
                  
                  do ivec = 1,3
                     force_elec(ivec) = force_elec(ivec) &
     &                    + (3.0d0 * ( charge_j * octi_dd_rr(ivec) - charge_i * octj_dd_rr(ivec)) &
     &                    + 6.0d0 * (octj_d_r_d_dipj(ivec) + octi_d_r_d_dipj(ivec)) & 
     &                    + 6.0d0 * (quadj_dd_octi(ivec) - quadi_dd_octj(ivec)))* gfunc(3) &
                     
     &                    + (- (dipi(ivec) * octj_ddd_rrr + 3.0d0 * dipi_d_r * octj_dd_rr(ivec)) & 
     &                    - (dipj(ivec) * octi_ddd_rrr + 3.0d0 * dipj_d_r * octi_dd_rr(ivec)) &
     &                    + 6.0d0 * (quadi_d_octj_dd_rr(ivec) - quadj_d_octi_dd_rr(ivec)) & 
     &                    + 12.0d0 * (quadi_d_r_d_octj_d_r(ivec) - quadj_d_r_d_octi_d_r(ivec))  & 
     &                    - 18.0d0 * (octi_dd_octj_d_r(ivec) + octj_dd_octi_d_r(ivec)))*gfunc(4) &
                     
     &                    + (3.0d0 * (octi_dd_rr(ivec)*quadj_dd_rr - octj_dd_rr(ivec)*quadi_dd_rr) &
     &                    + 2.0d0 * (octi_ddd_rrr * quadj_d_r(ivec) - octj_ddd_rrr * quadi_d_r(ivec)) & 
     &                    + 18.0d0 * (octi_d_r_d_octj_dd_rr(ivec) + octj_d_r_d_octi_dd_rr(ivec))) * gfunc(5) & 
                     
     &                    - 3.0d0*(octi_dd_rr(ivec) * octj_ddd_rrr + octj_dd_rr(ivec) * octi_ddd_rrr)*gfunc(6) 
                  end do 
               endif
               
               do n = 0,irank + jrank
                  do ivec = 1,3
                     force_elec(ivec) = force_elec(ivec) - ggg(n) * gfunc(n+1) * dif(ivec) 
                  end do 
               end do 
               
               xforce = force_elec(1) * fac(2)
               yforce = force_elec(2) * fac(2)
               zforce = force_elec(3) * fac(2)
               
               if(.not.same_atom) then 
                  
                  ffix = ffix + xforce + xforce_pair
                  ffiy = ffiy + yforce + yforce_pair
                  ffiz = ffiz + zforce + zforce_pair
                  
                  ffjx = ffjx - xforce - xforce_pair
                  ffjy = ffjy - yforce - yforce_pair
                  ffjz = ffjz - zforce - zforce_pair
               endif
               
               if(periodic) then                        
                  
                  f1dif = cveci(1,1) * xdif + cveci(1,2) * ydif + cveci(1,3) * zdif
                  f2dif = cveci(2,2) * ydif + cveci(2,3) * zdif
                  f3dif = cveci(3,3) * zdif
                  
                  dudcvec_coulomb(1,1) = dudcvec_coulomb(1,1) + xforce * f1dif
                  
                  dudcvec_coulomb(1,2) = dudcvec_coulomb(1,2) + xforce * f2dif
                  dudcvec_coulomb(2,2) = dudcvec_coulomb(2,2) + yforce * f2dif
                  
                  dudcvec_coulomb(1,3) = dudcvec_coulomb(1,3) + xforce * f3dif
                  dudcvec_coulomb(2,3) = dudcvec_coulomb(2,3) + yforce * f3dif
                  dudcvec_coulomb(3,3) = dudcvec_coulomb(3,3) + zforce * f3dif
                  
                  dudcvec_pair(1,1) = dudcvec_pair(1,1) + xforce_pair * f1dif
                  
                  dudcvec_pair(1,2) = dudcvec_pair(1,2) + xforce_pair * f2dif
                  dudcvec_pair(2,2) = dudcvec_pair(2,2) + yforce_pair * f2dif
                  
                  dudcvec_pair(1,3) = dudcvec_pair(1,3) + xforce_pair * f3dif
                  dudcvec_pair(2,3) = dudcvec_pair(2,3) + yforce_pair * f3dif
                  dudcvec_pair(3,3) = dudcvec_pair(3,3) + zforce_pair * f3dif
                  
                  xxdif = xxdif0 + dcellx - xdif
                  yydif = yydif0 + dcelly - ydif
                  zzdif = zzdif0 + dcellz - zdif
                  
                  call nearest(xxdif,yydif,zzdif,scvec,scveci)
                  
                  f1 = cveci(1,1) * xxdif + cveci(1,2) * yydif + cveci(1,3) * zzdif
                  f2 = cveci(2,2) * yydif + cveci(2,3) * zzdif
                  f3 = cveci(3,3) * zzdif
                  
                  dudcvec_intra(1,1) = dudcvec_intra(1,1) + (xforce + xforce_pair) * f1
                  
                  dudcvec_intra(1,2) = dudcvec_intra(1,2) + (xforce + xforce_pair) * f2
                  dudcvec_intra(2,2) = dudcvec_intra(2,2) + (yforce + yforce_pair) * f2
                  
                  dudcvec_intra(1,3) = dudcvec_intra(1,3) + (xforce + xforce_pair) * f3
                  dudcvec_intra(2,3) = dudcvec_intra(2,3) + (yforce + yforce_pair) * f3
                  dudcvec_intra(3,3) = dudcvec_intra(3,3) + (zforce + zforce_pair) * f3
                  
               endif
               
            endif
         endif
         
      end do 
      end do 
      end do 
      
      force(1,iatom,imol) = force(1,iatom,imol) + ffix
      force(2,iatom,imol) = force(2,iatom,imol) + ffiy
      force(3,iatom,imol) = force(3,iatom,imol) + ffiz
                                               
      force(1,jatom,jmol) = force(1,jatom,jmol) + ffjx
      force(2,jatom,jmol) = force(2,jatom,jmol) + ffjy
      force(3,jatom,jmol) = force(3,jatom,jmol) + ffjz

      end do jatom_loop
      end do jmol_loop
      end do iatom_loop
      end do imol_loop

      end subroutine real


      subroutine self()
!----------------------------------------------
!     CALCULATE THE SELF TERMS
!----------------------------------------------
      use common_data
      implicit none
      real(8) ::  a,b,c,d,e,f
      real(8) :: trace_quad
      integer natoms_in_mol,imol,iatom,imoltype
      integer ivec,jvec,kvec,ii

      a = 4.0d0 * pi / (4.0d0 * pi**1.5)

      do imol = 1,nmol
         imoltype = mol_type(imol)
         natoms_in_mol = mol_natoms(imoltype)
         do iatom = 1,natoms_in_mol

            c = 2.0d0 * ewald_eps * a
            fieldtens0(iatom,imol) = fieldtens0(iatom,imol) - c * lab_charge(iatom,imol)

            if(max_rank.ge.1) then 
               d = 2.0d0 * ((2.0d0/3.0d0) * ewald_eps **3 ) * a
               do ii = 1,3
                  fieldtens1(ii,iatom,imol) = fieldtens1(ii,iatom,imol) - d * lab_dipole(ii,iatom,imol)
               end do 
            endif

            if(max_rank.ge.2) then 
               e = (16.0d0/5.0d0) * ewald_eps**5 * a
               do ii = 1,5
                  fieldtens2(ii,iatom,imol) = fieldtens2(ii,iatom,imol) - e * lab_spherical2(ii,iatom,imol)
               end do 
            endif

            if(max_rank.ge.3) then 
               f = (96.0d0 / 7.0d0) * ewald_eps**7 * a 
               do ii = 1,7
                  fieldtens3(ii,iatom,imol) = fieldtens3(ii,iatom,imol) - f * lab_spherical3(ii,iatom,imol)
               end do 
            endif
         end do 
      end do 
      
      end subroutine self


      subroutine get_rigid_frac_forces()
!     calculates the rigid body fractional forces, assuming the rigid body forces have already been calculated
      use common_data
      use nr_mod
      implicit none
      integer imol,ivec
      real(8), dimension(3) :: vec1,vec2
      real(8), dimension(3,3) :: cvec_trans

      cvec_trans = transpose(cvec)

      do imol = 1,nmol
         vec1(1) = forcemol(1,imol)
         vec1(2) = forcemol(2,imol)
         vec1(3) = forcemol(3,imol)

         call mat3multvec(cvec_trans,vec1,vec2)

         fracforcemol(1,imol) = vec2(1)
         fracforcemol(2,imol) = vec2(2)
         fracforcemol(3,imol) = vec2(3)
      end do 


      end subroutine get_rigid_frac_forces
     

      subroutine calculate_rigid_body_forces()
!     calculates the rigid body forces from the Cartesian forces
      use common_data
      implicit none
      integer imol,iatom,natoms_in_mol,moltype
      real(8) :: fx,fy,fz,tx,ty,tz
      real(8), dimension(3,3,3) :: drotmat
      real(8), dimension(3,3) :: rotmat
      real(8), dimension(3) :: r0,r1
      real(8), dimension(3,3) :: drot
      real(8), dimension(3) :: dudeuler
      integer :: ivec,jvec,kvec,lvec,mvec,nvec,m,j,ii,jj,conformer
      real(8) :: phi,theta,psi
      real(8), dimension(3) :: d0,d1
      real(8) :: dfx,dfy,dfz
      real(8) :: q_field,o_field
      real(8), dimension(3,3) :: quad0
      real(8), dimension(3) :: qvec,fvec,ovec

      do imol = 1,nmol
         moltype = mol_type(imol)
         natoms_in_mol = mol_natoms(moltype)
         conformer = mol_conformer(imol)

         fx = 0.0d0 
         fy = 0.0d0 
         fz = 0.0d0 
         
         tx = 0.0d0 
         ty = 0.0d0 
         tz = 0.0d0 

         dudeuler = 0.0d0 

         phi = mol_euler(1,imol)
         theta = mol_euler(2,imol)
         psi = mol_euler(3,imol)

         call get_rotmat(phi,theta,psi,rotmat)
         call get_drotmat(phi,theta,psi,drotmat)

         do iatom = 1,natoms_in_mol
            
            fx = fx + force(1,iatom,imol)
            fy = fy + force(2,iatom,imol)
            fz = fz + force(3,iatom,imol)

            tx = tx + lab_relative_coord(2,iatom,imol)*force(3,iatom,imol) & 
     &              - lab_relative_coord(3,iatom,imol)*force(2,iatom,imol)
            ty = ty + lab_relative_coord(3,iatom,imol)*force(1,iatom,imol) & 
     &              - lab_relative_coord(1,iatom,imol)*force(3,iatom,imol)
            tz = tz + lab_relative_coord(1,iatom,imol)*force(2,iatom,imol) & 
     &              - lab_relative_coord(2,iatom,imol)*force(1,iatom,imol)


            if(max_rank.ge.1) then 
!     dipole contribution
               do jj = 1,3
                  qvec(jj) = lab_dipole(jj,iatom,imol)
                  fvec(jj) = fieldtens1(jj,iatom,imol) * fac(2)
               end do 

               tx = tx - 1.0d0 * (qvec(2) * fvec(3) - qvec(3) * fvec(2))
               ty = ty - 1.0d0 * (qvec(3) * fvec(1) - qvec(1) * fvec(3))
               tz = tz - 1.0d0 * (qvec(1) * fvec(2) - qvec(2) * fvec(1))
            endif
            if(max_rank.ge.2) then 
!     quadrupole contribution
               do ii = 1,3
                  do jj = 1,3
                     qvec(jj) = lab_quad(ii,jj,iatom,imol)
                     fvec(jj) = fieldtens1_1(ii,jj,iatom,imol) * fac(2)
                  end do 

                  tx = tx - 2.0d0 * (qvec(2) * fvec(3) - qvec(3) * fvec(2))
                  ty = ty - 2.0d0 * (qvec(3) * fvec(1) - qvec(1) * fvec(3))
                  tz = tz - 2.0d0 * (qvec(1) * fvec(2) - qvec(2) * fvec(1))

               end do 
            endif
            if(max_rank.ge.3) then 
!     octopole contribution

               do ii = 1,5
                  do jj = 1,3
                     ovec(jj) = lab_spherical2_1(ii,jj,iatom,imol)
                     fvec(jj) = fieldtens2_1(ii,jj,iatom,imol) * fac(2)
                  end do 
                  
                  tx = tx - 3.0d0 * (ovec(2) * fvec(3) - ovec(3) * fvec(2))
                  ty = ty - 3.0d0 * (ovec(3) * fvec(1) - ovec(1) * fvec(3))
                  tz = tz - 3.0d0 * (ovec(1) * fvec(2) - ovec(2) * fvec(1))
                  
               end do 
            endif

            r0(1) = site_coord0(1,iatom,conformer,moltype)
            r0(2) = site_coord0(2,iatom,conformer,moltype)
            r0(3) = site_coord0(3,iatom,conformer,moltype)

            do m = 1,3
               do ivec = 1,3
                  do jvec = 1,3
                     drot(ivec,jvec) = drotmat(ivec,jvec,m)
                  end do 
               end do 

               r1 = matmul(drot,r0) 
               
               dudeuler(m) = dudeuler(m) + force(1,iatom,imol) * r1(1) + &
     &                                     force(2,iatom,imol) * r1(2) + &
     &                                     force(3,iatom,imol) * r1(3)

               if(max_rank.ge.1) then 
                  do ivec = 1,3
                     d0(ivec) = site_dipole0(ivec,iatom,conformer,moltype)
                  end do 
                  
                  d1 = matmul(drot,d0)

                  do ivec = 1,3
                     dudeuler(m) = dudeuler(m) - fac(2) * d1(ivec) * dfield(ivec,iatom,imol)
                  end do 
               endif

               if(max_rank.ge.2) then 
                  
                  do ivec = 1,3
                     do jvec = 1,3
                        q_field = qfield(ivec,jvec,iatom,imol)
                        do kvec = 1,3
                           do lvec = 1,3
                              dudeuler(m) = dudeuler(m)  & 
     &                             - 2.0d0 * drot(ivec,kvec) * rotmat(jvec,lvec) &
     &                             * q_field * fac(2) * site_quad0(kvec,lvec,iatom,conformer,moltype) 
                           end do 
                        end do 
                     end do 
                  end do 
               endif

               if(max_rank.ge.3) then 

                  do ivec = 1,3
                     do jvec = 1,3
                        do kvec = 1,3
                           o_field = ofield(ivec,jvec,kvec,iatom,imol)
                           do lvec = 1,3
                              do mvec = 1,3
                                 do nvec = 1,3
                                    dudeuler(m) = dudeuler(m)  & 
     &                             - 3.0d0 * drot(ivec,lvec) * rotmat(jvec,mvec) * rotmat(kvec,nvec) & 
     &                             * o_field * fac(2) * site_oct0(lvec,mvec,nvec,iatom,conformer,moltype) 
                                 end do 
                              end do 
                           end do 
                        end do 
                     end do 
                  end do 
               endif

            end do !m loop
         end do !iatom loop


         forcemol(1,imol) = forcemol(1,imol) + fx
         forcemol(2,imol) = forcemol(2,imol) + fy
         forcemol(3,imol) = forcemol(3,imol) + fz

         torquemol(1,imol) = torquemol(1,imol) + tx
         torquemol(2,imol) = torquemol(2,imol) + ty
         torquemol(3,imol) = torquemol(3,imol) + tz

         dudeulermol(1,imol) = dudeuler(1)
         dudeulermol(2,imol) = dudeuler(2)
         dudeulermol(3,imol) = dudeuler(3)

      end do 
      end subroutine calculate_rigid_body_forces

      subroutine print_rigid_coordinates(ifile)
!     writes the rigid body coordinates to file
      use common_data
      use nr_mod
      implicit none
      integer :: imol,moltype,imode
      integer :: ifile
      real(8) :: x,y,z
      real(8) :: amag,bmag,cmag,alpha,beta,gamma
      real(8) :: f1,f2,f3,g1,g2,g3
      real(8), dimension(3) :: vec1,vec2

      if(.not.periodic) rigid_output_mode = 'XYZ'

      write(ifile,*) nmol
      if(periodic) then
         if(write_cvec.eq.'XYZ') then 
            write(ifile,*) 'XYZ',cvec(1,1),cvec(1,2),cvec(2,2),cvec(1,3),cvec(2,3),cvec(3,3) & 
     &           ,utotal * pefac,pressure_internal*fac(7),density*fac(8)
         else if(write_cvec.eq.'ANGLE') then
            call get_angle_cvec(cvec,amag,bmag,cmag,alpha,beta,gamma)
            write(ifile,*) 'ANGLE',amag,bmag,cmag,alpha,beta,gamma & 
     &           ,utotal * pefac,pressure_internal*fac(7),density*fac(8)
         endif
      else
         write(ifile,*) 
      endif
      write(ifile,*) 'RIGID ',rigid_output_mode
      do imol = 1,nmol
         moltype = mol_type(imol)
         if(periodic) then 
            f1 = rfrac(1,imol)
            f2 = rfrac(2,imol)
            f3 = rfrac(3,imol)
            f1 = f1 - anint(f1)
            f2 = f2 - anint(f2)
            f3 = f3 - anint(f3)
            if(rigid_output_mode.eq.'FRAC') then
               g1 = f1
               g2 = f2
               g3 = f3
            else if(rigid_output_mode.eq.'XYZ') then
               vec1(1) = f1
               vec1(2) = f2
               vec1(3) = f3
               call mat3multvec(cvec,vec1,vec2)               
               x = vec2(1)
               y = vec2(2)
               z = vec2(3)
               g1 = x
               g2 = y
               g3 = z
            endif
         else
            x = mol_com(1,imol)
            y = mol_com(2,imol)
            z = mol_com(3,imol)
            g1 = x
            g2 = y
            g3 = z
         endif
         write(ifile,17) molecule_name_list(moltype),mol_conformer(imol),g1,g2,g3 &
     & ,mol_euler(1,imol),mol_euler(2,imol),mol_euler(3,imol)

 17      format('MOLECULE  ',a16,' CONFORMER',i4,' COORD',6f12.6)

      end do 
      call flush(ifile)
      
      end subroutine print_rigid_coordinates

      subroutine print_coordinates(ifile)
!     writes the cartesian coordinates to file

      use common_data
      implicit none
      integer imol,iatom
      real(8) :: x,y,z
      character(len = 3) :: atomic_name
      real(8) :: amag,bmag,cmag,alpha,beta,gamma
      real(8) :: s,dx,dy,dz,dmag
      integer natoms_in_mol,moltype,atno
      integer ifile
      logical okflag

      write(ifile,*) nmassatoms
      if(periodic) then
         if(write_cvec.eq.'XYZ') then 
            write(ifile,*) 'XYZ',cvec(1,1),cvec(1,2),cvec(2,2),cvec(1,3),cvec(2,3),cvec(3,3) & 
     &           ,utotal * pefac,pressure_internal*fac(7),density*fac(8)
         else if(write_cvec.eq.'ANGLE') then
            call get_angle_cvec(cvec,amag,bmag,cmag,alpha,beta,gamma)
            write(ifile,*) 'ANGLE',amag,bmag,cmag,alpha,beta,gamma & 
     &           ,utotal * pefac,pressure_internal*fac(7),density*fac(8)
         endif
         else
            write(ifile,*) utotal * pefac
         endif
      do imol = 1,nmol

         moltype = mol_type(imol)
         natoms_in_mol = mol_natoms(moltype)
         do iatom = 1,natoms_in_mol
            x = lab_coord(1,iatom,imol)
            y = lab_coord(2,iatom,imol)
            z = lab_coord(3,iatom,imol)
            atomic_name = site_name(moltype,iatom)
            if(atomic_name.ne.'X') then 
               write(ifile,*) atomic_name,x,y,z
            endif
         end do 
      end do 

      call flush(ifile)
      end subroutine print_coordinates

      subroutine generate_cartesians()
!     generates the cartesians from the XYZ and the euler angles
      use common_data
      implicit none
      integer natoms_in_mol,moltype,iatom,imol
      real(8) :: phi,theta,psi
      real(8), dimension(3,3) :: rotmat
      real(8), dimension(3) :: r0,r1,d0,d1
      real(8) :: xcom,ycom,zcom
      real(8) :: x,y,z
      real(8), dimension(3,3) :: quad0,quad1
      real(8), dimension(3,3,3) :: oct0,oct1
      real(8) :: qxx,qxy,qxz,qyy,qyz,qzz
      real(8) :: qxxx,qxxy,qxyy,qyyy,qxxz,qxyz,qyyz,qxzz,qyzz,qzzz
      integer dim,imode,ivec,jvec,kvec,ii,jj
      integer conformer

      do imol = 1,nmol
         moltype = mol_type(imol)
         natoms_in_mol = mol_natoms(moltype)
         conformer = mol_conformer(imol)
         dim = dimensionality(moltype)

         xcom = mol_com(1,imol)
         ycom = mol_com(2,imol)
         zcom = mol_com(3,imol)

         if(dim.eq.1) then

            lab_coord(1,1,imol) = xcom
            lab_coord(2,1,imol) = ycom
            lab_coord(3,1,imol) = zcom
            
            lab_relative_coord(1,1,imol) = 0.0d0
            lab_relative_coord(2,1,imol) = 0.0d0
            lab_relative_coord(3,1,imol) = 0.0d0
            
         else
            phi = mol_euler(1,imol)
            theta = mol_euler(2,imol)
            psi = mol_euler(3,imol)
            
            do iatom = 1,natoms_in_mol
               do ivec = 1,3
                  r0(ivec) = site_coord0(ivec,iatom,conformer,moltype)
               end do 
               
               call get_rotmat(phi,theta,psi,rotmat)
               r1 = matmul(rotmat,r0)
               
               lab_coord(1,iatom,imol) = r1(1) + xcom
               lab_coord(2,iatom,imol) = r1(2) + ycom
               lab_coord(3,iatom,imol) = r1(3) + zcom
            
               lab_relative_coord(1,iatom,imol) = r1(1)
               lab_relative_coord(2,iatom,imol) = r1(2)
               lab_relative_coord(3,iatom,imol) = r1(3)

               lab_charge(iatom,imol) = site_charge0(iatom,conformer,moltype)

               if(max_rank.ge.1) then 
                  do ivec = 1,3
                     d0(ivec) = site_dipole0(ivec,iatom,conformer,moltype)
                  end do 
                  call get_rotmat(phi,theta,psi,rotmat)
                  d1 = matmul(rotmat,d0)
                  do ivec = 1,3
                     lab_dipole(ivec,iatom,imol)  = d1(ivec)
                  end do 
               endif
               if(max_rank.ge.2) then 
                  do ivec = 1,3
                     do jvec = 1,3
                        quad0(ivec,jvec) = site_quad0(ivec,jvec,iatom,conformer,moltype)
                     end do 
                  end do 
                  call rotate_2tensor(rotmat,quad0,quad1)
                  do ivec = 1,3
                     do jvec = 1,3
                        lab_quad(ivec,jvec,iatom,imol) = quad1(ivec,jvec)
                     end do 
                  end do 

!     find the quadrupole in spherical harmonics

                  qxx = quad1(1,1)
                  qxy = quad1(1,2)
                  qxz = quad1(1,3)
                  qyy = quad1(2,2)
                  qyz = quad1(2,3)
                  qzz = quad1(3,3)

                  call get_lab_spherical2(iatom,imol,qxx,qxy,qxz,qyy,qyz,qzz)
                  
               endif
               if(max_rank.ge.3) then 
                  do ivec = 1,3
                     do jvec = 1,3
                        do kvec = 1,3
                           oct0(ivec,jvec,kvec) = site_oct0(ivec,jvec,kvec,iatom,conformer,moltype)
                        end do 
                     end do 
                  end do 
                  call rotate_3tensor(rotmat,oct0,oct1)
                  do ivec = 1,3
                     do jvec = 1,3
                        do kvec = 1,3
                           lab_oct(ivec,jvec,kvec,iatom,imol) = oct1(ivec,jvec,kvec)
                        end do 
                     end do 
                  end do 


!     find the octopole in spherical harmonics

                  qxxx = oct1(1,1,1)
                  qxxy = oct1(1,1,2)
                  qxyy = oct1(1,2,2)
                  qyyy = oct1(2,2,2)
                  qxxz = oct1(1,1,3)
                  qxyz = oct1(1,2,3)
                  qyyz = oct1(2,2,3)
                  qxzz = oct1(1,3,3)
                  qyzz = oct1(2,3,3)
                  qzzz = oct1(3,3,3)

                  call get_lab_spherical3(iatom,imol,qxxx,qxxy,qxyy,qyyy,qxxz,qxyz,qyyz,qxzz,qyzz,qzzz)

               endif
            
            end do 
         endif
         
      end do 

      end subroutine generate_cartesians

      subroutine get_rotmat(phi,theta,psi,rotmat)
!     generates the rotation matrix from the euler angles

      use common_data
      implicit none
      real(8), dimension(3,3) :: rotmat
      real(8) :: phi,theta,psi
      real(8) :: s_phi,c_phi,s_theta,c_theta,s_psi,c_psi

      s_phi = dsin(phi)
      c_phi = dcos(phi)

      s_theta = dsin(theta)
      c_theta = dcos(theta)

      s_psi = dsin(psi)
      c_psi = dcos(psi)

      rotmat(1,1) = c_phi * c_psi - s_phi * c_theta * s_psi
      rotmat(1,2) = s_phi * c_psi + c_phi * c_theta * s_psi
      rotmat(1,3) = s_theta * s_psi
      
      rotmat(2,1) = -c_phi * s_psi - s_phi * c_theta * c_psi
      rotmat(2,2) = -s_phi * s_psi + c_phi * c_theta * c_psi
      rotmat(2,3) = s_theta * c_psi

      rotmat(3,1) = s_phi * s_theta
      rotmat(3,2) = -c_phi * s_theta
      rotmat(3,3) = c_theta

      end subroutine get_rotmat


      subroutine get_drotmat(phi,theta,psi,drotmat)
!     generates the derivatives of the rotation matrix from the euler angles

      use common_data
      implicit none
      real(8), dimension(3,3,3) :: drotmat
      real(8) :: phi,theta,psi
      real(8) :: s_phi,c_phi,s_theta,c_theta,s_psi,c_psi

      s_phi = dsin(phi)
      c_phi = dcos(phi)

      s_theta = dsin(theta)
      c_theta = dcos(theta)

      s_psi = dsin(psi)
      c_psi = dcos(psi)
      

!     derivatives wrt phi

      drotmat(1,1,1) = -s_phi * c_psi - c_phi * c_theta * s_psi
      drotmat(1,2,1) = c_phi * c_psi - s_phi * c_theta * s_psi
      drotmat(1,3,1) = 0.0d0
      
      drotmat(2,1,1) = s_phi * s_psi - c_phi * c_theta * c_psi   
      drotmat(2,2,1) = -c_phi * s_psi - s_phi * c_theta *c_psi    
      drotmat(2,3,1) = 0.0d0 

      drotmat(3,1,1) = c_phi * s_theta
      drotmat(3,2,1) = s_phi * s_theta
      drotmat(3,3,1) = 0.0d0 

!     derivatives wrt theta

      drotmat(1,1,2) = s_phi * s_theta * s_psi
      drotmat(1,2,2) = - c_phi * s_theta * s_psi
      drotmat(1,3,2) = c_theta * s_psi
      
      drotmat(2,1,2) = s_phi * s_theta * c_psi
      drotmat(2,2,2) = - c_phi * s_theta *c_psi
      drotmat(2,3,2) = c_theta * c_psi

      drotmat(3,1,2) = s_phi * c_theta
      drotmat(3,2,2) = -c_phi * c_theta
      drotmat(3,3,2) = -s_theta

!     derivatives wrt psi

      drotmat(1,1,3) = - c_phi * s_psi - s_phi * c_theta * c_psi
      drotmat(1,2,3) = - s_phi * s_psi + c_phi * c_theta * c_psi
      drotmat(1,3,3) = s_theta * c_psi
      
      drotmat(2,1,3) = -c_phi * c_psi + s_phi * c_theta * s_psi
      drotmat(2,2,3) = -s_phi * c_psi - c_phi * c_theta * s_psi
      drotmat(2,3,3) = -s_theta * s_psi

      drotmat(3,1,3) = 0.0d0 
      drotmat(3,2,3) = 0.0d0 
      drotmat(3,3,3) = 0.0d0 

      end subroutine get_drotmat

      subroutine get_rotmat_for_axis(uvec,theta,axismat)
      implicit none
      real(8) :: theta,c,s
      real(8), dimension(3,3) :: axismat
      real(8) :: uu,u2
      real(8), dimension(3) :: uvec
      integer :: ivec,jvec

!     normalise 

      u2 = uvec(1)**2 + uvec(2)**2 + uvec(3)**2
      uu = dsqrt(u2)
      uvec = uvec / uu 
      
      c = dcos(theta)
      s = dsin(theta)

      do ivec = 1,3
         do jvec = 1,3
            axismat(ivec,jvec) = (1.0d0 - c) * uvec(ivec)*uvec(jvec)
            if(ivec.eq.jvec) axismat(ivec,jvec) = axismat(ivec,jvec) + c
         end do 
      end do 
      
      axismat(1,2) = axismat(1,2) - uvec(3) * s
      axismat(2,1) = axismat(2,1) + uvec(3) * s

      axismat(1,3) = axismat(1,3) + uvec(2) * s
      axismat(3,1) = axismat(3,1) - uvec(2) * s

      axismat(2,3) = axismat(2,3) - uvec(1) * s
      axismat(3,2) = axismat(3,2) + uvec(1) * s

      end subroutine get_rotmat_for_axis

      subroutine read_control(control_list,param_list,ncontrol)
!     reads the control file

      use common_data
      use parse_text
      implicit none
      character(len = 32), dimension(1024) :: control_list
      real(8), dimension(1024,8) :: param_list
      integer ncontrol

      integer ios
      integer atno,site_index,moltype,k
      real(8) :: press,temp,schance
      integer imol
      character(len=2048) :: buffer
      character(len = 2048), dimension(64) :: textlist
      integer line
      real(8) :: radius,ewald_error,rmin,rdamp,delta
      real(8) :: delta_frac,delta_xyz,delta_cvec,delta_angle,dmass
      real(8) :: tolerance
      character(len=10) :: totest,toprint
      character(len=8) :: cvec_mode
      character(len=32) :: filename
      integer imode,nitems
      logical file_exists

!     read the control file

!     default values
      temp = 100.0d0 
      input_energy_unit = 'KJ'
      pefac = fac(3)

      ios = 0
      line = 0 
      ncontrol = 0 
      control_list = ""
      param_list = 0
      do while (ios == 0)
         read(40, '(A)', iostat=ios) buffer

         if (ios == 0) then
            line = line + 1

            call split_text_to_list(buffer,textlist,nitems)

            select case(upcase(textlist(1)))
            case('INPUT_FILE')
               read(textlist(2),*,iostat=ios) filename
               inquire(file = filename,exist = file_exists)
               if(.not.file_exists) call file_doesnt_exist(filename)
               open(10,file = filename)
            case('MODEL_FILE')
               read(textlist(2),*,iostat=ios) filename
               inquire(file = filename,exist = file_exists)
               if(.not.file_exists) call file_doesnt_exist(filename)
               open(30,file = filename)
            case('GEOMETRY_FILE')
               read(textlist(2),*,iostat=ios) filename
               inquire(file = filename,exist = file_exists)
               if(.not.file_exists) call file_doesnt_exist(filename)
               open(20,file = filename)
            case('MULTIPOLE_FILE')
               read(textlist(2),*,iostat=ios) filename
               inquire(file = filename,exist = file_exists)
               if(.not.file_exists) call file_doesnt_exist(filename)
               open(29,file = filename)
            case('BATCH_FILE')
               read(textlist(2),*,iostat=ios) filename
               inquire(file = filename,exist = file_exists)
               if(.not.file_exists) call file_doesnt_exist(filename)
               batch_file = filename
            case('PROB_FILE')
               read(textlist(2),*,iostat=ios) filename
               inquire(file = filename,exist = file_exists)
               if(.not.file_exists) call file_doesnt_exist(filename)
               open(43,file = filename)
            case('PRESSURE')
               read(textlist(2),*,iostat=ios) press
               ncontrol = ncontrol + 1
               control_list(ncontrol) = "PRESSURE"
               param_list(ncontrol,1) = press
            case('MINTOL')
               read(textlist(2),*,iostat=ios) tolerance
               ncontrol = ncontrol + 1
               control_list(ncontrol) = 'TOLERANCE'
               param_list(ncontrol,1) = tolerance
            case('EWALD_ERROR')
              read(textlist(2),*,iostat=ios) ewald_error
              ncontrol = ncontrol + 1
              control_list(ncontrol) = "EWALD_ERROR"
              param_list(ncontrol,1) = ewald_error
           case('BATCH')
              ncontrol = ncontrol + 1
              control_list(ncontrol) = 'BATCH'
           case('BATCH_XYZ')
              ncontrol = ncontrol + 1
              control_list(ncontrol) = 'BATCH_XYZ'
            case('RIGID')
               ncontrol = ncontrol + 1
               if(textlist(2).eq.'OUTPUT') then
                  if(textlist(3).eq.'FRAC') then
                     control_list(ncontrol) = 'RIGID OUTPUT FRAC'
                  else if(textlist(3).eq.'XYZ') then
                     control_list(ncontrol) = 'RIGID OUTPUT XYZ'
                  endif
               endif
            case('BASIN')
               if(textlist(2).eq.'TEMP') then
                  read(textlist(3),*,iostat = ios) temp
                  ncontrol = ncontrol + 1
                  param_list(ncontrol,1) = temp
                  control_list(ncontrol) = 'BASIN TEMP'
               else if(textlist(2).eq.'RUN') then
                  ncontrol = ncontrol + 1
                  control_list(ncontrol) = 'BASIN RUN'
               else if(textlist(2).eq.'RESTART') then
                  ncontrol = ncontrol + 1
                  control_list(ncontrol) = 'BASIN RESTART'
               else if(textlist(2).eq.'DELTA_FRAC') then
                  ncontrol = ncontrol + 1
                  control_list(ncontrol) = 'DELTA_FRAC'
                  read(textlist(3),*,iostat = ios) delta_frac
                  param_list(ncontrol,1) = delta_frac
               else if(textlist(2).eq.'DELTA_XYZ') then
                  ncontrol = ncontrol + 1
                  control_list(ncontrol) = 'DELTA_XYZ'
                  read(textlist(3),*,iostat = ios) delta_xyz
                  param_list(ncontrol,1) = delta_xyz
               else if(textlist(2).eq.'DELTA_CVEC') then
                  ncontrol = ncontrol + 1
                  control_list(ncontrol) = 'DELTA_CVEC'
                  read(textlist(3),*,iostat = ios) delta_cvec
                  param_list(ncontrol,1) = delta_cvec
               else if(textlist(2).eq.'DELTA_ANGLE') then
                  ncontrol = ncontrol + 1
                  control_list(ncontrol) = 'DELTA_ANGLE'
                  read(textlist(3),*,iostat = ios) delta_angle
                  param_list(ncontrol,1) = delta_angle
               else if(textlist(2).eq.'RMIN') then
                  ncontrol = ncontrol + 1
                  control_list(ncontrol) = 'BASIN RMIN'
                  read(textlist(3),*,iostat=ios) rmin
                  param_list(ncontrol,1) = rmin
               else if(textlist(2).eq.'SWAP') then 
                  ncontrol = ncontrol + 1
                  control_list(ncontrol) = 'BASIN SWAP'
                  read(textlist(3),*,iostat = ios) schance
                  param_list(ncontrol,1) = schance
               else if(textlist(2).eq.'CONFORMER_SWAP') then 
                  ncontrol = ncontrol + 1
                  control_list(ncontrol) = 'BASIN CONFORMER_SWAP'
                  read(textlist(3),*,iostat = ios) schance
                  param_list(ncontrol,1) = schance
               endif
            case('FINDMIN')
               ncontrol = ncontrol + 1
               control_list(ncontrol) = 'FINDMIN'
            case('ENERGY')
               ncontrol = ncontrol + 1
               control_list(ncontrol) = "ENERGY"
            case('TEST')
               read(textlist(2),*,iostat=ios) totest
               ncontrol = ncontrol + 1
               control_list(ncontrol) = totest
               delta = 1.d-6
               if(nitems.gt.2) read(textlist(3),*,iostat = ios) delta
               param_list(ncontrol,1) = delta
            case('STOP')
               ncontrol = ncontrol + 1
               control_list(ncontrol) = 'STOP'
            case('RCUT')
               read(textlist(2),*,iostat=ios) radius
               ncontrol = ncontrol + 1
               control_list(ncontrol) = 'RCUT'
               param_list(ncontrol,1) = radius
            case('RDAMP_CUT')
               read(textlist(2),*,iostat = ios) radius
               ncontrol = ncontrol + 1
               control_list(ncontrol) = 'RDAMP_CUT'
               param_list(ncontrol,1) = radius
            case('INPUT')
               if(textlist(2).eq.'UNITS') then
                  if(textlist(3).eq.'KJ') then
                     input_energy_unit = 'KJ'
                     pefac = fac(3)
                     input_energy_unit = 'kJ/mol'
                  else if(textlist(3).eq.'KCAL') then
                     input_energy_unit = 'KCAL'
                     pefac = fac(4)
                     input_energy_unit = 'kcal/mol'
                  else
                     write(*,*) 'ERROR: INVALID INPUT UNITS TYPE'
                  endif
               endif
            case('OUTPUT')
               if(textlist(2).eq.'UNITS') then
                  ncontrol = ncontrol + 1
                  if(textlist(3).eq.'KJ') then
                     control_list(ncontrol) = 'OUTPUT KJ'
                  else if(textlist(3).eq.'KCAL') then
                     control_list(ncontrol) = 'OUTPUT KCAL'
                  else
                     write(*,*) 'ERROR: INVALID INPUT UNITS TYPE'
                  endif
               endif
            case('CVEC')
               if(textlist(2).eq.'INPUT') then

               read(textlist(3),*,iostat = ios) cvec_mode
               cvec_mode = upcase(cvec_mode)

               if(cvec_mode.eq.'XYZ') then
                  read_cvec = 'xyz'
               else if(cvec_mode.eq.'ANGLE') then 
                  read_cvec = 'angle'
               endif
               else if(textlist(2).eq.'OUTPUT') then
                  ncontrol = ncontrol + 1
                  read(textlist(3),*,iostat = ios) cvec_mode
                  if(cvec_mode.eq.'XYZ') then
                     control_list(ncontrol) = 'CVEC OUTPUT XYZ'
                  else if(cvec_mode.eq.'ANGLE') then
                     control_list(ncontrol) = 'CVEC OUTPUT ANGLE'                        
                  endif
               endif
            case('PRINT')
               read(textlist(2),*,iostat=ios) toprint
               ncontrol = ncontrol + 1
               if(toprint.eq.'RIGID') then
               control_list(ncontrol) = 'PRINT_RIGID'
               else if(toprint.eq.'XYZ') then
                  control_list(ncontrol) = 'PRINT_XYZ'
               else if(toprint.eq.'MULTIPOLES') then 
                  control_list(ncontrol) = 'PRINT_MULTIPOLES'
               else
                  write(*,*) 'INVALID PRINT OPTION'
                  stop
               endif
            case default
               write(*,21) line
               write(*,*)
 21            format('SKIPPING INVALID LABEL IN CONTROL FILE AT LINE', i4)
            end select
         end if
      end do

      if(input_energy_unit.eq.'KJ') then 
         write(*,*)
         write(*,*) 'INPUT ENERGY UNIT IS kJ/mol'
      else if(input_energy_unit.eq.'KCAL') then 
         write(*,*) 
         write(*,*) 'INPUT ENERGY UNIT IS kcal/mol'
      endif

      end subroutine read_control
      
      subroutine read_molecular_geometry
      use common_data
      use parse_text
      implicit none
      integer natoms_in_mol,moltype,ios,site_index,nmass_atoms_in_mol,iatom
      character(len = 16) molecule_name
      character(len = 3) atomic_name
      character(len=2048) :: buffer
      character(len = 32) :: token
      character(len = 2048), dimension(64) :: textlist
      real(8) :: x,y,z,dmass
      real(8) :: dx,dy,dz
      real(8) :: conformer_energy
      integer nitems,conformer,m,max_moltypes,nitem
      logical accept

      max_moltypes = 0
      conformer = 0 
      do while(.true.)
         read(20, '(A)', iostat=ios) buffer
         if(ios.ne.0) exit

         call split_text_to_list(buffer,textlist,nitems)

!     defaults
         conformer = 1
         conformer_energy = 0 
         
         nitem = 0 
         do while(.true.) 
            nitem = nitem + 1
            if(nitem.gt.nitems) exit

            read(textlist(nitem),*,iostat = ios) token

            if(ios.ne.0) exit
            token = upcase(token)

            select case(token)
            case('MOLECULE')
               nitem = nitem + 1
               read(textlist(nitem),*) molecule_name
               molecule_name = upcase(molecule_name)
            case('NATOMS')
               nitem = nitem + 1
               read(textlist(nitem),*) natoms_in_mol
            case('CONFORMER')
               nitem = nitem + 1
               read(textlist(nitem),*) conformer
            case('ENERGY')
               nitem = nitem + 1
               read(textlist(nitem),*) conformer_energy
            case default
               write(*,*) 'ERROR IN MOLECULAR GEOMETRY FILE: &
     &UNIDENTIFIED TOKEN ', token
               stop
            end select
         end do 

         if(ios.ne.0) exit

!     check to see if natoms_max is big enough

         if(natoms_in_mol.gt.natom_max) then 
            write(*,*) 'ERROR: NATOM_MAX IN SIZE.DAT NEEDS &
     &TO BE AT LEAST ',natoms_in_mol
            stop
         endif

!     check to see if the molecule is on the list
         accept = .false.
         do m = 1,max_moltypes
            if(molecule_name_list(m).eq.molecule_name) then 
               accept = .true.
               exit
            endif
         end do 
         if(accept) then 
            moltype = m
            mol_nconformers(moltype) = mol_nconformers(moltype) + 1
         else
            max_moltypes = max_moltypes + 1
            moltype = max_moltypes
            mol_nconformers(moltype) = 1
         endif

!     check array sizes in size.dat
         if(max_moltypes.gt.nmoltype_max) then
            write(*,*) 'ERROR: NMOL_MAX IN SIZE.DAT NEEDS TO BE &
   &AT LEAST ',max_moltypes
            stop
         endif
         if(mol_nconformers(moltype).gt.nconformer_max) then 
            write(*,*) 'ERROR: NCONFORMER_MAX IN SIZE.DAT NEEDS TO BE &
   &AT LEAST ',mol_nconformers(moltype)
            stop
         endif
         if(max_moltypes.gt.nmoltype_max) then 
            write(*,*) 'ERROR: NMOLTYPE_MAX IN SIZE.DAT NEEDS TO BE &
   &AT LEAST ',max_moltypes
            stop
         endif

         molecule_name_list(moltype) = molecule_name
         mol_natoms(moltype) = natoms_in_mol
         mol_conformer_energy(conformer,moltype) = conformer_energy / pefac

!     set the dimensionality of the molecule
         if(natoms_in_mol.eq.1) then
            dimensionality(moltype) = 0
         else
            dimensionality(moltype) = 3
         endif

!     It used to be that 1-atom molecules didn't need Euler angles, 
!     but, now we have multipoles, it's more complicated, 
!     so, for now, set the dimensionality of everything to 3
         dimensionality = 3


         site_index = 0 

         do iatom = 1,natoms_in_mol
            read(20,*,iostat = ios) atomic_name,x,y,z
            if(ios.ne.0) exit
            site_index = site_index + 1
            
            site_coord0(1,site_index,conformer,moltype) = x
            site_coord0(2,site_index,conformer,moltype) = y
            site_coord0(3,site_index,conformer,moltype) = z

            site_name(moltype,site_index) = atomic_name
         end do 

      end do 

      nmoltypes = moltype

!     read in model parameters
      call read_multipoles
      call read_model

!     translate the reference molecules to their centers of mass
      call translate_to_center_of_mass()

      call read_input_file(ios)

      end subroutine read_molecular_geometry
      

      subroutine read_input_file(ios)
      use common_data
      implicit none
      integer ios

!     reads INPUT_FILE

!     read in the number of molecules and the cell-vectors from INPUT_FILE
      call read_input_header(ios)

      if(ios.ne.0) return
!     read in the structure from INPUT_FILE
      if(read_mode.eq.'RIGID') then 
         call read_rigid_input
      else if(read_mode.eq.'XYZ') then 
         call read_xyz_input
      else if(read_mode.eq.'XYZ_MULTIPOLE') then 
         call read_xyz_multipoles_input
      endif

      end subroutine read_input_file
      

      subroutine read_input_header(ios)
      use common_data
      use parse_text
      implicit none
      character(len = 2048) :: cvec_buffer,buffer
      character(len = 2048), dimension(64) :: textlist      
      integer k,ios,nitems

      read(10,*,iostat = ios) nmol
      if(ios.ne.0) return
      read(10,73) cvec_buffer
      
 73   format(a256)
      call get_cvec(cvec_buffer)

!     read in the type of input file, either XYZ or rigid

      read(10,73) buffer
      buffer = adjustl(buffer)
      call split_text_to_list(buffer,textlist,nitems)
      if(upcase(textlist(1)).eq.'RIGID') then
         read_mode = 'RIGID'
         if(upcase(textlist(2)).eq.'XYZ') then
            rigid_input_mode = 'XYZ'
         else if(upcase(textlist(2)).eq.'FRAC') then
            rigid_input_mode = 'FRAC'
         else
            write(*,*) 'ERROR: INVALID RIGID INPUT MODE'
            stop
         endif
      else if(upcase(textlist(1)).eq.'XYZ') then
         read_mode = 'XYZ'
         if(upcase(textlist(2)).eq.'MULTIPOLE') then
            read_mode = 'XYZ_MULTIPOLE'
         endif
      else
         write(*,*) 'ERROR: INVALID READ MODE'
         stop
      endif


      end subroutine read_input_header

      subroutine get_cvec(buffer)
      use nr_mod
      use common_data
      use parse_text
      implicit none
      character(len = 2048) :: buffer
      real(8) :: amag,bmag,cmag,alpha,beta,gamma
      character(len = 2048), dimension(64) :: textlist
      character(len = 8) :: type
      logical okflag
      integer nitems

      periodic = .false.
      call split_text_to_list(buffer,textlist,nitems)
!     for blank line, which is an acceptable way to specify no pbc
      if(nitems.ne.1) then 
         read(textlist(1),*) type
         type = upcase(type)
         select case(type)
      case('XYZ') 
         periodic = .true.
         read(textlist(2),*) cvec(1,1)
         read(textlist(3),*) cvec(1,2)
         read(textlist(4),*) cvec(2,2)
         read(textlist(5),*) cvec(1,3)
         read(textlist(6),*) cvec(2,3)
         read(textlist(7),*) cvec(3,3)
      case('ANGLE')
         periodic = .true.
         read(textlist(2),*) amag            
         read(textlist(3),*) bmag
         read(textlist(4),*) cmag
         read(textlist(5),*) alpha
         read(textlist(6),*) beta
         read(textlist(7),*) gamma
         call get_cartesian_cvec(amag,bmag,cmag,alpha,beta,gamma,cvec)
      case('NOPBC')
         periodic = .false.
      end select
      endif
      
      if(periodic)  call mat3inverse(cvec,cveci,okflag)
      
      end subroutine get_cvec


      subroutine read_model
      use common_data
      use nr_mod
      use parse_text
      implicit none
      character(len=2048) :: buffer
      character(len = 16) :: molname
      character(len=2048), dimension(64) :: textlist
      integer j,m,nitem,nitems,nterms_default,ios,nterms,moltype,iatom
      integer natoms_in_mol,ns
      integer, dimension(16) :: power_list,power_list_default
      real(8), dimension(100,100) :: spline_data
      real(8), dimension(16) :: coeff_list
      integer n,nb,i,imax
      real(8) :: r0,u0,dudr0,dudrmax,uu,dudr
      logical:: error,match,has_mass,has_index
      real(8) :: r,rmin,rdamp,dmass
      real(8) :: rspline_max,rspline_min
      real(8), dimension(256) :: xspline,yspline,spline_coeff
      integer index,i_index,j_index,i_index_list,j_index_list,iparam,k
      integer nmass_atoms_in_mol
      character(len = 10) :: bondname
      character(len = 3) :: atname
      character(len = 32) :: token

      open(73,file = 'upair.dat')
      open(74,file = 'uspline.dat')
      
      nterms_default = 0 
      model_index_max = 0 

!     defaults
      ns = 100
      call set_nspline(ns)

      do while(.true.)
         read(30,'(A)',iostat = ios) buffer
         buffer = adjustl(buffer)
         call split_text_to_list(buffer,textlist,nitems)

         if(upcase(textlist(1)).eq.'MOLECULE') then 
            read(textlist(2),*) molname
            molname = upcase(molname)

            moltype = 0 
            do k = 1,nmoltypes
               if(molname.eq.molecule_name_list(k)) moltype = k
            end do 
            if(moltype.eq.0) then 
              write(*,*) 'ERROR IN MODEL FILE: &
     & MOLECULE NAME ',molname,' NOT RECOGNISED'
              stop               
            endif

            nmass_atoms_in_mol = 0 
            natoms_in_mol = mol_natoms(moltype)

            do iatom = 1,natoms_in_mol

               read(30, '(A)', iostat=ios) buffer
               if(ios.ne.0) then 
                  write(*,*) 'ERROR READING IN MODEL FILE'
                  stop
               endif

               call split_text_to_list(buffer,textlist,nitems)

               has_mass = .false.
               has_index = .false.

               nitem = 0 
               read(textlist(1),*) atname
               nitem = nitem + 1
               do while(.true.) 
                  nitem = nitem + 1
                  if(nitem.gt.nitems) exit

                  read(textlist(nitem),*) token
                  token = upcase(token)

                  select case(token)
                    case('TYPE')
                       nitem = nitem + 1
                       read(textlist(nitem),*) index
                       has_index = .true.
                    case('MASS')
                       nitem = nitem + 1
                       read(textlist(nitem),*) dmass
                       has_mass = .true.
                    case default
                      write(*,*) 'ERROR IN MODEL FILE: &
     &UNIDENTIFIED TOKEN ', token
                      stop
                  end select
               end do 

               if(.not.has_mass) then 
                  write(*,*) 'ERROR IN MODEL FILE: &
     &MASS MUST BE SPECIFIED FOR ALL ATOMS'
                  stop
               else if(.not.has_index) then 
                  write(*,*) 'ERROR IN MODEL FILE: &
     &INDEX MUST BE SPECIFIED FOR ALL ATOMS'
                  stop
               endif

               pair_index(iatom,moltype) = index
               mass(iatom,moltype) = dmass
               if(atname.ne.'X') nmass_atoms_in_mol = nmass_atoms_in_mol + 1
               if(index.gt.model_index_max) model_index_max = index
            end do 

             mol_nmassatoms(moltype) = nmass_atoms_in_mol
         else if(textlist(1).eq.'INTER') then
            j = 0
            do m = 2,nitems
               j = j + 1
               read(textlist(m),*) power_list_default(j)
            end do 
            nterms_default = j
            exit
         else if(textlist(1).eq.'DAMP') then
            read(textlist(2),*) rdamp
            call set_damp(rdamp)
            write(*,*)
            write(*,21) rdamp
 21         format('SETTING DAMP TO ',f14.6,' Angstroms')
         else if(textlist(1).eq.'NSPLINE') then 
            read(textlist(2),*) ns
            call set_nspline(ns)
            write(*,*) 
            write(*,22) nspline
 22         format('SETTING NSPLINE TO ',i6)
         else
            exit
            stop
         endif
      end do 

      nb = 0
      nmodel_param = 0
 
      do i_index = 1,model_index_max
         do j_index = i_index,model_index_max
            nb = nb + 1
            model_list(nb,1) = i_index
            model_list(nb,2) = j_index
!     defaults
            spline_data(1:4,nb) = 0.0d0
            model_nterms(nb) = nterms_default
            do n = 1,nterms_default
               model_power(n,nb) = power_list_default(n)
            end do 
            nmodel_param = nmodel_param + nterms_default
         end do 
      end do 

      nbondlist = nb

      nb = 0 
      iparam = 0 
      do while(.true.)
         read(30,'(A)',iostat = ios) buffer
         if(ios.ne.0) exit
         call split_text_to_list(buffer,textlist,nitems)
         if(textlist(1).eq.'PAIR') then
            read(textlist(2),*) i_index
            read(textlist(3),*) j_index
            nb = nb + 1
!     find out which bond this is
            match = .false.
            nb = 0 
            do m = 1,nbondlist
               i_index_list = model_list(m,1)
               j_index_list = model_list(m,2)
               if(i_index.eq.i_index_list.and.j_index.eq.j_index_list) match = .true.
               if(j_index.eq.i_index_list.and.i_index.eq.j_index_list) match = .true.               
               if(match) then
                  nb = m
                  exit
               endif
            end do 

            if(nb.ne.0) then 
               bond_name(nb) = bondname
            else
               write(*,*) 'ERROR: Bond ',i_index,j_index,&
     &' in MODEL file / INTER section is out of range'
               stop
            endif

         else if(textlist(1).eq.'POWER') then 
            j = 0
            do m = 2,nitems
               j = j + 1
               read(textlist(m),*) power_list(j)
               model_power(j,nb) = power_list(j)
            end do 
            nterms = j
            model_nterms(nb) = nterms
            nmodel_param = nmodel_param + nterms - nterms_default
         else if(textlist(1).eq.'COEFF') then 
            j = 0
            do m = 2,nitems
               j = j + 1
               read(textlist(m),*) coeff_list(j)
               model_coeff(j,nb) = coeff_list(j) / pefac
            end do 

         else if(upcase(textlist(1)).eq.'SPLINE') then 

            nitem = 1
            do while(.true.)
               nitem = nitem + 1
               if(nitem.gt.nitems) exit
               
               read(textlist(nitem),*) token
               token = upcase(token)

               select case(token)
                 case('RMIN')
                    nitem = nitem + 1
                    read(textlist(nitem),*) rspline_min
                 case('RMAX')
                    nitem = nitem + 1
                    read(textlist(nitem),*) rspline_max
                 case('R0')
                    nitem = nitem + 1
                    read(textlist(nitem),*) r0
                 case('U0')
                    nitem = nitem + 1
                    read(textlist(nitem),*) u0
                 case('DUDR0')
                    nitem = nitem + 1
                    read(textlist(nitem),*) dudr0
                 case default
                   write(*,*) 'ERROR IN MODEL FILE: &
     &UNIDENTIFIED TOKEN ', token
                   stop
               end select
            end do 

            spline_data(1,nb) = rspline_min
            spline_data(2,nb) = rspline_max
            spline_data(3,nb) = r0
            spline_data(4,nb) = u0 / pefac
            spline_data(5,nb) = dudr0 / pefac
         else if(textlist(1).eq.'   ') then
         else
            write(*,*) 'ERROR: DID NOT UNDERSTAND CARD ',trim(textlist(1)),&
     &' IN MODEL_FILE'
            stop
         endif
      end do 

!     assign a unique index to each coefficient
      iparam = 0 
      do nb = 1,nbondlist
         do j = 1,model_nterms(nb)
            iparam = iparam + 1
            model_param(j,nb) = iparam
         end do 
      end do 

!     CALCULATE SPLINE TERMS

      do nb = 1,nbondlist

         rspline_min = spline_data(1,nb) 
         rspline_max = spline_data(2,nb)
         rsplinemax(nb) = rspline_max
         r0 = spline_data(3,nb) 
         u0 = spline_data(4,nb) 
         dudr0 = spline_data(5,nb) 
         
         if(rspline_max.eq.0) cycle
         
         write(73,*) '#BOND',nb
         write(74,*) '#BOND',nb
         write(73,*)
         write(74,*)

         n = 1
         xspline(1) = r0
         yspline(1) = u0

         do i = 1,nspline - 1
            n = n + 1
            r = rspline_min + (rspline_max - rspline_min) * dble(i) / dble(nspline - 1)
            uu = 0.0d0 
            do j = 1,nterms
               uu = uu + model_coeff(j,nb) / r**power_list(j)
            end do 
            write(73,*) r,uu*pefac
            xspline(n)  = r
            yspline(n) = uu
         end do 
         dudrmax = 0.0d0 
         do j = 1,nterms
            dudrmax = dudrmax - (model_power(j,nb) * (model_coeff(j,nb)) / rspline_max**power_list(j))/rspline_max
         end do 
         call spline(xspline,yspline,nspline,dudr0,dudrmax,spline_coeff)         
         do n = 1,nspline
            bond_splinecoeff(n,nb) = spline_coeff(n)
            bond_xspline(n,nb) = xspline(n)
            bond_yspline(n,nb) = yspline(n)
         end do 
         
!     print out spline fits
         imax = 256
         rmin = 0.1d0
         do i = 1,imax
            r = rmin + (rspline_max - rmin) * dble(i) / dble(imax)
            call splint(xspline,yspline,spline_coeff,nspline,r,uu,dudr)      
            write(74,*) r,uu*pefac
         end do 
         call flush(73)
         call flush(74)
      end do 

      end subroutine read_model

      subroutine read_multipoles
      use common_data
      use parse_text
      implicit none
      integer imtype,imoltype,jmoltype,natoms_in_mol,iatom,index
      integer k,moltype,ios,conformer,nmult
      real(8) :: epsilon,sigma,dmass
      real(8) :: param_val,power
      integer, dimension(16) :: power_list,power_list_default
      real(8) :: rdamp
      real(8) :: q20,q21c,q21s,q22c,q22s
      real(8) :: q30,q31c,q31s,q32c,q32s,q33c,q33s
      real(8) :: trace
      character(len = 3) atomic_name
      character(len = 16) molecule_name,imol_name,jmol_name,aaa
      character(len=2048) :: buffer
      character(len=2048), dimension(64) :: textlist
      integer:: natoms_in_imol,natoms_in_jmol,jatom,imol,jmol,nterms,nterms_default
      integer j,i,n,imax,nbond,nb,m,mm,iparam,nitem,nitems
      integer i1,i2
      integer ivec,jvec,kvec
      real(8), dimension(16) :: multipole_vec
      real(8) :: qxx,qxy,qxz,qyy,qyz,qzz
      real(8) :: qxxx,qxxy,qxyy,qyyy,qxxz,qxyz,qyyz,qxzz,qyzz,qzzz      
      character(len = 32) :: token

      ios = 0 

!     read in molecular geometry file

      do while(.true.)
         read(29,'(A)',iostat = ios) buffer
         if(ios.ne.0) exit

         buffer = adjustl(buffer)
         call split_text_to_list(buffer,textlist,nitems)

!     defaults
         conformer = 1
         nitem = 0 
         do while(.true.) 
            nitem = nitem + 1
            if(nitem.gt.nitems) exit
            read(textlist(nitem),*,iostat = ios) token
            if(ios.ne.0) exit
            token = upcase(token)
            select case(token)
            case('MOLECULE')
               nitem = nitem + 1
               read(textlist(nitem),*) molecule_name
               molecule_name = upcase(molecule_name)
            case('CONFORMER')
               nitem = nitem + 1
               read(textlist(nitem),*) conformer
            case default
               write(*,*) 'ERROR IN MULTIPOLE FILE: &
     &UNIDENTIFIED TOKEN ', token
               stop
            end select 
         end do 

         if(ios.ne.0) exit
         moltype = 0 
         do k = 1,nmoltypes
            if(molecule_name.eq.molecule_name_list(k)) moltype = k
         end do 

         if(moltype.eq.0) then 
            write(*,*) 'ERROR IN MOLECULAR GEOMETRY FILE: &
     & MOLECULE NAME ',molecule_name,' NOT RECOGNISED'
            stop
         endif
         if(conformer.gt.mol_nconformers(moltype)) then 
            write(*,*) 'ERROR IN MULTIPOLE FILE: CONFORMER NUMBER ',conformer, &
     & ' > MAX NUMBER OF ',mol_nconformers(moltype),&
     & ' CONFORMERS FOR MOLECULE ',molecule_name
            stop
         endif

         natoms_in_mol = mol_natoms(moltype)

         do iatom = 1,natoms_in_mol
            read(29,'(A)',iostat = ios) buffer
            call split_text_to_list(buffer,textlist,nitems)
            
            read(textlist(1),*) atomic_name
            nmult = 0 
            do k = 2,nitems
               nmult = nmult + 1
               read(textlist(k),*) multipole_vec(nmult)
            end do 

            call assign_multipoles(nmult,iatom,conformer,moltype,multipole_vec)

         end do 
      end do 

      end subroutine read_multipoles

      subroutine read_rigid_input
      use common_data
      use nr_mod
      use parse_text
      implicit none
      character(len = 3) atomic_name
      character(len = 32) :: token
      integer natoms_in_mol
      real(8) :: x,y,z,phi,theta,psi
      real(8) :: amag,bmag,cmag,alpha,beta,gamma
      real(8) :: g1,g2,g3
      real(8), dimension(3) :: vec1,vec2
      character(len = 16) molecule_name
      character(len=2048) :: buffer
      character(len=2048), dimension(64) :: textlist
      integer ios,imol,k,moltype,nnn,imode,nitem,nitems,conformer
      logical okflag

      nmoloftype = 0  
      natoms = 0
      nmassatoms = 0 
      do imol = 1,nmol
         
         read(10,'(A)',iostat = ios) buffer
         call split_text_to_list(buffer,textlist,nitems)

!     default
         conformer = 1

         nitem = 0 
         do while(.true.)
            nitem = nitem + 1
            if(nitem.gt.nitems) exit

            read(textlist(nitem),*,iostat = ios) token
            if(ios.lt.0) then 
               write(*,*) 'ERROR: REACHED END OF INPUT FILE. WRONG NUMBER OF MOLECULES?'
               stop
            endif
            token = upcase(token)

            select case(token)
            case('MOLECULE')
               nitem = nitem + 1
               read(textlist(nitem),*) molecule_name
               molecule_name = upcase(molecule_name)
            case('CONFORMER')
               nitem = nitem + 1
               read(textlist(nitem),*) conformer
            case('COORD')
               nitem = nitem + 1
               read(textlist(nitem),*) g1
               nitem = nitem + 1
               read(textlist(nitem),*) g2
               nitem = nitem + 1
               read(textlist(nitem),*) g3
               nitem = nitem + 1
               read(textlist(nitem),*) phi
               nitem = nitem + 1
               read(textlist(nitem),*) theta
               nitem = nitem + 1
               read(textlist(nitem),*) psi
            case default
               write(*,*) 'ERROR IN INPUT FILE: &
     &UNIDENTIFIED TOKEN ', token
               stop
            end select
         end do 


!     determine the molecule's moltype index by comparing the name to the 
!     list of names

         moltype = 0
         do k = 1,nmoltypes
            if(molecule_name.eq.molecule_name_list(k)) moltype = k
         end do 
         if(moltype.eq.0) then 
            write(*,*) 'ERROR IN INPUT FILE: MOLECULE NAME ',molecule_name,' NOT RECOGNISED'
            stop
         endif
         if(conformer.gt.mol_nconformers(moltype)) then 
            write(*,*) 'ERROR IN INPUT FILE: CONFORMER NUMBER ',conformer, &
     & ' > MAX NUMBER OF ',mol_nconformers(moltype),&
     & ' CONFORMERS FOR MOLECULE ',molecule_name
            stop
         endif

         mol_conformer(imol) = conformer

         nmoloftype(moltype) = nmoloftype(moltype) + 1

         natoms = natoms + mol_natoms(moltype)
         nmassatoms = nmassatoms + mol_nmassatoms(moltype)

         if(periodic) then
            if(rigid_input_mode.eq.'FRAC') then
               rfrac(1,imol) = g1
               rfrac(2,imol) = g2
               rfrac(3,imol) = g3
            else if(rigid_input_mode.eq.'XYZ') then
               vec1(1) = g1
               vec1(2) = g2
               vec1(3) = g3
               call mat3multvec(cveci,vec1,vec2)
               rfrac(1,imol) = vec2(1)
               rfrac(2,imol) = vec2(2)
               rfrac(3,imol) = vec2(3)
            else
               write(*,*) 'ERROR: INVALID RIGID INPUT MODE'
               stop
            endif
         else
            mol_com(1,imol) = g1
            mol_com(2,imol) = g2
            mol_com(3,imol) = g3
         endif
         mol_euler(1,imol) = phi
         mol_euler(2,imol) = theta
         mol_euler(3,imol) = psi
         mol_type(imol) = moltype
      end do 
      end subroutine read_rigid_input

      subroutine read_xyz_input
      use common_data
      use nr_mod
      use parse_text
      implicit none
      integer ios
      real(8), dimension(3,500) :: rr
      real(8) :: comx,comy,comz,totmass
      real(8), dimension(3) :: vec1,vec2
      real(8) :: x1,y1,z1,x2,y2,z2,x3,y3,z3
      real(8) :: ax,ay,az
      real(8) :: x21,y21,z21
      real(8), dimension(3,3) :: axes0,axes,axes0i,rotmat
      real(8) :: phi,theta,psi,stheta
      real(8) :: xdif,ydif,zdif
      real(8) :: rdis2,error2,error,error_min
      real(8), allocatable, dimension(:) :: error_conf
      real(8), allocatable, dimension(:,:) :: euler_conf
      real(8), dimension(3) :: r0,r1
      character(len = 16) molecule_name
      character(len = 3) atom_name
      character(len=2048) :: buffer
      character(len = 32) :: token
      character(len = 2048), dimension(64) :: textlist
      integer imol,moltype,natoms_in_mol,iatom
      integer k,anchor1,anchor2,anchor3,conformer,conf_best,nitems,nitem,ivec
      integer iconf_min,iconf_max
      logical okflag,find_conformer

      allocate(euler_conf(3,nconformer_max))
      allocate(error_conf(nconformer_max))

      write(*,*) 
      write(*,*) 'READING XYZ INPUT COORDINATES'

!     anchor1,anchor2,anchor3 are the indices of the particles which are used to construct the axes.

!     default values
      anchor1 = 1
      anchor2 = 2
      anchor3 = 3

      nmoloftype = 0  
      nmassatoms = 0 
      do imol = 1,nmol

         find_conformer = .true.

         read(10, '(A)', iostat=ios) buffer
         call split_text_to_list(buffer,textlist,nitems)

!     default
         conformer = 1
         nitem = 0 
         do while(.true.)
            nitem = nitem + 1
            if(nitem.gt.nitems) exit

            read(textlist(nitem),*) token
            token = upcase(token)

            select case(token)
            case('MOLECULE')
               nitem = nitem + 1
               read(textlist(nitem),*) molecule_name
               molecule_name = upcase(molecule_name)
            case('CONFORMER')
               nitem = nitem + 1
               read(textlist(nitem),*) conformer
               find_conformer = .false.
            case('ANCHORS')
               nitem = nitem + 1
               read(textlist(nitem),*) anchor1
               nitem = nitem + 1
               read(textlist(nitem),*) anchor2
               nitem = nitem + 1
               read(textlist(nitem),*) anchor3
            case default 
               write(*,*) 'ERROR IN INPUT FILE: &
     &UNIDENTIFIED TOKEN ', token
               stop
            end select
         end do 


!     determine the molecule's moltype index by comparing the name to the 
!     list of names

         moltype = 0
         do k = 1,nmoltypes
            if(molecule_name.eq.molecule_name_list(k)) moltype = k
         end do 

         if(moltype.eq.0) then 
            write(*,*) 'ERROR: molecule name ',molecule_name,&
    & 'in input file does not correspond to molecule in list'
            stop
         endif
         if(conformer.gt.mol_nconformers(moltype)) then 
            write(*,*) 'ERROR IN INPUT FILE: CONFORMER NUMBER ',conformer, &
     & ' > MAX NUMBER OF ',mol_nconformers(moltype),&
     & ' CONFORMERS FOR MOLECULE ',molecule_name
            stop
         endif
         if(anchor1.eq.anchor2.or.anchor1.eq.anchor3.or.anchor2.eq.anchor3) then 
            write(*,*) 'ERROR IN INPUT FILE: ANCHORS ',anchor1,anchor2,anchor3, &
     & ' CONTAIN DUPLICATES'
            stop
         endif

         mol_type(imol) = moltype
         nmoloftype(moltype) = nmoloftype(moltype) + 1
         natoms_in_mol = mol_nmassatoms(moltype)
         nmassatoms = nmassatoms + natoms_in_mol
         
         if(natoms_in_mol.lt.max(anchor1,anchor2,anchor3)) then 
            write(*,*) 'ERROR IN INPUT FILE: AT LEAST ONE&
     & ANCHOR IS LARGER THAN, ',natoms_in_mol,', THE NUMBER OF&
     & ATOMS IN THIS MOLECULE'
            stop
         endif

         if(conformer.gt.mol_nconformers(moltype)) then 
            write(*,*) 'ERROR: CONFORMER NUMBER',conformer,&
     &' IS LARGER THAN NUMBER OF CONFORMERS FOR THIS MOLECULE'
            stop
         endif

         comx = 0.0d0 
         comy = 0.0d0 
         comz = 0.0d0 
         totmass = 0.0d0 
         do iatom = 1,natoms_in_mol
            read(10,*,iostat = ios) atom_name,rr(1,iatom),rr(2,iatom),rr(3,iatom)
            if(ios.ne.0) then 
               write(*,*) 'ERROR READING INPUT FILE. &
     & PROBLEM WITH NUMBER OF MOLECULES?'
               stop
            endif
            comx = comx + rr(1,iatom) * mass(iatom,moltype)
            comy = comy + rr(2,iatom) * mass(iatom,moltype)
            comz = comz + rr(3,iatom) * mass(iatom,moltype)
            totmass = totmass + mass(iatom,moltype)
         end do 
         comx = comx / totmass
         comy = comy / totmass
         comz = comz / totmass

         if(periodic) then
!     find fractional coordinates of the COM
            mol_com(1,imol) = comx
            mol_com(2,imol) = comy
            mol_com(3,imol) = comz
            
            vec1(1) = comx
            vec1(2) = comy
            vec1(3) = comz

            call mat3multvec(cveci,vec1,vec2)

            rfrac(1,imol) = vec2(1)
            rfrac(2,imol) = vec2(2)
            rfrac(3,imol) = vec2(3)
         else
            mol_com(1,imol) = comx
            mol_com(2,imol) = comy
            mol_com(3,imol) = comz
         endif
         
!     reference vector used for diatomics to produce a valid third vector

         ax = 0.0d0
         ay = 0.0d0
         az = 1.0d0

!     find axes of molecule from first three atoms

         x1 = rr(1,anchor1)
         y1 = rr(2,anchor1)
         z1 = rr(3,anchor1)
         
         x2 = rr(1,anchor2)
         y2 = rr(2,anchor2)
         z2 = rr(3,anchor2)

         if(mol_nmassatoms(moltype).eq.2) then

            x21 = x2 - x1
            y21 = y2 - y1
            z21 = z2 - z1

            x3 = ay * z21 - az * y21
            y3 = az * x21 - ax * z21
            z3 = ax * y21 - ay * x21

         else
            x3 = rr(1,anchor3)
            y3 = rr(2,anchor3)
            z3 = rr(3,anchor3)
         endif

         call get_axes(x1,y1,z1,x2,y2,z2,x3,y3,z3,axes)

         error_min = 1.d+23

!     decide whether to search for best conformer, or use the one in the input file

         if(find_conformer) then 
            iconf_min = 1
            iconf_max = mol_nconformers(moltype)
         else
            iconf_min = conformer
            iconf_max = conformer
         endif

         
         conformer_loop: do conformer = iconf_min,iconf_max
         
!     find the axes for the reference molecule
         
         x1 = site_coord0(1,anchor1,conformer,moltype)
         y1 = site_coord0(2,anchor1,conformer,moltype)
         z1 = site_coord0(3,anchor1,conformer,moltype)
         
         x2 = site_coord0(1,anchor2,conformer,moltype)
         y2 = site_coord0(2,anchor2,conformer,moltype)
         z2 = site_coord0(3,anchor2,conformer,moltype)
         
         if(mol_nmassatoms(moltype).eq.2) then
            
            x21 = x2 - x1
            y21 = y2 - y1
            z21 = z2 - z1
            
            x3 = ay * z21 - az * y21
            y3 = az * x21 - ax * z21
            z3 = ax * y21 - ay * x21
            
         else
            x3 = site_coord0(1,anchor3,conformer,moltype)
            y3 = site_coord0(2,anchor3,conformer,moltype)
            z3 = site_coord0(3,anchor3,conformer,moltype)
         endif
         
         call get_axes(x1,y1,z1,x2,y2,z2,x3,y3,z3,axes0)
         
         call mat3inverse(axes0,axes0i,okflag)
         rotmat = matmul(axes,axes0i)
         
         call get_euler_from_rotmat(rotmat,phi,theta,psi)
         
         euler_conf(1,conformer) = phi
         euler_conf(2,conformer) = theta
         euler_conf(3,conformer) = psi
         
!     calculate the error with this conformer 
         
         call get_rotmat(phi,theta,psi,rotmat)         
         
         error2 = 0.0d0 
         do iatom = 1,natoms_in_mol
            do ivec = 1,3
               r0(ivec) = site_coord0(ivec,iatom,conformer,moltype)
            end do 
            r1 = matmul(rotmat,r0)
            xdif = rr(1,iatom) - comx
            ydif = rr(2,iatom) - comy
            zdif = rr(3,iatom) - comz
            error2 = error2 + (xdif - r1(1))**2 + (ydif - r1(2))**2 + (zdif - r1(3))**2
         end do 
         error = dsqrt(error2/dble(natoms_in_mol))
         error_conf(conformer) = error
         if(error.lt.error_min) then 
            error_min = error
            conf_best = conformer
         endif
         
      end do conformer_loop
      
      if(find_conformer) write(*,*) 'MOLECULE',imol,' BEST CONFORMER = '&
     &     ,conf_best,' RMS ERROR = ',error_conf(conf_best)
      
      mol_euler(1,imol) = euler_conf(1,conf_best)
      mol_euler(2,imol) = euler_conf(2,conf_best)
      mol_euler(3,imol) = euler_conf(3,conf_best)
      
      mol_conformer(imol) = conf_best
      
      end do 

      end subroutine read_xyz_input

      subroutine read_xyz_multipoles_input
      use common_data
      use parse_text
      implicit none
      integer natoms_in_mol,nitems,moltype,k,ios,imol,conformer,iatom
      character(len=2048) :: buffer
      character(len = 16) molecule_name
      character(len = 2048), dimension(64) :: textlist
      real(8), dimension(:,:), allocatable :: rr
      real(8), dimension(3) :: vec1,vec2
      real(8) :: totmass
      real(8) :: comx,comy,comz
      real(8), dimension(16) :: multi_vec
      real(8) :: qxx,qxy,qxz,qyy,qyz,qzz
      real(8) :: qxxx,qxxy,qxyy,qyyy,qxxz,qxyz,qyyz,qxzz,qyzz,qzzz
      character(len = 3) atom_name
      integer ivec,jvec,kvec,nmult

      allocate(rr(3,natom_max))

      max_rank = -1
      site_rank = -1 

      nmoloftype = 0  
      nmassatoms = 0 
      do imol = 1,nmol
         read(10, '(A)', iostat=ios) buffer
         call split_text_to_list(buffer,textlist,nitems)
         read(textlist(1),*) molecule_name
         molecule_name = upcase(molecule_name)
         conformer = 1
         if(nitems.gt.1) read(textlist(2),*) conformer

!     determine the molecule's moltype index by comparing the name to the 
!     list of names

         moltype = 0 
         do k = 1,nmoltypes
            if(molecule_name.eq.molecule_name_list(k)) moltype = k
         end do 
         if(moltype.eq.0) then 
            write(*,*) 'ERROR: molecule name ',molecule_name,&
    & 'in input file does not correspond to molecule in list'
            stop
         endif
         mol_type(imol) = moltype
         mol_conformer(imol) = conformer
         nmoloftype(moltype) = nmoloftype(moltype) + 1
         natoms_in_mol = mol_nmassatoms(moltype)
         nmassatoms = nmassatoms + natoms_in_mol

         comx = 0.0d0 
         comy = 0.0d0 
         comz = 0.0d0 
         totmass = 0.0d0 
         do iatom = 1,natoms_in_mol

            read(10,'(A)',iostat = ios) buffer
            call split_text_to_list(buffer,textlist,nitems)

            read(textlist(1),*,iostat = ios) atom_name
            read(textlist(2),*,iostat = ios) rr(1,iatom)
            read(textlist(3),*,iostat = ios) rr(2,iatom)
            read(textlist(4),*,iostat = ios) rr(3,iatom)
            
            nmult = 0 
            do k = 5,nitems
               nmult = nmult + 1
               read(textlist(k),*) multi_vec(nmult)
            end do 

            lab_coord(1,iatom,imol) = rr(1,iatom)
            lab_coord(2,iatom,imol) = rr(2,iatom)
            lab_coord(3,iatom,imol) = rr(3,iatom)

            comx = comx + rr(1,iatom) * mass(iatom,moltype)
            comy = comy + rr(2,iatom) * mass(iatom,moltype)
            comz = comz + rr(3,iatom) * mass(iatom,moltype)

            totmass = totmass + mass(iatom,moltype)

!     here we're only interested in the lab multipoles. Subroutine call will assign 
!     multipoles to site arrays. So, follow up by assigning them to lab arrays.

            call assign_multipoles(nmult,iatom,conformer,moltype,multi_vec)

            lab_charge(iatom,imol) = site_charge0(iatom,conformer,moltype)

            if(max_rank.ge.1) then 
               do ivec = 1,3
                  lab_dipole(ivec,iatom,imol) = site_dipole0(ivec,iatom,conformer,moltype)
                  if(max_rank.ge.2) then 
                     do jvec = 1,3
                        lab_quad(ivec,jvec,iatom,imol) = site_quad0(ivec,jvec,iatom,conformer,moltype)
                        if(max_rank.ge.3) then 
                           do kvec = 1,3
                              lab_oct(ivec,jvec,kvec,iatom,imol) = site_oct0(ivec,jvec,kvec,iatom,conformer,moltype)
                           end do 
                        endif
                     end do 
                  endif
               end do 
            endif

            if(max_rank.ge.2) then 
               qxx = lab_quad(1,1,iatom,imol)
               qxy = lab_quad(1,2,iatom,imol)
               qxz = lab_quad(1,3,iatom,imol)
               qyy = lab_quad(2,2,iatom,imol)
               qyz = lab_quad(2,3,iatom,imol)
               qzz = lab_quad(3,3,iatom,imol)
               call get_lab_spherical2(iatom,imol,qxx,qxy,qxz,qyy,qyz,qzz)
            endif

            if(max_rank.ge.3) then 
               qxxx = lab_oct(1,1,1,iatom,imol)
               qxxy = lab_oct(1,1,2,iatom,imol)
               qxyy = lab_oct(1,2,2,iatom,imol)
               qyyy = lab_oct(2,2,2,iatom,imol)
               qxxz = lab_oct(1,1,3,iatom,imol)
               qxyz = lab_oct(1,2,3,iatom,imol)
               qyyz = lab_oct(2,2,3,iatom,imol)
               qxzz = lab_oct(1,3,3,iatom,imol)
               qyzz = lab_oct(2,3,3,iatom,imol)
               qzzz = lab_oct(3,3,3,iatom,imol)
               call get_lab_spherical3(iatom,imol,qxxx,qxxy,qxyy,qyyy,qxxz,qxyz,qyyz,qxzz,qyzz,qzzz)
            endif

         end do 

         comx = comx / totmass
         comy = comy / totmass
         comz = comz / totmass

         mol_com(1,imol) = comx
         mol_com(2,imol) = comy
         mol_com(3,imol) = comz

!     calculate the relative coordinates

         do iatom = 1,natoms_in_mol
            lab_relative_coord(1,iatom,imol) = lab_coord(1,iatom,imol) - comx
            lab_relative_coord(2,iatom,imol) = lab_coord(2,iatom,imol) - comy
            lab_relative_coord(3,iatom,imol) = lab_coord(3,iatom,imol) - comz
         end do 

      end do 

      end subroutine read_xyz_multipoles_input


      subroutine get_axes(x1,y1,z1,x2,y2,z2,x3,y3,z3,axes)
      implicit none
      real(8) :: x1,y1,z1,x2,y2,z2,x3,y3,z3
      real(8) :: xx,xy,xz,yx,yy,yz,zx,zy,zz
      real(8) :: x21,y21,z21,x31,y31,z31
      real(8) :: r2xr3x,r2xr3y,r2xr3z
      real(8) :: r
      real(8), dimension(3,3) :: axes

      x21 = x2 - x1
      y21 = y2 - y1
      z21 = z2 - z1

      x31 = x3 - x1
      y31 = y3 - y1
      z31 = z3 - z1


!     x axis

      xx = x21
      xy = y21
      xz = z21

      r = dsqrt(xx**2 + xy**2 + xz**2)

      xx = xx/r
      xy = xy/r
      xz = xz/r

!     y axis given by r21 X r31

      r2xr3x = (y21 * z31 - y31 * z21)
      r2xr3y = (z21 * x31 - z31 * x21)
      r2xr3z = (x21 * y31 - x31 * y21)

      yx = r2xr3x
      yy = r2xr3y
      yz = r2xr3z

      r = dsqrt(yx**2 + yy**2 + yz**2)
      yx = yx / r
      yy = yy / r
      yz = yz / r

!     z axis given by x X y

      zx = xy * yz - xz * yy
      zy = xz * yx - xx * yz
      zz = xx * yy - xy * yx

      axes(1,1) = xx
      axes(2,1) = xy
      axes(3,1) = xz

      axes(1,2) = yx
      axes(2,2) = yy
      axes(3,2) = yz

      axes(1,3) = zx
      axes(2,3) = zy
      axes(3,3) = zz

      end subroutine get_axes


      subroutine translate_to_center_of_mass()
!     translates the molecules to their centers of masses
      use common_data      
      implicit none
      integer imoltype,iatom,conformer
      real(8) xcom,ycom,zcom,totmass

      do imoltype = 1,nmoltypes
         do conformer = 1,mol_nconformers(imoltype)
            xcom = 0.0d0 
            ycom = 0.0d0 
            zcom = 0.0d0 
            totmass = 0.0d0 
            
            do iatom = 1,mol_natoms(imoltype)
               if(site_name(imoltype,iatom).ne.'X') then
                  xcom = xcom + site_coord0(1,iatom,conformer,imoltype) * mass(iatom,imoltype)
                  ycom = ycom + site_coord0(2,iatom,conformer,imoltype) * mass(iatom,imoltype)
                  zcom = zcom + site_coord0(3,iatom,conformer,imoltype) * mass(iatom,imoltype)
                  totmass = totmass + mass(iatom,imoltype)
               endif
            end do 
            
            xcom = xcom / totmass
            ycom = ycom / totmass
            zcom = zcom / totmass
            
            do iatom = 1,mol_natoms(imoltype)
               site_coord0(1,iatom,conformer,imoltype) = site_coord0(1,iatom,conformer,imoltype) - xcom
               site_coord0(2,iatom,conformer,imoltype) = site_coord0(2,iatom,conformer,imoltype) - ycom
               site_coord0(3,iatom,conformer,imoltype) = site_coord0(3,iatom,conformer,imoltype) - zcom
         end do 
      end do 
      end do 
      
      end subroutine translate_to_center_of_mass
      
      subroutine init_random_seed()
      integer :: i, n, clock
      integer, dimension(:), allocatable :: seed

      call random_seed(size = n)
      allocate(seed(n))

      call system_clock(count=clock)

      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      call random_seed(put = seed)

      deallocate(seed)
      end subroutine init_random_seed

      subroutine CALC_ELEC()
!----------------------------------------------
!     CALCULATES THE ELECTROSTATIC ENERGY
!----------------------------------------------
      use common_data
      implicit none
      integer :: imol,iatom,natoms_in_mol,moltype
      integer :: ivec,jvec,kvec
      real(8) :: dix,diy,diz,trace
      real(8) :: qxx,qxy,qxz,qyy,qyz,qzz
      real(8) :: qxxx,qxxy,qxyy,qyyy,qxxz,qxyz,qyyz,qxzz,qyzz,qzzz
      real(8) :: q20,q21c,q21s,q22c,q22s
      real(8) :: q30,q31c,q31s,q32c,q32s,q33c,q33s

      ucharge = 0.0d0
      udipole = 0.0d0 
      uquad = 0.0d0 
      uoct = 0.0d0 

     call symmetrise_fields()

      do imol = 1,nmol
         moltype = mol_type(imol)
         natoms_in_mol = mol_natoms(moltype)
         do iatom = 1,natoms_in_mol
            ucharge = ucharge + 0.5d0 * lab_charge(iatom,imol) * cfield(iatom,imol) * fac(2) 
            do ivec = 1,3
               udipole = udipole + 0.5d0 * lab_dipole(ivec,iatom,imol) * dfield(ivec,iatom,imol) * fac(2)
               do jvec = 1,3
                  uquad = uquad + 0.5d0 * lab_quad(ivec,jvec,iatom,imol) * qfield(ivec,jvec,iatom,imol) * fac(2)
                  do kvec = 1,3
                     uoct = uoct + 0.5d0 * lab_oct(ivec,jvec,kvec,iatom,imol) * ofield(ivec,jvec,kvec,iatom,imol) * fac(2)
                  end do 
               end do 
            end do 
         end do 
      end do 

      end subroutine CALC_ELEC
      

      subroutine PRESS()
      use common_data
      implicit none
      real(8) :: vol
      real(8) :: ax,ay,az,bx,by,bz,cx,cy,cz

      vol = cvec(1,1) * cvec(2,2) * cvec(3,3)
      upv = upv + pressure * vol

      dudcvec_pv(1,1) = cveci(1,1) * pressure * vol
      dudcvec_pv(2,2) = cveci(2,2) * pressure * vol
      dudcvec_pv(3,3) = cveci(3,3) * pressure * vol

      end subroutine press

      subroutine calc_press_and_dens
      use common_data
      implicit none
      real(8) :: vol
      real(8), dimension(3,3) :: dudcvec0

      dudcvec0 = dudcvec - dudcvec_pv

      vol = cvec(1,1)*cvec(2,2)*cvec(3,3)
      
      pressure_internal = -(dudcvec0(1,1) * cvec(1,1) &
     &                    + dudcvec0(1,2) * cvec(1,2) &
     &                    + dudcvec0(2,2) * cvec(2,2) &
     &                    + dudcvec0(1,3) * cvec(1,3) &
     &                    + dudcvec0(2,3) * cvec(2,3) &
     &                    + dudcvec0(3,3) * cvec(3,3)) / (3.0d0 * vol)
      
      density = cell_mass / vol

      end subroutine calc_press_and_dens

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


      function gauss(rand,sigma)
      implicit none
      logical :: accept
      real(8) :: gauss,sigma,x,y,rand,ee

      accept = .false.
      do while(.not.accept)
         call random_number(rand)
         x = (rand - 0.5d0) * 10.d0 * sigma
         call random_number(rand)
         y = rand
         ee = dexp(-x*x/(2.d0 * sigma * sigma))
         if(y < ee) accept = .true.

      end do 
      gauss = x

      end function gauss
      

      subroutine check_for_collisions(rmin,too_close,dismin,printrdf)
      use common_data
      use nr_mod
      implicit none
      integer imol,jmol,imoltype,jmoltype,iatom,jatom
      integer ibin,ibinmax
      integer natoms_in_imol,natoms_in_jmol
      real(8), dimension(0:2048) :: rdf
      real(8) :: rmin,dismin,dismin2,dr
      real(8) :: xi,yi,zi,xj,yj,zj
      real(8) :: rdis,rdis2,rmin2
      real(8) :: xdif,ydif,zdif
      real(8) :: dcellx,dcelly,dcellz
      real(8) :: dcyz,dcxz,dcxy
      real(8) :: ax,ay,az,bx,by,bz,cx,cy,cz
      real(8) :: vol
      integer :: jcell1max,jcell2max,jcell3max
      integer :: jcell1,jcell2,jcell3
      real(8), dimension(3,3) :: scvec,scveci
      real(8) :: astar_mag,bstar_mag,cstar_mag
      logical :: same_mol,too_close,okflag,same_cell,same_atom,loop,printrdf

      rmin2 = rmin*rmin
      dismin2 = 100.0d0 

      ibinmax = 1024
      rdf = 0.0d0 

!     finds the size of the supercell which contains the cut-off sphere

      if(periodic) then 

!     finds the size of the supercell which contains the cut-off sphere
         vol = cvec(1,1) * cvec(2,2) * cvec(3,3)         
         
!     magnitude of the reciprocal lattice vectors
         astar_mag = dsqrt(cveci(1,1)**2 + cveci(1,2)**2 + cveci(1,3)**2)
         bstar_mag = dsqrt(cveci(2,1)**2 + cveci(2,2)**2 + cveci(2,3)**2)
         cstar_mag = dsqrt(cveci(3,1)**2 + cveci(3,2)**2 + cveci(3,3)**2)

         jcell1max = int(2.0d0 * rcut * astar_mag) + 1
         jcell2max = int(2.0d0 * rcut * bstar_mag) + 1
         jcell3max = int(2.0d0 * rcut * cstar_mag) + 1

!     supercell cell-vectors

         scvec(:,1) = (jcell1max) * cvec(:,1)
         scvec(:,2) = (jcell2max) * cvec(:,2)
         scvec(:,3) = (jcell3max) * cvec(:,3)
         
         call mat3inverse(scvec,scveci,okflag)      
      else
         jcell1max = 1
         jcell2max = 1
         jcell3max = 1
      endif

      too_close = .false.
      do imol = 1,nmol
         imoltype = mol_type(imol)
         natoms_in_imol = mol_natoms(imoltype)
         do iatom = 1,natoms_in_imol

            xi = lab_coord(1,iatom,imol)
            yi = lab_coord(2,iatom,imol)
            zi = lab_coord(3,iatom,imol)

            do jmol = imol,nmol
               jmoltype = mol_type(jmol)
               natoms_in_jmol = mol_natoms(jmoltype)
               do jatom = 1,natoms_in_jmol

                  do jcell3 = 0,jcell3max-1
                     dcellz = jcell3 * cvec(3,3)
                     dcyz = jcell3 * cvec(2,3)
                     dcxz = jcell3 * cvec(1,3)
                     do jcell2 = 0,jcell2max-1
                        dcelly = jcell2 * cvec(2,2) + dcyz
                        dcxy = jcell2 * cvec(1,2)
                        do jcell1 = 0,jcell1max-1
                           dcellx = jcell1 * cvec(1,1) + dcxy + dcxz

                           xj = lab_coord(1,jatom,jmol)
                           yj = lab_coord(2,jatom,jmol)
                           zj = lab_coord(3,jatom,jmol)


                           same_atom = .false.
                           if(imol.eq.jmol.and.iatom.eq.jatom) same_atom = .true.

!     decide whether to loop
                           loop = .false.
                           same_mol = .false.
                  
                           if(periodic) then
                              same_cell = .false.
                              if(jcell1.eq.0.and.jcell2.eq.0.and.jcell3.eq.0) same_cell = .true.
                              if(same_cell.and.(imol.eq.jmol)) same_mol = .true.
                              
                              if((jmol.gt.imol).or.(imol.eq.jmol.and.jatom.gt.iatom)) loop = .true.
                              
                              if(same_atom) then
                                 loop = .false.
                                 if(.not. same_cell) loop = .true.
                              endif
                           else
                              if(jmol.gt.imol) loop = .true.
                           endif
                           if(loop) then
                              
                              xdif = xj - xi + dcellx
                              ydif = yj - yi + dcelly
                              zdif = zj - zi + dcellz
                              
                              if(periodic) then
                                 call nearest(xdif,ydif,zdif,scvec,scveci)
                              endif
                              
                              if(same_atom) then
                                 if((xdif.gt.0).or.(xdif.eq.0.and.ydif.gt.0)&
     &                                .or.(xdif.eq.0.and.ydif.eq.0.and.zdif.gt.0)) then 
                              else
                                 cycle
                              endif
                           endif
                           
                           
                           rdis2 = xdif**2 + ydif**2 + zdif**2
                           
                           if(.not.same_mol.and.rdis2.le.rcut*rcut) then 
                              rdis = dsqrt(rdis2)
                              ibin = anint((rdis / rcut) * ibinmax)
                              rdf(ibin) = rdf(ibin) + 1.0d0 
                           
                              if(rdis2.le.dismin2) then
                                 dismin2 = rdis2
                              endif
                              if(rdis2.lt.rmin2) then
                                 too_close = .true.
                              endif
                           endif

                           
                        endif
                        
                     end do
                  end do 
               end do 
               
               
            end do 
         end do 
      end do 
      end do 
      dismin = dsqrt(dismin2)

      if(printrdf) then 
         write(*,*) 'WRITING RDF'
         write(77,*) '#rdf',ibinmax-1,utotal * pefac,pressure_internal*fac(7),density*fac(8)
         dr = rcut / dble(ibinmax)
         do ibin = 1,ibinmax-1
            rdis = rcut * dble(ibin) / dble(ibinmax)
            write(77,*) rdis,rdf(ibin)/(2.0d0 * pi * rdis*rdis * dr * dble(natoms ** 2)/vol)
         end do 
         write(77,*)
         call flush(77)
      endif

      end subroutine check_for_collisions

      
      subroutine nearest(xdif,ydif,zdif,cmat,cmati)
      implicit none
      real(8),dimension(3,3) :: cmat,cmati
      real(8) :: xdif,ydif,zdif
      real(8),dimension(3) :: vec1,vec2

!     applies the minimum image convention for triclinic cells
!     note, assumes that both cmat and cmati are upper triangular

      vec1(1) = anint(cmati(1,1) * xdif + cmati(1,2) * ydif + cmati(1,3) * zdif)
      vec1(2) = anint(cmati(2,2) * ydif + cmati(2,3) * zdif)
      vec1(3) = anint(cmati(3,3) * zdif)

      xdif = xdif - (cmat(1,1) * vec1(1) + cmat(1,2) * vec1(2) + cmat(1,3) * vec1(3))
      ydif = ydif - (cmat(2,2) * vec1(2) + cmat(2,3) * vec1(3))
      zdif = zdif - (cmat(3,3) * vec1(3))

      end subroutine nearest

      subroutine pair_longrange
!     calculates the long range correction to the lennard jones interaction
      use common_data
      implicit none
      integer imoltype,jmoltype,iatom,jatom,natoms_in_imol,natoms_in_jmol
      real(8) :: vol
      real(8) :: prefac,mult
      integer :: nmoloftype_i,nmoloftype_j
      integer :: nterms,i_index_list,j_index_list,m,imoltype_list,jmoltype_list,i_index,j_index,j,m_index
      real(8) :: upair_lr
      real(8) :: rc
      real(8) :: coeff
      integer :: power,nb,mm
      logical :: match

!     because of the switch function, use 0.95 rcut instead of rcut
      rc = 0.95d0 * rcut

      vol = cvec(1,1) * cvec(2,2) * cvec(3,3)
      prefac = 4.0d0 * pi / vol

      upair_lr = 0.0d0 
      do imoltype = 1,nmoltypes
         natoms_in_imol = mol_natoms(imoltype)
         nmoloftype_i = nmoloftype(imoltype)
         do iatom = 1,natoms_in_imol
            
            i_index = pair_index(iatom,imoltype)
            
            do jmoltype = imoltype,nmoltypes
               natoms_in_jmol = mol_natoms(jmoltype)
               nmoloftype_j = nmoloftype(jmoltype)
               do jatom = 1,natoms_in_jmol
                  
                  if(imoltype.ne.jmoltype.or.jatom.ge.iatom) then 
                     mult = 1.0d0 
                     if(imoltype.eq.jmoltype.and.iatom.eq.jatom) mult = 0.50d0 
                     
                     j_index = pair_index(jatom,jmoltype)
                     match = .false.
                     mm = 0 
                     do m = 1,nbondlist

                        i_index_list = model_list(m,1)
                        j_index_list = model_list(m,2)
                           
                        if(i_index.eq.i_index_list.and.j_index.eq.j_index_list) match = .true.
                        if(j_index.eq.i_index_list.and.i_index.eq.j_index_list) match = .true.
                        if(match) then
                           nb = m
                           exit
                        endif
                     end do 
                     
                     if(nb.ne.0) then 
                        nterms = model_nterms(nb)
                     else
                        nterms = 0
                     endif

                     do j = 1,nterms
                        power = model_power(j,nb)
                        coeff = model_coeff(j,nb)
                        
                        upair_lr = upair_lr + mult * prefac * dble(nmoloftype_i * nmoloftype_j) &
     &                       * (coeff * rc**(3-power))/(power - 3.0d0)                            
                     end do 
                     
                  endif
               end do
            end do
         end do 
      end do 

      upair = upair + upair_lr

      dudcvec_pair(1,1) = dudcvec_pair(1,1) - cveci(1,1) * upair_lr 
      dudcvec_pair(2,2) = dudcvec_pair(2,2) - cveci(2,2) * upair_lr 
      dudcvec_pair(3,3) = dudcvec_pair(3,3) - cveci(3,3) * upair_lr 

      end subroutine pair_longrange

      subroutine set_swap_conformer_chance(schance)
      use common_data
      implicit none
      real(8) :: schance
      swap_conformer_chance = schance
      end subroutine set_swap_conformer_chance

      subroutine set_swap_chance(schance)
      use common_data
      implicit none
      real(8) :: schance
      swap_chance = schance
      end subroutine set_swap_chance

      subroutine set_damp(rdamp)
      use common_data
      implicit none
      real(8) :: rdamp
      damp_width = rdamp
!     default value, if not set by user
      if(.not.rdamp_userset) rdamp_cutoff = 4.0d0 * damp_width
      eps_damp = 1.0d0/(dsqrt(2.0d0) * damp_width)
      eps_damp_sqrtpii = 1.0d0/(eps_damp*dsqrt(pi))
      end subroutine set_damp

      subroutine set_nspline(ns)
      use common_data
      implicit none
      integer ns

      nspline = ns

      if(nspline.gt.nspline_max) then 
         write(*,*) 'ERROR: nspline ',nspline,&
     & ' > nspline_max in size.dat = ',nspline_max
         stop
      endif
      end subroutine set_nspline

      subroutine set_pressure(press)
      use common_data
      implicit none
      real(8) :: press

      pressure = press / fac(7)
      end subroutine set_pressure

      subroutine set_tolerance(tolerance)
      use common_data
      implicit none
      real(8) :: tolerance
      mintol = tolerance / pefac
      end subroutine set_tolerance

      subroutine get_cell_mass()
      use common_data
      implicit none
      integer imol,natoms_in_mol,iatom,imoltype

      cell_mass = 0.0d0 
      do imol = 1,nmol
         imoltype = mol_type(imol)
         natoms_in_mol = mol_natoms(imoltype)
         do iatom = 1,natoms_in_mol
            cell_mass = cell_mass + mass(iatom,imoltype)
         end do 
      end do 

      end subroutine get_cell_mass
      
      subroutine set_r_closest(r)
      use common_data
      implicit none
      real(8) :: r

      r_closest = r
      end subroutine set_r_closest

      subroutine get_cartesian_cvec(amag,bmag,cmag,alpha,beta,gamma,cvec0)
      use common_data
      implicit none
      real(8) :: alpha,beta,gamma,amag,bmag,cmag
      real(8) :: alpha_rad,beta_rad,gamma_rad
      real(8), dimension(6) :: param
      real(8) :: angfac
      real(8), dimension(3,3) :: cvec0

      angfac = 180.0d0 / pi

      alpha_rad = alpha / angfac
      beta_rad = beta / angfac
      gamma_rad = gamma / angfac

      cvec0 = 0.0d0

      cvec0(1,1) = amag
      cvec0(1,2) = bmag * dcos(gamma_rad)
      cvec0(2,2) = dsqrt(bmag**2 - cvec0(1,2)**2)
      cvec0(1,3) = cmag * dcos(beta_rad)
      cvec0(2,3) = (bmag * cmag * dcos(alpha_rad) - cvec0(1,2) * cvec0(1,3))/cvec0(2,2)
      cvec0(3,3) = dsqrt(cmag**2 - cvec0(1,3)**2 - cvec0(2,3)**2)

      end subroutine get_cartesian_cvec

      subroutine get_angle_cvec(cvec0,amag,bmag,cmag,alpha,beta,gamma)
      use common_data      
      implicit none
      real(8) :: angfac
      real(8) :: ax,ay,az,bx,by,bz,cx,cy,cz
      real(8) :: amag,bmag,cmag,alpha,beta,gamma
      real(8), dimension(6) :: param,param2
      real(8), dimension(3,3) :: cvec0

      angfac = 180.0d0 / pi

      ax = cvec0(1,1)
      ay = 0.0d0
      az = 0.0d0

      bx = cvec0(1,2)
      by = cvec0(2,2)
      bz = 0.0d0

      cx = cvec0(1,3)
      cy = cvec0(2,3)
      cz = cvec0(3,3)

      amag = dsqrt(ax**2 + ay**2 + az**2)
      bmag = dsqrt(bx**2 + by**2 + bz**2)
      cmag = dsqrt(cx**2 + cy**2 + cz**2)

      alpha = dacos((bx*cx + by*cy + bz*cz)/(bmag*cmag)) * angfac
      beta = dacos((cx*ax + cy*ay + cz*az)/(cmag*amag)) * angfac
      gamma = dacos((ax*bx + ay*by + az*bz)/(amag*bmag)) * angfac

      end subroutine get_angle_cvec

      subroutine print_forces_and_torques_and_derivs(ifile)
      use common_data
      implicit none
      integer j,imol,ivec,ifile

      do imol = 1,nmol
         write(ifile,*) imol,forcemol(1,imol)*pefac,forcemol(2,imol)*pefac,forcemol(3,imol)*pefac&
     &                  ,torquemol(1,imol)*pefac,torquemol(2,imol)*pefac,torquemol(3,imol)*pefac
      end do 

      end subroutine print_forces_and_torques_and_derivs


      subroutine calc_pair_from_polynomial(rdis,ri,dupair,rc0,rc1,switch,dswitchdr,nb,nterms)
      use common_data
      implicit none
      real(8) :: upair0,dupair0,dupair,rc0,rc1
      real(8) :: rdis,ri,r2i,rpi
      real(8) :: sigma6,sigma12,epsilon
      real(8) :: switch,dswitchdr
      integer :: j,nterms,power,nb,ibin,iparam
      real(8) :: coeff

      r2i = ri * ri
      dupair = 0.0d0 

      do j = 1,nterms
         power = model_power(j,nb)
         coeff = model_coeff(j,nb)
         iparam = model_param(j,nb)

         rpi = ri**power

         upair0 = coeff * rpi
!     dupair = -1/r dupair/dr
         dupair0 = power * coeff * rpi * r2i

         if(periodic) then
            if(rdis > rc0 .and. rdis < rc1) then 
               upair = upair + upair0 * switch
               dupair = dupair - upair0 * dswitchdr * ri + dupair0 * switch 
            else if(rdis < rc0) then
               upair = upair + upair0
               dupair = dupair + dupair0
            endif
         else
            upair = upair + upair0
            dupair = dupair + dupair0
         endif
      end do 

      end subroutine calc_pair_from_polynomial


      subroutine calc_pair_from_spline(rdis,ri,dupair,rc0,rc1,switch,dswitchdr,nb)
      use common_data
      use nr_mod

      implicit none
      real(8) :: rdis,ri,dupair,rc0,rc1,switch,dswitchdr
      real(8), dimension(256) :: spline_coeff,xspline,yspline
      real(8) :: u,dudr,upair0,dupair0
      integer n,nb

      do n = 1,nspline
         spline_coeff(n) = bond_splinecoeff(n,nb)
         xspline(n) = bond_xspline(n,nb)
         yspline(n) = bond_yspline(n,nb)
      end do 
      call splint(xspline,yspline,spline_coeff,nspline,rdis,u,dudr)      

      upair0 = u 
!     dupair = -1/r dupair/dr
      dupair0 = -dudr *ri
      if(periodic) then 
         if(rdis > rc0 .and. rdis < rc1) then 
            upair = upair + upair0 * switch
            dupair = dupair - upair0 * dswitchdr * ri + dupair0 * switch 
         else if(rdis < rc0) then
            upair = upair + upair0
            dupair = dupair + dupair0
         else if(rdis > rc1) then 
         endif
      else
         upair = upair + upair0
         dupair = dupair + dupair0
      endif

      end subroutine calc_pair_from_spline


      subroutine batch_xyz
      use common_data
      implicit none
      character(len = 4) :: name
      real(8) :: x,y,z
      integer :: iatom,moltype,natoms_in_mol,imol
      real(8), dimension(3) :: com
      real(8) :: totmass,rdis
      real(8) :: x1,y1,z1,x2,y2,z2,x3,y3,z3
      real(8), dimension(3,3) :: axes,axes0,axes0i,rotmat
      real(8) :: phi,theta,psi,stheta
      logical okflag
      integer ios,istep,ibin,n,nitems,i,nb,nterms
      character(len = 2048) :: cvec_buffer,buffer
      character(len = 2048), dimension(64) :: textlist      
      logical mindistance

      write(*,*) 'SUBROUTINE BATCH_XYZ'

      open(150,file = 'hist.dat')
      open(151,file = 'mindistance.dat')

      calc_fm = .true.
      use_spline = .false.
      model_coeff = 0.0
      rhist = 0.0d0 
      rhmax = 10.0d0 
      rhmin = 0.0d0 
      nhbins = 64

      open(78,file = 'model_forces.dat')
      write(78,*) nmodel_param,nmol,nbondlist
      do nb = 1,nbondlist
         nterms = model_nterms(nb)
         write(78,*) nb,model_list(nb,1),model_list(nb,2),(model_power(n,nb),n=1,nterms)
      end do 

      open(10,file = batch_file)
      istep = 0 
      do while(.true.)

!     call to read_input_file to read in the data and get the euler angles
         call read_input_file(ios)
         if(ios.ne.0) exit

         istep = istep + 1

         call energy()
         write(78,*) istep,utotal * pefac
         write(*,*) istep,utotal * pefac

         call calculate_rigid_body_forces()
         call print_forces_and_torques_and_derivs(78)
         call print_coordinates(140)
      end do 


!     write out the radial histograms. Useful for the force-matching stage.
      do n = 1,nbondlist
         write(150,*) n
         mindistance = .true.
         do ibin = 0,nhbins
            rdis = rhmin + (rhmax - rhmin)*dble(ibin)/dble(nhbins)
            write(150,*) rdis,rhist(n,ibin)
            if(mindistance.and.rhist(n,ibin).gt.10) then 
               write(151,18) model_list(n,1),model_list(n,2),rdis
 18            format(2i4,f10.5)
               mindistance = .false.
               endif
         end do 
         write(150,*)
      end do 

      stop
      end subroutine batch_xyz

      subroutine  batch
      use common_data
      implicit none
      character(len = 10) :: mintype
      integer n,ios
      real(8) :: dismin
      logical printrdf,too_close

      open(10,file = batch_file)
     
      open(140,file = 'batch.energy')
      open(141,file = 'batch.xyz')
      open(142,file = 'batch.rigid')

      printrdf = .true.
      open(77,file = 'rdf.dat')

      write(*,*) 'batch_file = ',batch_file

      do while(.true.)
         call read_input_file(ios)
         if(ios.ne.0) stop

         if(periodic) then
            call generate_cartesians_from_frac()
         else
            call generate_cartesians()
         endif
         mintol = 0.0d0 
         call minimize()
         call get_cell_mass
         if(periodic) then
            call calc_press_and_dens
            write(*,*) utotal * pefac,pressure_internal*fac(7),density*fac(8)
            write(140,*) utotal * pefac,pressure_internal*fac(7),density*fac(8)
            call flush(140)
         else
            write(*,*) utotal * pefac
            write(140,*) utotal * pefac
         endif
         call print_coordinates(141)
         call print_rigid_coordinates(142)
         call check_for_collisions(r_closest,too_close,dismin,printrdf)
         end do 
      stop

      end subroutine batch

      subroutine get_random_euler(phi,theta,psi,random)
      implicit none
      real(8) :: random
      real(8) :: x1,y1,z1,x2,y2,z2,x3,y3,z3
      real(8) :: dot,r1,r2
      real(8) :: rotmat(3,3)
      real(8) :: phi,theta,psi,stheta

      call make_rand_vec(random,x1,y1,z1)
      call make_rand_vec(random,x2,y2,z2)

!     project the r2 vector such that it's normal to r1, and in the 12 plane

      dot = x1*x2 + y1*y2 + z1*z2

      x2 = x2 - dot * x1
      y2 = y2 - dot * y1
      z2 = z2 - dot * z1

      r1 = dsqrt(x1**2 + y1**2 + z1**2)
      x1 = x1 / r1
      y1 = y1 / r1
      z1 = z1 / r1

      r2 = dsqrt(x2**2 + y2**2 + z2**2)
      x2 = x2 / r2
      y2 = y2 / r2
      z2 = z2 / r2

!     take cross product to get the components of r3
      x3 = y1*z2 - z1*y2
      y3 = z1*x2 - x1*z2
      z3 = x1*y2 - y1*x2

!     now that we have r1,r2,r3, generate the rotation matrix
      rotmat(1,1) = x1
      rotmat(1,2) = y1
      rotmat(1,3) = z1

      rotmat(2,1) = x2
      rotmat(2,2) = y2
      rotmat(2,3) = z2

      rotmat(3,1) = x3
      rotmat(3,2) = y3
      rotmat(3,3) = z3

      call get_euler_from_rotmat(rotmat,phi,theta,psi)

      end subroutine get_random_euler


      subroutine make_rand_vec(random,x,y,z)
      implicit none
      real(8) :: random
      real(8) :: pi
      real(8) :: theta
      real(8) :: x,y,z
      pi = 4.0d0 * atan(1.0d0)

      call random_number(random)
      theta = random * 2.0d0 * pi
      call random_number(random)
      z = (random-0.5d0)*2.0d0 

      x = dsqrt(1.0d0-z**2)*dcos(theta)
      y = dsqrt(1.0d0-z**2)*dsin(theta)

      end subroutine make_rand_vec

      subroutine get_optimum_lattice_vectors(new_lattice)
      use common_data
      implicit none
      real(8), dimension(3,3) :: cvec_new,rotmat,rotmat0,rotmat1
      real(8) :: phi,theta,psi,stheta
      real(8) :: comx,comy,comz
      logical new_lattice
      integer imol

!     firstly find the cartesians

      call generate_cartesians_from_frac()
      call find_new_cvec(cvec_new,rotmat,new_lattice)
      if(.not.new_lattice) return
      
      do imol = 1,nmol
         phi = mol_euler(1,imol)
         theta = mol_euler(2,imol)
         psi = mol_euler(3,imol) 
         
         call get_rotmat(phi,theta,psi,rotmat0)
         rotmat1 = matmul(rotmat,rotmat0)
         
         call get_euler_from_rotmat(rotmat1,phi,theta,psi)
         
         mol_euler(1,imol) = phi
         mol_euler(2,imol) = theta
         mol_euler(3,imol) = psi
         
         comx = mol_com(1,imol)
         comy = mol_com(2,imol)
         comz = mol_com(3,imol)
         
         mol_com(1,imol) = rotmat(1,1)*comx + rotmat(1,2)*comy + rotmat(1,3)*comz
         mol_com(2,imol) = rotmat(2,1)*comx + rotmat(2,2)*comy + rotmat(2,3)*comz
         mol_com(3,imol) = rotmat(3,1)*comx + rotmat(3,2)*comy + rotmat(3,3)*comz

      end do 
      cvec = cvec_new

!     finally, transform back into fractional coordinates.
      call get_frac_coords()

      end subroutine get_optimum_lattice_vectors

      subroutine find_new_cvec(cvec_new,rotmat,new_lattice)
      use common_data
      implicit none
      real(8), dimension(3,3) :: cvec_new,rotmat
      logical new_lattice,change
      real(8) :: angfac
      real(8) :: ax,ay,az,bx,by,bz,cx,cy,cz,amag,bmag,cmag,alpha,beta,gamma
      real(8) :: fx,fy,fz,gx,gy,gz,fmag,gmag,dot_ab,dot_bc,dot_ca
      real(8) :: xx,xy,xz,yx,yy,yz,zx,zy,zz
      real(8) :: axnorm,aynorm,aznorm,bxnorm,bynorm,bznorm
      integer nchange

      angfac = 180.0d0/pi

      ax = cvec(1,1)
      ay = cvec(2,1) 
      az = cvec(3,1)

      bx = cvec(1,2)
      by = cvec(2,2) 
      bz = cvec(3,2) 

      cx = cvec(1,3)
      cy = cvec(2,3) 
      cz = cvec(3,3)

      change = .true.
      nchange = 0
      do while(change)
         change = .false.
         call change_lattice_vector(ax,ay,az,bx,by,bz,change)
         call change_lattice_vector(ax,ay,az,-bx,-by,-bz,change)
         call change_lattice_vector(ax,ay,az,cx,cy,cz,change)
         call change_lattice_vector(ax,ay,az,-cx,-cy,-cz,change)
       
         call change_lattice_vector(bx,by,bz,cx,cy,cz,change)
         call change_lattice_vector(bx,by,bz,-cx,-cy,-cz,change)
         call change_lattice_vector(bx,by,bz,ax,ay,az,change)
         call change_lattice_vector(bx,by,bz,-ax,-ay,-az,change)
       
         call change_lattice_vector(cx,cy,cz,ax,ay,az,change)
         call change_lattice_vector(cx,cy,cz,-ax,-ay,-az,change)
         call change_lattice_vector(cx,cy,cz,bx,by,bz,change)
         call change_lattice_vector(cx,cy,cz,-bx,-by,-bz,change)
         if(change) nchange = nchange + 1
      end do 

      if(nchange.eq.0) then 
         new_lattice = .false.
         return
      else
         new_lattice = .true.
      endif

      amag = dsqrt(ax**2 + ay**2 + az**2)
      bmag = dsqrt(bx**2 + by**2 + bz**2)
      cmag = dsqrt(cx**2 + cy**2 + cz**2)

      dot_ab = (ax * bx + ay * by + az * bz)/(amag * bmag)
      dot_bc = (bx * cx + by * cy + bz * cz)/(bmag * cmag)
      dot_ca = (cx * ax + cy * ay + cz * az)/(cmag * amag)

      alpha = dacos(dot_bc) * angfac
      beta = dacos(dot_ca) * angfac
      gamma = dacos(dot_ab) * angfac

      call get_cartesian_cvec(amag,bmag,cmag,alpha,beta,gamma,cvec_new)

      axnorm = ax / amag
      aynorm = ay / amag
      aznorm = az / amag

      bxnorm = bx / bmag
      bynorm = by / bmag
      bznorm = bz / bmag

      fx = bxnorm - dot_ab * axnorm
      fy = bynorm - dot_ab * aynorm
      fz = bznorm - dot_ab * aznorm

      fmag = dsqrt(fx**2 + fy**2 + fz**2)
      
      xx = axnorm
      xy = aynorm
      xz = aznorm
      
      yx = fx / fmag
      yy = fy / fmag
      yz = fz / fmag

      gx = xy * yz - xz * yy
      gy = xz * yx - xx * yz
      gz = xx * yy - xy * yx

      gmag = dsqrt(gx**2 + gy**2 + gz**2)

      zx = gx / gmag
      zy = gy / gmag
      zz = gz / gmag

      rotmat(1,1) = xx
      rotmat(1,2) = xy
      rotmat(1,3) = xz

      rotmat(2,1) = yx
      rotmat(2,2) = yy
      rotmat(2,3) = yz

      rotmat(3,1) = zx
      rotmat(3,2) = zy
      rotmat(3,3) = zz

      end subroutine find_new_cvec

      subroutine change_lattice_vector(ax,ay,az,bx,by,bz,change)
      implicit none
      real(8) :: amag,bmag,amagnew
      real(8) :: axnew,aynew,aznew
      real(8) :: ax,ay,az,bx,by,bz
      logical change

      amag = dsqrt(ax**2 + ay**2 + az**2)
      bmag = dsqrt(bx**2 + by**2 + bz**2)

      axnew = ax + bx
      aynew = ay + by
      aznew = az + bz
      
      amagnew = dsqrt(axnew**2 + aynew**2 + aznew**2)

      if(amagnew.lt.amag) then
         ax = axnew
         ay = aynew
         az = aznew
         change = .true.
      endif

      end subroutine change_lattice_vector

      subroutine rotate_2tensor(rotmat,mat0,mat1)
      implicit none
      real(8), dimension(3,3) :: rotmat,mat0,mat1
      integer ivec,jvec,kvec,lvec
      real(8) :: sum

      mat1 = 0.0d0 
      do ivec = 1,3
         do jvec = 1,3
            sum = 0.0d0 
            do kvec = 1,3
               do lvec = 1,3
                  sum = sum + rotmat(ivec,kvec) * rotmat(jvec,lvec) * mat0(kvec,lvec) 
               end do 
            end do 
            mat1(ivec,jvec) = sum 
         end do 
      end do 
      
      end subroutine rotate_2tensor


      subroutine rotate_3tensor(rotmat,mat0,mat1)
      implicit none
      real(8), dimension(3,3) :: rotmat
      real(8), dimension(3,3,3) :: mat0,mat1
      integer ivec,jvec,kvec,lvec,mvec,nvec
      real(8) :: sum
      
      mat1 = 0.0d0 
      do ivec = 1,3
         do jvec = 1,3
            do kvec = 1,3
               sum = 0.0d0 
               do lvec = 1,3
                  do mvec = 1,3
                     do nvec = 1,3
                        sum = sum + rotmat(ivec,lvec) * rotmat(jvec,mvec) * rotmat(kvec,nvec) * mat0(lvec,mvec,nvec) 
                     end do 
                  end do 
               end do 
               mat1(ivec,jvec,kvec) = sum 
            end do 
         end do 
      end do 
      
      end subroutine rotate_3tensor

      subroutine print_multipoles
      use common_data
      implicit none
      integer moltype,natoms_in_mol,imol,iatom
      integer ivec,jvec,kvec,k

      real(8) :: charge_mol
      real(8), dimension(3) :: dipole_mol,solid1
      real(8), dimension(3,3) :: quad_mol
      real(8), dimension(3,3,3) :: oct_mol
      real(8), dimension(3) :: rvec
      real(8) :: q00_dma
      real(8) :: q10_dma,q11c_dma,q11s_dma
      real(8) :: qxx,qxy,qxz,qyy,qyz,qzz
      real(8) :: q20,q21c,q21s,q22c,q22s
      real(8) :: q20_dma,q21c_dma,q21s_dma,q22c_dma,q22s_dma
      real(8) :: q30,q31c,q31s,q32c,q32s,q33c,q33s
      real(8) :: q30_dma,q31c_dma,q31s_dma,q32c_dma,q32s_dma,q33c_dma,q33s_dma
      real(8) :: qxx_g,qxy_g,qxz_g,qyy_g,qyz_g,qzz_g,q20_g,q21c_g,q21s_g,q22c_g,q22s_g
      real(8) :: qxxx,qxxy,qxyy,qyyy,qxxz,qxyz,qyyz,qxzz,qyzz,qzzz
      real(8) :: qxxx_g,qxxy_g,qxyy_g,qyyy_g,qxxz_g,qxyz_g,qyyz_g,qxzz_g,qyzz_g,qzzz_g
      real(8) :: qxx2,qxy2,qxz2,qyy2,qyz2,qzz2
      real(8), dimension(5) :: solid2
      real(8), dimension(7) :: solid3
      real(8), dimension(0:10) :: conv,stone_conv

      write(*,*) '*****************************************************'
      write(*,*) '***   MOLECULAR MULTIPOLE MOMENTS IN &
     &GAUSSIAN AND GDMA FORMAT'
      write(*,*) '*****************************************************'

!     conversion factors for multipoles from e bohr^n to internal units
      conv(0) = 1.0d0 
      do k = 1,10
         conv(k) = conv(k-1) * fac(5)
      end do 

!     conversion factors between Stone's normalisation and ours
      stone_conv(0) = 1.0d0
      stone_conv(1) = 1.0d0
      stone_conv(2) = dsqrt(6.0d0) / 6.0d0
      stone_conv(3) = dsqrt(10.0d0)/ 30.0d0

      if(periodic) then
         call generate_cartesians_from_frac()
      else
         call generate_cartesians()
      endif

      do imol = 1,nmol
         charge_mol = 0.0d0 
         dipole_mol = 0.0d0 
         quad_mol = 0.0d0
         oct_mol = 0.0d0 

         moltype = mol_type(imol)
         natoms_in_mol = mol_natoms(moltype)

         do iatom = 1,natoms_in_mol
            charge_mol = charge_mol + lab_charge(iatom,imol)
            do ivec = 1,3
               rvec(ivec) = lab_relative_coord(ivec,iatom,imol)
            end do 
            do ivec = 1,3
               dipole_mol(ivec) = dipole_mol(ivec) & 
     &              + lab_charge(iatom,imol) * rvec(ivec) & 
     &              + lab_dipole(ivec,iatom,imol)  
               do jvec = 1,3
                  quad_mol(ivec,jvec) = quad_mol(ivec,jvec) &
     &              + 0.5d0 * lab_charge(iatom,imol) * rvec(ivec) * rvec(jvec) & 
     &              + 0.5d0 * (lab_dipole(ivec,iatom,imol) * rvec(jvec) + lab_dipole(jvec,iatom,imol) * rvec(ivec)) & 
     &              + lab_quad(ivec,jvec,iatom,imol) 
                  do kvec = 1,3
                     oct_mol(ivec,jvec,kvec) = oct_mol(ivec,jvec,kvec) &
     &                    + (1.0d0 / 6.0d0) * lab_charge(iatom,imol) * rvec(ivec) * rvec(jvec) * rvec(kvec) & 
     &                    + (1.0d0 / 6.0d0) * (lab_dipole(ivec,iatom,imol) * rvec(jvec) * rvec(kvec) &
     &                    +                    lab_dipole(jvec,iatom,imol) * rvec(kvec) * rvec(ivec) &
     &                    +                    lab_dipole(kvec,iatom,imol) * rvec(ivec) * rvec(jvec)) &
     &                    + (2.0d0 / 6.0d0) * (lab_quad(ivec,jvec,iatom,imol) * rvec(kvec) & 
     &                    +                    lab_quad(jvec,kvec,iatom,imol) * rvec(ivec) & 
     &                    +                    lab_quad(kvec,ivec,iatom,imol) * rvec(jvec)) & 
     &                    + lab_oct(ivec,jvec,kvec,iatom,imol) 
                  end do 
               end do 
            end do 
         end do 

         qxx = quad_mol(1,1)
         qxy = quad_mol(1,2)
         qxz = quad_mol(1,3)
         qyy = quad_mol(2,2)
         qyz = quad_mol(2,3)
         qzz = quad_mol(3,3)

!     remove traces 
         call convert_quad_to_spherical(qxx,qxy,qxz,qyy,qyz,qzz,solid2)
         call convert_quad_to_cartesian(solid2,qxx,qxy,qxz,qyy,qyz,qzz)

         qxx_g = 2.0d0 * qxx * fac(10)
         qxy_g = 2.0d0 * qxy * fac(10) 
         qxz_g = 2.0d0 * qxz * fac(10) 
         qyy_g = 2.0d0 * qyy * fac(10)
         qyz_g = 2.0d0 * qyz * fac(10)
         qzz_g = 2.0d0 * qzz * fac(10) 

         qxxx = oct_mol(1,1,1)
         qxxy = oct_mol(1,1,2)
         qxyy = oct_mol(1,2,2)
         qyyy = oct_mol(2,2,2)
         qxxz = oct_mol(1,1,3)
         qxyz = oct_mol(1,2,3)
         qyyz = oct_mol(2,2,3)
         qxzz = oct_mol(1,3,3)
         qyzz = oct_mol(2,3,3)
         qzzz = oct_mol(3,3,3)

!     remove traces 
         call  convert_oct_to_spherical(qxxx,qxxy,qxyy,qyyy,qxxz,qxyz,qyyz,qxzz,qyzz,qzzz,solid3)
         call convert_oct_to_cartesian(solid3,qxxx,qxxy,qxyy,qyyy,qxxz,qxyz,qyyz,qxzz,qyzz,qzzz)

         qxxx_g = 6.0d0 * qxxx * fac(10)
         qxxy_g = 6.0d0 * qxxy * fac(10)
         qxyy_g = 6.0d0 * qxyy * fac(10)
         qyyy_g = 6.0d0 * qyyy * fac(10)
         qxxz_g = 6.0d0 * qxxz * fac(10)
         qxyz_g = 6.0d0 * qxyz * fac(10)
         qyyz_g = 6.0d0 * qyyz * fac(10)
         qxzz_g = 6.0d0 * qxzz * fac(10)
         qyzz_g = 6.0d0 * qyzz * fac(10)
         qzzz_g = 6.0d0 * qzzz * fac(10)
         
         write(*,*) 'molecule ',imol
         write(*,*) 'Gaussian format'
         write(*,19) 'C = ',charge_mol
         write(*,19) 'Dx = ',dipole_mol(1)*fac(10),'Dy = ',dipole_mol(2)*fac(10),'Dz = ',dipole_mol(3)*fac(10)
         write(*,19) 'Qxx = ',qxx_g,'Qyy = ',qyy_g,'Qzz = ',qzz_g
         write(*,19) 'Qxy = ',qxy_g,'Qxz = ',qxz_g,'Qyz = ',qyz_g
         write(*,19) 'Qxxx = ',qxxx_g,'Qyyy = ',qyyy_g,'Qzzz = ',qzzz_g
         write(*,19) 'Qxxy = ',qxxy_g,'Qxxz = ',qxxz_g,'Qxzz = ',qxzz_g
         write(*,19) 'Qyyz = ',qxxy_g,'Qyz = ',qxyz_g

 19      format(a8,f14.6,a8,f14.6,a8,f14.6)

         q00_dma = charge_mol

         solid1 = dipole_mol

         solid1  = solid1 / (conv(1) * stone_conv(1))
         solid2  = solid2 / (conv(2) * stone_conv(2))
         solid3  = solid3 / (conv(3) * stone_conv(3))

         q11c_dma = solid1(1)
         q11s_dma = solid1(2)
         q10_dma =  solid1(3)

         q20_dma =  solid2(1)
         q21c_dma = solid2(2)
         q21s_dma = solid2(3)
         q22c_dma = solid2(4)
         q22s_dma = solid2(5)

         q30_dma = solid3(1) 
         q31c_dma = solid3(2)
         q31s_dma = solid3(3)
         q32c_dma = solid3(4)
         q32s_dma = solid3(5)
         q33c_dma = solid3(6)
         q33s_dma = solid3(7)

         write(*,*) 'GDMA format'

         write(*,19) 'Q00 = ',q00_dma
         write(*,19) 'Q10 = ',q10_dma,'Q11c = ',q11c_dma,'Q11s = ',q11s_dma
         write(*,19) 'Q20 = ',q20_dma,'Q21c = ',q21c_dma,'Q21s = ',q21s_dma
         write(*,19) 'Q22c = ',q22c_dma,'Q22s = ',q22s_dma
         write(*,19) 'Q30 = ',q30_dma,'Q31c = ',q31c_dma,'Q31s = ',q31s_dma
         write(*,19) 'Q32c = ',q32c_dma,'Q32s = ',q32s_dma,'Q33c = ',q33c_dma
         write(*,19) 'Q33s = ',q33s_dma

         write(*,*) 
      end do 


      end subroutine print_multipoles

      subroutine convert_quad_to_spherical(xx,xy,xz,yy,yz,zz,solid2)
      implicit none
      real(8) :: xx,xy,xz,yy,yz,zz,rr
      real(8), dimension(5) :: solid2
      real(8) :: fac1 = dsqrt(6.0d0)/6.0d0, fac2 = dsqrt(2.0d0), fac3 = dsqrt(2.0d0)/2.0d0

!     converts quadrupoles to spherical form

      rr = xx + yy + zz

!     20
      solid2(1) = fac1 * (3.0d0 * zz - rr)
!     21c
      solid2(2) = fac2 * xz
!     21s
      solid2(3) = fac2 * yz
!     22c
      solid2(4) = fac3 * (xx - yy)
!     22s
      solid2(5) = fac2 * xy

      end subroutine convert_quad_to_spherical

      subroutine convert_quad_to_cartesian(solid2,xx,xy,xz,yy,yz,zz)
      implicit none
      real(8) :: xx,xy,xz,yy,yz,zz
      real(8), dimension(5) :: solid2
      real(8) :: fac1 = dsqrt(6.0d0)/6.0d0, fac2 = dsqrt(2.0d0)/2.0d0, fac3 = dsqrt(6.0d0)/3.0d0

!     gives traceless quadrupole moments in Cartesian form

      xx = -fac1 * solid2(1) + fac2 * solid2(4)
      yy = -fac1 * solid2(1) - fac2 * solid2(4)
      zz = solid2(1) *fac3

      xy = fac2 * solid2(5)
      xz = fac2 * solid2(2)
      yz = fac2 * solid2(3)

      end subroutine convert_quad_to_cartesian


      subroutine convert_oct_to_cartesian(solid3,&
     &           xxx,xxy,xyy,yyy,xxz,xyz,yyz,xzz,yzz,zzz)
      implicit none
      real(8) :: xxx,xxy,xyy,yyy,xxz,xyz,yyz,xzz,yzz,zzz
      real(8), dimension(7) :: solid3
      real(8) :: fac1 = dsqrt(15.0d0)/10.0d0, fac2 = dsqrt(15.0d0)/30.0d0, fac3 = dsqrt(10.0d0)/10.0d0 & 
     &     ,fac4 = dsqrt(6.0d0)/6.0d0, fac5 = dsqrt(15.0d0)/15.0d0, fac6 = dsqrt(2.0d0 / 5.0d0)

!     gives traceless octopole moments in Cartesian form

      xxx = 0.5d0 * solid3(6) - fac1 * solid3(2)
      xxy = 0.5d0 * solid3(7) - fac2 * solid3(3)
      xyy = -0.5d0 * solid3(6) - fac2  * solid3(2)
      yyy = - 0.5d0 * solid3(7) - fac1 * solid3(3)
      xxz = fac4 * solid3(4) - fac3 * solid3(1)
      xyz = fac4 * solid3(5)
      yyz = -fac4 * solid3(4) - fac3 * solid3(1)
      xzz = 2.0d0 * fac5 * solid3(2)
      yzz = 2.0d0 * fac5 * solid3(3)
      zzz = fac6 * solid3(1)

      end subroutine convert_oct_to_cartesian

      subroutine convert_oct_to_spherical(xxx,xxy,xyy,yyy,xxz,xyz,yyz,xzz,yzz,zzz,&
     &     solid3)
      implicit none
      real(8) :: xxx,xxy,xyy,yyy,xxz,xyz,yyz,xzz,yzz,zzz,xrr,yrr,zrr
      real(8), dimension(7) :: solid3
      real(8) :: fac1 = dsqrt(10.0d0)/10.0d0, fac2 = dsqrt(15.0d0)/10.0d0, fac3 = dsqrt(3.0d0/2.0d0) & 
     &     , fac4 = dsqrt(6.0d0)

!     converts octopoles to spherical form
!     factor of 6 because oct_internal = (1/3!) sum_i C_i r_i,alpha r_i,beta r_i,gamma

      xrr = xxx + xyy + xzz
      yrr = xxy + yyy + yzz
      zrr = xxz + yyz + zzz

!     30
      solid3(1) = fac1 * (5 * zzz - 3.0d0 * zrr)
!     31c
      solid3(2) = fac2 * (5.0d0 * xzz - xrr)
!     31s
      solid3(3) = fac2 * (5.0d0 * yzz - yrr)
!     32c
      solid3(4) = fac3 * (xxz - yyz)
!     32s
      solid3(5) = fac4 * xyz
!     33c
      solid3(6) = 0.5d0 * (xxx - 3.0d0 * xyy)
!     33s
      solid3(7) = 0.5d0 * (3.0d0 * xxy - yyy)
      end subroutine convert_oct_to_spherical


      subroutine convert_tens2_to_cartesian(solid2,cart2)
      implicit none
      real(8), dimension(5) :: solid2
      real(8), dimension(3,3) :: cart2
      real(8) :: xx,xy,xz,yy,yz,zz
      
      call convert_quad_to_cartesian(solid2,xx,xy,xz,yy,yz,zz)

      cart2(1,1) = xx
      cart2(1,2) = xy
      cart2(1,3) = xz
      cart2(2,2) = yy
      cart2(2,3) = yz
      cart2(3,3) = zz

      cart2(2,1) = cart2(1,2)
      cart2(3,1) = cart2(1,3)
      cart2(3,2) = cart2(2,3)

      end subroutine convert_tens2_to_cartesian

      subroutine convert_tens3_to_cartesian(solid3,cart3)
      implicit none
      real(8), dimension(7) :: solid3
      real(8), dimension(3,3,3) :: cart3
      real(8) :: xxx,xxy,xyy,yyy,xxz,xyz,yyz,xzz,yzz,zzz
      integer ivec,jvec,kvec

      call convert_oct_to_cartesian(solid3,xxx,xxy,xyy,yyy,xxz,xyz,yyz,xzz,yzz,zzz)

      cart3(1,1,1) = xxx
      cart3(1,1,2) = xxy
      cart3(1,2,2) = xyy
      cart3(2,2,2) = yyy
      cart3(1,1,3) = xxz
      cart3(1,2,3) = xyz
      cart3(2,2,3) = yyz
      cart3(1,3,3) = xzz
      cart3(2,3,3) = yzz
      cart3(3,3,3) = zzz

      do kvec = 1,3
         do jvec = 1,kvec
            do ivec = 1,jvec
               cart3(ivec,kvec,jvec) = cart3(ivec,jvec,kvec)
               cart3(jvec,ivec,kvec) = cart3(ivec,jvec,kvec)
               cart3(jvec,kvec,ivec) = cart3(ivec,jvec,kvec)
               cart3(kvec,ivec,jvec) = cart3(ivec,jvec,kvec)
               cart3(kvec,jvec,ivec) = cart3(ivec,jvec,kvec)
            end do 
         end do 
      end do 


      end subroutine convert_tens3_to_cartesian

      subroutine convert_tens2_1_to_cartesian(solid2_1,cart3)
      implicit none
      real(8) :: xxx,xxy,xyy,yyy,xxz,xyz,yyz,xzz,yzz,zzz
      real(8), dimension(5) :: solid2
      real(8), dimension(7) :: solid3
      real(8), dimension(5,3) :: solid2_1
      real(8), dimension(3,3,3) :: cart3
      real(8) :: xx,xy,xz,yy,yz,zz
      integer ii,jj,ivec,jvec,kvec

!     doesn't assume that the resulting cartesian tensor has symmetry

      cart3 = 0.0d0 

      do ii = 1,3
         do jj = 1,5
            solid2(jj) = solid2_1(jj,ii)
         end do 

         call convert_quad_to_cartesian(solid2,xx,xy,xz,yy,yz,zz)

         cart3(ii,1,1) = xx
         cart3(ii,1,2) = xy
         cart3(ii,1,3) = xz
         cart3(ii,2,2) = yy
         cart3(ii,2,3) = yz
         cart3(ii,3,3) = zz
      end do 


      do ivec = 1,3
         do kvec = 1,3
            do jvec = 1,kvec
               cart3(ivec,kvec,jvec) = cart3(ivec,jvec,kvec)
            end do 
         end do
      end do 



      end subroutine convert_tens2_1_to_cartesian


      subroutine get_euler_from_rotmat(rotmat,phi,theta,psi)
      implicit none
      real(8), dimension(3,3) :: rotmat
      real(8) :: phi,theta,psi,stheta,spsi

      phi = datan2(rotmat(3,1),-rotmat(3,2))
      psi = datan2(rotmat(1,3),rotmat(2,3))

      spsi = dsin(psi)

      if(dabs(spsi).ge.1.d-5) then 
         stheta = rotmat(1,3) / spsi
      else
         stheta = rotmat(2,3) / dcos(psi)
      endif
      theta = datan2(stheta,rotmat(3,3))

      end subroutine get_euler_from_rotmat

      subroutine symmetrise_fields()
      use common_data
      implicit none
      real(8) :: qxx,qxy,qxz,qyy,qyz,qzz
      real(8) :: qxxx,qxxy,qxyy,qyyy,qxxz,qxyz,qyyz,qxzz,qyzz,qzzz
      real(8) :: q20,q21c,q21s,q22c,q22s
      real(8) :: q30,q31c,q31s,q32c,q32s,q33c,q33s
      real(8), dimension(5) :: solid2
      real(8), dimension(7) :: solid3
      real(8), dimension(5,3) :: solid2_1
      real(8) :: xxx,xxy,xyy,yyy,xxz,xyz,yyz,xzz,yzz,zzz
      real(8) :: xxx2_1,xxy2_1,xyy2_1,yyy2_1,xxz2_1,xyz2_1,yyz2_1,xzz2_1,yzz2_1,zzz2_1
      real(8), dimension(3,3,3) :: cart3
      real(8), dimension(3,3) :: cart2
      integer :: imol,iatom,natoms_in_mol,moltype
      integer ivec,jvec,kvec,ii,jj

      
!     makes sure that the fields have the right symmetry, and also removes their traces

      do imol = 1,nmol
         moltype = mol_type(imol)
         natoms_in_mol = mol_natoms(moltype)
         do iatom = 1,natoms_in_mol

            cfield(iatom,imol) = cfield(iatom,imol) + fieldtens0(iatom,imol)
            
            do ii = 1,3
               dfield(ii,iatom,imol) = dfield(ii,iatom,imol) + fieldtens1(ii,iatom,imol)
            end do 

            do ii = 1,5
               solid2(ii) = fieldtens2(ii,iatom,imol)
            end do 
            call convert_tens2_to_cartesian(solid2,cart2)

            do ivec = 1,3
               do jvec = 1,3
                  qfield(ivec,jvec,iatom,imol) = qfield(ivec,jvec,iatom,imol) + cart2(ivec,jvec) & 
     &                 + fieldtens1_1(ivec,jvec,iatom,imol)
               end do 
            end do 

            do ii = 1,5
               do jj = 1,3
                  solid2_1(ii,jj) = fieldtens2_1(ii,jj,iatom,imol)
               end do 
            end do 
            call convert_tens2_1_to_cartesian(solid2_1,cart3)
            
            do ivec = 1,3
               do jvec = 1,3
                  do kvec = 1,3
                     ofield(ivec,jvec,kvec,iatom,imol) = ofield(ivec,jvec,kvec,iatom,imol) + cart3(ivec,jvec,kvec)
                  end do 
               end do 
            end do 

            do ii = 1,7
               solid3(ii) = fieldtens3(ii,iatom,imol)
            end do 

            call convert_tens3_to_cartesian(solid3,cart3)

            do ivec = 1,3
               do jvec = 1,3
                  do kvec = 1,3
                     ofield(ivec,jvec,kvec,iatom,imol) = ofield(ivec,jvec,kvec,iatom,imol) + cart3(ivec,jvec,kvec)
                  end do 
               end do 
            end do 


            
!     symmetrise the fields
            do ivec = 1,3
               if(max_rank.ge.2) then 
                  do jvec = 1,ivec
                     qfield(ivec,jvec,iatom,imol) = (qfield(ivec,jvec,iatom,imol) + qfield(jvec,ivec,iatom,imol))/2.0d0 
                     qfield(jvec,ivec,iatom,imol) = qfield(ivec,jvec,iatom,imol)
                     if(max_rank.ge.3) then 
                        do kvec = 1,jvec
                           ofield(ivec,jvec,kvec,iatom,imol) = &
     &                           (ofield(ivec,jvec,kvec,iatom,imol) + ofield(ivec,kvec,jvec,iatom,imol) &
     &                          + ofield(jvec,ivec,kvec,iatom,imol) + ofield(jvec,kvec,ivec,iatom,imol) &
     &                          + ofield(kvec,ivec,jvec,iatom,imol) + ofield(kvec,jvec,ivec,iatom,imol))/6.0d0 
                           ofield(ivec,kvec,jvec,iatom,imol) = ofield(ivec,jvec,kvec,iatom,imol)
                           ofield(jvec,ivec,kvec,iatom,imol) = ofield(ivec,jvec,kvec,iatom,imol)
                           ofield(jvec,kvec,ivec,iatom,imol) = ofield(ivec,jvec,kvec,iatom,imol)
                           ofield(kvec,ivec,jvec,iatom,imol) = ofield(ivec,jvec,kvec,iatom,imol)
                           ofield(kvec,jvec,ivec,iatom,imol) = ofield(ivec,jvec,kvec,iatom,imol)
                        end do 
                     endif
                  end do
               endif
            end do 
            
!     remove the traces by converting into spherical harmonics and then back to cartesians
!     doesn't affect the energies, but still nice to have the traceless fields. Also pack the fields 
!     back into fieldtens for later calculation of the torques

            if(max_rank.ge.1) then 
               do ii = 1,3
                  fieldtens1(ii,iatom,imol) = dfield(ii,iatom,imol)
               end do 
            endif
            if(max_rank.ge.2) then 

               qxx = qfield(1,1,iatom,imol)
               qxy = qfield(1,2,iatom,imol)
               qxz = qfield(1,3,iatom,imol)
               qyy = qfield(2,2,iatom,imol)
               qyz = qfield(2,3,iatom,imol)
               qzz = qfield(3,3,iatom,imol)
               
               call convert_quad_to_spherical(qxx,qxy,qxz,qyy,qyz,qzz,solid2)
               call convert_quad_to_cartesian(solid2,qxx,qxy,qxz,qyy,qyz,qzz)
               
               qfield(1,1,iatom,imol) = qxx
               qfield(1,2,iatom,imol) = qxy
               qfield(1,3,iatom,imol) = qxz
               qfield(2,2,iatom,imol) = qyy
               qfield(2,3,iatom,imol) = qyz
               qfield(3,3,iatom,imol) = qzz
               
               do jvec = 1,3
                  do ivec = 1,jvec
                     qfield(jvec,ivec,iatom,imol) = qfield(ivec,jvec,iatom,imol)
                  end do 
               end do 

               do ivec = 1,3
                  do jvec = 1,3
                     fieldtens1_1(ivec,jvec,iatom,imol) = qfield(ivec,jvec,iatom,imol)
                  end do 
               end do 
               
            endif
            
            if(max_rank.ge.3) then 
               qxxx = ofield(1,1,1,iatom,imol)
               qxxy = ofield(1,1,2,iatom,imol)
               qxyy = ofield(1,2,2,iatom,imol)
               qyyy = ofield(2,2,2,iatom,imol)
               qxxz = ofield(1,1,3,iatom,imol)
               qxyz = ofield(1,2,3,iatom,imol)
               qyyz = ofield(2,2,3,iatom,imol)
               qxzz = ofield(1,3,3,iatom,imol)
               qyzz = ofield(2,3,3,iatom,imol)
               qzzz = ofield(3,3,3,iatom,imol)
               
               call  convert_oct_to_spherical(qxxx,qxxy,qxyy,qyyy,qxxz,qxyz,qyyz,qxzz,qyzz,qzzz,solid3)
               call convert_oct_to_cartesian(solid3,qxxx,qxxy,qxyy,qyyy,qxxz,qxyz,qyyz,qxzz,qyzz,qzzz)

               call make_cart3(qxxx,qxxy,qxyy,qyyy,qxxz,qxyz,qyyz,qxzz,qyzz,qzzz,cart3)

               do ivec = 1,3
                  do jvec = 1,3
                     do kvec = 1,3
                        ofield(ivec,jvec,kvec,iatom,imol) = cart3(ivec,jvec,kvec)
                     end do 
                  end do 
               end do 
               
               call get_spherical3(qxxx,qxxy,qxyy,qyyy,qxxz,qxyz,qyyz,qxzz,qyzz,qzzz,solid3,solid2_1)

               do ii = 1,3
                  do jj = 1,5
                     fieldtens2_1(jj,ii,iatom,imol) = solid2_1(jj,ii)
                  end do 
               end do 

            endif
         end do 
      end do 

      end subroutine symmetrise_fields


      subroutine get_spherical3(xxx,xxy,xyy,yyy,xxz,xyz,yyz,xzz,yzz,zzz,solid3,solid2_1)
      implicit none
      real(8), dimension(5) :: solid2
      real(8), dimension(7) :: solid3
      real(8), dimension(5,3) :: solid2_1
      real(8) :: xxx,xxy,xyy,yyy,xxz,xyz,yyz,xzz,yzz,zzz
      real(8) :: xx,xy,xz,yy,yz,zz
      integer ii,jj 

!     gets the spherical tensors 3,2_1

      call  convert_oct_to_spherical(xxx,xxy,xyy,yyy,xxz,xyz,yyz,xzz,yzz,zzz,solid3)
            
      xx = xxx
      xy = xxy
      xz = xxz
      yy = xyy
      yz = xyz
      zz = xzz
      
      call convert_quad_to_spherical(xx,xy,xz,yy,yz,zz,solid2)
      
      do ii = 1,5
         solid2_1(ii,1) = solid2(ii)
      end do 
      
      xx = xxy
      xy = xyy
      xz = xyz
      yy = yyy
      yz = yyz
      zz = yzz
      
      call convert_quad_to_spherical(xx,xy,xz,yy,yz,zz,solid2)
      
      do ii = 1,5
         solid2_1(ii,2) = solid2(ii)
      end do 
      
      xx = xxz
      xy = xyz
      xz = xzz
      yy = yyz
      yz = yzz
      zz = zzz
      
      call convert_quad_to_spherical(xx,xy,xz,yy,yz,zz,solid2)

      do ii = 1,5
         solid2_1(ii,3) = solid2(ii)
      end do 

      end subroutine get_spherical3

      subroutine make_cart2(xx,xy,xz,yy,yz,zz,cart2)
      real(8) :: xx,xy,xz,yy,yz,zz
      real(8), dimension(3,3) :: cart2
      cart2(1,1) = xx
      cart2(2,2) = yy
      cart2(3,3) = zz

      cart2(1,2) = xy
      cart2(1,3) = xz
      cart2(2,3) = yz

      cart2(2,1) = xy
      cart2(3,1) = xz
      cart2(3,2) = yz
      end subroutine make_cart2


      subroutine make_cart3(xxx,xxy,xyy,yyy,xxz,xyz,yyz,xzz,yzz,zzz,cart3)
      real(8) :: xxx,xxy,xyy,yyy,xxz,xyz,yyz,xzz,yzz,zzz
      real(8), dimension(3,3,3) :: cart3

      cart3(1,1,1) = xxx
      cart3(1,1,2) = xxy
      cart3(1,2,2) = xyy
      cart3(2,2,2) = yyy
      cart3(1,1,3) = xxz
      cart3(1,2,3) = xyz
      cart3(2,2,3) = yyz
      cart3(1,3,3) = xzz
      cart3(2,3,3) = yzz
      cart3(3,3,3) = zzz
      
      do kvec = 1,3
         do jvec = 1,kvec
            do ivec = 1,jvec

               cart3(ivec,kvec,jvec) = cart3(ivec,jvec,kvec)
               cart3(jvec,ivec,kvec) = cart3(ivec,jvec,kvec)
               cart3(jvec,kvec,ivec) = cart3(ivec,jvec,kvec)
               cart3(kvec,ivec,jvec) = cart3(ivec,jvec,kvec)
               cart3(kvec,jvec,ivec) = cart3(ivec,jvec,kvec)

            end do 
         end do 
      end do 

      end subroutine make_cart3


      subroutine CALC_CONFORMER_ENERGY
      use common_data
      implicit none
      integer conformer,imoltype,imol

      do imol = 1,nmol
         imoltype = mol_type(imol)
         conformer = mol_conformer(imol)
         uconformer = uconformer + mol_conformer_energy(conformer,imoltype)
      end do
      end subroutine CALC_CONFORMER_ENERGY

      subroutine set_rdampcut(radius)
      use common_data
      implicit none
      real(8) :: radius
      rdamp_cutoff = radius
!     flag indicates whether the rdamp_cutoff variable was set by the user
!     or whether the default value should be used.
      rdamp_userset = .true.
      end subroutine set_rdampcut

      subroutine assign_multipoles(nmult,iatom,conformer,moltype,multipole_vec)
      use common_data
      implicit none
      integer nmult,iatom,conformer,moltype
      real(8), dimension(16) :: multipole_vec
      real(8) :: qxx,qxy,qxz,qyy,qyz,qzz
      real(8) :: qxxx,qxxy,qxyy,qyyy,qxxz,qxyz,qyyz,qxzz,qyzz,qzzz
      real(8), dimension(5) :: solid2
      real(8), dimension(7) :: solid3
      real(8), dimension(0:10) :: conv,stone_conv
      real(8) :: charge,tot_charge
      real(8), dimension(3,3) :: cart2
      real(8), dimension(3,3,3) :: cart3
      integer ivec,jvec,kvec
      integer k,n,rank
      real(8) :: dx,dy,dz

!     conversion factors between Stone's normalisation and ours
      stone_conv(0) = 1.0d0 
      stone_conv(1) = 1.0d0 
      stone_conv(2) = dsqrt(6.0d0) / 6.0d0
      stone_conv(3) = dsqrt(10.0d0)/ 30.0d0 

      tot_charge = 0.d0 

!     conversion factors for multipoles from e Bohr^n to internal units
      
      conv(0) = 1.0d0 
      do k = 1,10
         conv(k) = conv(k-1) * fac(5)
      end do 

!     determine the rank
      do n = nmult,1,-1
         if(multipole_vec(n).ne.0) exit
      end do 

      if(n.gt.0) rank = 0
      if(n.gt.1) rank = 1
      if(n.gt.4) rank = 2
      if(n.gt.9) rank = 3

      if(rank.ge.0) then 
         if(max_rank.lt.0) max_rank = 0
         if(site_rank(iatom,moltype).lt.0) site_rank(iatom,moltype) = 0 
         charge = multipole_vec(1)
         site_charge0(iatom,conformer,moltype) = charge
         tot_charge = tot_charge + charge
      endif
      
      if(rank.ge.1) then 
         if(max_rank.lt.1) max_rank = 1
         if(site_rank(iatom,moltype).lt.1) site_rank(iatom,moltype) = 1

         dx = multipole_vec(2)
         dy = multipole_vec(3)
         dz = multipole_vec(4)
         
!     convert to internal units
         site_dipole0(1,iatom,conformer,moltype) = conv(1) * dx 
         site_dipole0(2,iatom,conformer,moltype) = conv(1) * dy 
         site_dipole0(3,iatom,conformer,moltype) = conv(1) * dz 
      else
         site_dipole0(:,iatom,conformer,moltype) = 0.0d0 
      endif

      if(rank.ge.2) then 
!     read in the quadrupoles
         if(max_rank.lt.2) max_rank = 2
         if(site_rank(iatom,moltype).lt.2) site_rank(iatom,moltype) = 2

         solid2(1) = multipole_vec(5)
         solid2(2) = multipole_vec(6)
         solid2(3) = multipole_vec(7)
         solid2(4) = multipole_vec(8)
         solid2(5) = multipole_vec(9)
         
         call convert_quad_to_cartesian(solid2,qxx,qxy,qxz,qyy,qyz,qzz)
         
         call make_cart2(qxx,qxy,qxz,qyy,qyz,qzz,cart2)
         
         cart2 = cart2 * conv(2) * stone_conv(2) 
         
         do ivec = 1,3
            do jvec = 1,3
               site_quad0(ivec,jvec,iatom,conformer,moltype) = cart2(ivec,jvec)
            end do 
         end do 
         
      else
         site_quad0(:,:,iatom,conformer,moltype) = 0.0d0 
      endif

      if(rank.ge.3) then 
!     read in the octopoles
         if(max_rank.lt.3) max_rank = 3
         if(site_rank(iatom,moltype).lt.3) site_rank(iatom,moltype) = 3
         
         solid3(1) = multipole_vec(10)
         solid3(2) = multipole_vec(11)
         solid3(3) = multipole_vec(12)
         solid3(4) = multipole_vec(13)
         solid3(5) = multipole_vec(14)       
         solid3(6) = multipole_vec(15)       
         solid3(7) = multipole_vec(16)    
         
         call convert_oct_to_cartesian(solid3,qxxx,qxxy,qxyy,qyyy,qxxz,qxyz,qyyz,qxzz,qyzz,qzzz)
         
         call make_cart3(qxxx,qxxy,qxyy,qyyy,qxxz,qxyz,qyyz,qxzz,qyzz,qzzz,cart3)
         
         cart3 = cart3 * conv(3) * stone_conv(3) 
         
         do ivec = 1,3 
            do jvec = 1,3
               do kvec = 1,3
                  site_oct0(ivec,jvec,kvec,iatom,conformer,moltype) = cart3(ivec,jvec,kvec)
               end do 
            end do 
         end do 
         
      else
         site_oct0(:,:,:,iatom,conformer,moltype) = 0.0d0 
      endif
      
      end subroutine assign_multipoles


      subroutine get_lab_spherical2(iatom,imol,qxx,qxy,qxz,qyy,qyz,qzz)
      use common_data
      implicit none
      real(8) :: qxx,qxy,qxz,qyy,qyz,qzz
      real(8), dimension(5) :: solid2
      integer :: ii
      integer iatom,imol

      call convert_quad_to_spherical(qxx,qxy,qxz,qyy,qyz,qzz,solid2)
      
      do ii = 1,5
         lab_spherical2(ii,iatom,imol) = solid2(ii)
      end do 

      end subroutine get_lab_spherical2


      subroutine get_lab_spherical3(iatom,imol,qxxx,qxxy,qxyy,qyyy,qxxz,qxyz,qyyz,qxzz,qyzz,qzzz)
      use common_data
      implicit none
      real(8) :: qxxx,qxxy,qxyy,qyyy,qxxz,qxyz,qyyz,qxzz,qyzz,qzzz
      real(8), dimension(7) :: solid3
      real(8), dimension(5,3) :: solid2_1
      integer :: ii,jj
      integer iatom,imol

      call get_spherical3(qxxx,qxxy,qxyy,qyyy,qxxz,qxyz,qyyz,qxzz,qyzz,qzzz,solid3,solid2_1)
      do ii = 1,7
         lab_spherical3(ii,iatom,imol) = solid3(ii)
      end do 
      do ii = 1,5
         do jj = 1,3
            lab_spherical2_1(ii,jj,iatom,imol) = solid2_1(ii,jj)
         end do 
      end do 
      
      end subroutine get_lab_spherical3

      subroutine print_data
!     dummy subroutine for nr_mod call-back during minimization. But can be used to 
!     print results during a relaxation.
      end subroutine print_data

