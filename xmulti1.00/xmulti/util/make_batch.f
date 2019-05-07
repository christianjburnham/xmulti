      real(8) :: x,y,z
      real(8), dimension(16) :: multipole_vec
      character(len = 3) :: name
      character(len = 20) :: molecule_name
      open(10,file = 'pf.xyz')
      open(12,file = 'multipoles.dat')
      molecule_name = 'mannitol'
      natoms_in_mol = 26

      amag = 9.32046329401  
      bmag = 8.41214155705	 
      cmag = 15.35747855488   
      alpha = 94.97262   
      beta = 95.90939   
      gamma = 96.2887

      read(12,*) 
      do while(.true.)
         read(10,*,iostat = ios) natoms
         nmol = natoms / natoms_in_mol
         if(ios.lt.0) exit
         read(12,*,iostat = ios) 
         read(10,*)
         write(*,*) nmol
         write(*,*) 'ANGLE',amag,bmag,cmag,alpha,beta,gamma
         write(*,*) 'XYZ MULTIPOLE'
         do imol = 1,nmol
            write(*,*) molecule_name
            do iatom = 1,natoms_in_mol
            read(10,*) name,x,y,z
            read(12,*) name,(multipole_vec(k),k=1,16)
            write(*,*) name,x,y,z,(multipole_vec(k),k=1,16)
            end do 
         end do
      end do
      end

