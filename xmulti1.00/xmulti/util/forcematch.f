      use parse_text
      implicit none
      common/datafm/nmol_fm,nparam,nstep
      common/lists/model_power,model_list,model_nterms,min_distance,nbondlist
      common/err/error
      common/tfac/tfac
      common/derr/derror
      common/arrays/target_force,target_torque,model_force0,model_torque0,dmodel_force,dmodel_torque,param
      common/fit/errorc,nstepsfit
      integer nmol_fm,nparam,nstep,nstepsfit,nbondlist
      real(8), dimension(3,200,20000) :: target_force,target_torque,model_force0,model_torque0
      real(8), dimension(3,200,50,20000) :: dmodel_force,dmodel_torque
      real(8), dimension(200) :: param
      real(8) :: error,errorc
      real(8), dimension(200) :: derror
      character(len = 2048) :: buffer
      character(len = 2048), dimension(64) :: textlist
      integer imol,ios,n1,n2,i,j
      integer iparam,nitems,n,nb
      integer, dimension(200,100) :: model_power
      integer, dimension(200) :: model_nterms
      integer, dimension(200,2) :: model_list
      
      real(8) :: fx,fy,fz,tx,ty,tz
      real(8) :: dismin,tfac
      real(8), dimension(200) :: min_distance

      open(10,file = 'target_forces.dat')
      open(12,file = 'model_forces.dat')
      open(23,file = 'mindistance.dat')
      open(30,file = 'forcematch.control')

      nstepsfit = 200

      do while(.true.)
         read(30,'(A)',iostat = ios) buffer
         call split_text_to_list(buffer,textlist,nitems)
         if(ios.lt.0) exit
         if(textlist(1).eq.'nstepsfit') then 
            read(textlist(2),*) nstepsfit
         else if(textlist(1).eq.'tfac') then 
               read(textlist(2),*) tfac
         endif
      end do 

      read(12,*) nparam,nmol_fm,nbondlist
      do n = 1,nbondlist
         read(23,*) i,j,dismin
         min_distance(n) = dismin
         read(12,'(A)') buffer
            call split_text_to_list(buffer,textlist,nitems)
            read(textlist(2),*) model_list(n,1)
            read(textlist(3),*) model_list(n,2)
            model_nterms(n) = nitems - 3
            nb = 0 
            do j = 4,nitems
               nb = nb + 1
               read(textlist(j),*) model_power(n,nb)
            end do 
      end do 


      nstep = 0
!     read in the forces
      do while(.true.)
         read(10,*,iostat = ios) n1
         read(12,*,iostat = ios) n2
         if(ios.lt.0) exit
         nstep = nstep + 1
!         write(*,*) 'n1,n2 = ',n1,n2

         do imol = 1,nmol_fm
            read(10,*) i,fx,fy,fz,tx,ty,tz

            target_force(1,imol,nstep) = fx
            target_force(2,imol,nstep) = fy
            target_force(3,imol,nstep) = fz

            target_torque(1,imol,nstep) = tx
            target_torque(2,imol,nstep) = ty
            target_torque(3,imol,nstep) = tz

            read(12,*) i,fx,fy,fz,tx,ty,tz

            model_force0(1,imol,nstep) = fx
            model_force0(2,imol,nstep) = fy
            model_force0(3,imol,nstep) = fz

            model_torque0(1,imol,nstep) = tx
            model_torque0(2,imol,nstep) = ty
            model_torque0(3,imol,nstep) = tz

         end do 

         do imol = 1,nmol_fm
            do iparam = 1,nparam
               read(12,*) i,j,fx,fy,fz,tx,ty,tz
               dmodel_force(1,iparam,imol,nstep) = fx
               dmodel_force(2,iparam,imol,nstep) = fy
               dmodel_force(3,iparam,imol,nstep) = fz

               dmodel_torque(1,iparam,imol,nstep) = tx
               dmodel_torque(2,iparam,imol,nstep) = ty
               dmodel_torque(3,iparam,imol,nstep) = tz
            end do 
         end do 

      end do 

      call read_param
!      call error_test
      call get_error
      write(*,*) 'RMS error = ',dsqrt(error),dsqrt(errorc)
      call minimize

      call print

      end

      subroutine read_param()
      implicit none
      common/arrays/target_force,target_torque,model_force0,model_torque0,dmodel_force,dmodel_torque,param
      common/datafm/nmol_fm,nparam,nstep
      real(8), dimension(3,200,20000) :: target_force,target_torque,model_force0,model_torque0
      real(8), dimension(3,200,50,20000) :: dmodel_force,dmodel_torque
      real(8), dimension(200) :: param
      integer iparam,j,ios
      integer nmol_fm,nparam,nstep

      open(55,file = 'param.dat')
      param = 0

      read(55,*,iostat = ios)
      do iparam = 1,nparam
         read(55,*,iostat = ios) j,param(iparam)
         if(ios.lt.0) exit
      end do 

      end subroutine read_param

      subroutine error_test()
      implicit none
      common/err/error
      common/derr/derror
      common/datafm/nmol_fm,nparam,nstep
      common/arrays/target_force,target_torque,model_force0,model_torque0,dmodel_force,dmodel_torque,param
      common/fit/errorc,nstepsfit
      integer nmol_fm,nparam,nstep,nstepsfit
      real(8), dimension(200) :: derror

      integer iparam
      real(8), dimension(200) :: derror0

      real(8), dimension(3,200,20000) :: target_force,target_torque,model_force0,model_torque0
      real(8), dimension(3,200,50,20000) :: dmodel_force,dmodel_torque
      real(8), dimension(200) :: param
      real(8) :: error,errorc
      real(8) :: error0,tiny

      tiny = 1.d-8

      call get_error()
!      write(*,*) 'error = ',error
      derror0 = derror
      error0 = error

      do iparam = 1,nparam
!         write(*,*) 'jjj',iparam,derror(iparam)
      end do 
!      stop

      do iparam = 1,nparam
         param(iparam) = param(iparam) + tiny 
         call get_error()
         write(*,*) iparam,derror0(iparam),(error-error0)/tiny
         param(iparam) = param(iparam) - tiny
      end do 

      stop
      end subroutine error_test

      subroutine minimize()
      use common_data
      use nr_mod
      implicit none
      real(8) :: xvec(60000)
      real(8) :: ftol,fret
      integer :: iter
      common/arrays/target_force,target_torque,model_force0,model_torque0,dmodel_force,dmodel_torque,param
      common/fit/errorc,nstepsfit
      common/datafm/nmol_fm,nparam,nstep
      common/err/error
      integer nmol_fm,nparam,nstep
      real(8), dimension(3,200,20000) :: target_force,target_torque,model_force0,model_torque0
      real(8), dimension(3,200,50,20000) :: dmodel_force,dmodel_torque
      real(8), dimension(200) :: param
      real(8) :: error,errorc
      real(8), dimension(200) :: derror
      integer n,iparam,nstepsfit


      open(40,file = 'param.out')

      n = 0 
      do iparam = 1,nparam
         n = n + 1
         xvec(n) = param(iparam)
      end do 

      ftol = 0.0d0 
      linmin_param = 1.0d0 
      call frprmn(xvec,n,ftol,iter,fret)

      write(*,*) 'RMS error = ',dsqrt(error),dsqrt(errorc)
      write(40,*) 'error = ',error

      call print_param


      end subroutine minimize

      subroutine print_param
      implicit none
      common/datafm/nmol_fm,nparam,nstep
      common/lists/model_power,model_list,model_nterms,min_distance,nbondlist
      common/arrays/target_force,target_torque,model_force0,model_torque0,dmodel_force,dmodel_torque,param
      integer nmol_fm,nparam,nstep
      real(8), dimension(3,200,20000) :: target_force,target_torque,model_force0,model_torque0
      real(8), dimension(3,200,50,20000) :: dmodel_force,dmodel_torque
      real(8), dimension(200) :: param
      integer, dimension(200) :: model_nterms
      integer, dimension(200,2) :: model_list
      integer, dimension(200,100) :: model_power
      integer nbondlist
      integer iparam,n,nterms,nb
      real(8), dimension(200) :: ppp
      real(8), dimension(200) :: min_distance

      do iparam = 1,nparam
         write(40,*) iparam,param(iparam)
      end do 

      iparam = 0 
      do n = 1,nbondlist
         nterms = model_nterms(n)
         do nb = 1,nterms
            iparam = iparam + 1
            ppp(nb) = param(iparam)
         end do 
         write(40,14) model_list(n,1),model_list(n,2)
 14      format('PAIR  ',2i16)
         write(40,15) (model_power(n,nb),nb=1,nterms)
 15      format('POWER ',10i16)
         write(40,16) (ppp(nb),nb=1,nterms)
 16      format('COEFF ',10f16.5)
         write(40,17) min_distance(n),12.0d0,1.0d0,400.0d0,-300.0d0
 17      format('SPLINE',5f16.5)
      end do 

      end subroutine print_param



      real(8) function func(xvec)
      implicit none
      common/arrays/target_force,target_torque,model_force0,model_torque0,dmodel_force,dmodel_torque,param
      common/datafm/nmol_fm,nparam,nstep
      common/err/error
      integer nmol_fm,nparam,nstep
      real(8), dimension(3,200,20000) :: target_force,target_torque,model_force0,model_torque0
      real(8), dimension(3,200,50,20000) :: dmodel_force,dmodel_torque
      real(8), dimension(200) :: param
      real(8) :: error
      real(8), dimension(200) :: derror

      integer iparam,n
      real(8) :: xvec(60000)

      n = 0 
      do iparam = 1,nparam
         n = n + 1
         param(n) = xvec(n)
      end do 
      call get_error

      func = error

      end function func

      subroutine dfunc(xvec,gvec)
      implicit none
      common/arrays/target_force,target_torque,model_force0,model_torque0,dmodel_force,dmodel_torque,param
      common/datafm/nmol_fm,nparam,nstep
      common/derr/derror
      integer nmol_fm,nparam,nstep
      real(8), dimension(3,200,20000) :: target_force,target_torque,model_force0,model_torque0
      real(8), dimension(3,200,50,20000) :: dmodel_force,dmodel_torque
      real(8), dimension(200) :: param
      real(8) :: error
      real(8), dimension(200) :: derror

      integer iparam,n
      real(8) :: xvec(60000),gvec(60000)

      n = 0 
      do iparam = 1,nparam
         n = n + 1
         gvec(n) = derror(iparam)
      end do 
      
      end subroutine dfunc



     subroutine get_error()
      implicit none
      common/datafm/nmol_fm,nparam,nstep
      common/err/error
      common/derr/derror
      common/arrays/target_force,target_torque,model_force0,model_torque0,dmodel_force,dmodel_torque,param
      common/fit/errorc,nstepsfit
      common/tfac/tfac
      integer nmol_fm,nparam,nstep,nstepsfit
      integer imol,iparam,n

      real(8), dimension(3,200,20000) :: target_force,target_torque,model_force0,model_torque0
      real(8), dimension(3,200,50,20000) :: dmodel_force,dmodel_torque
      real(8), dimension(3,200) :: model_force,model_torque
      real(8), dimension(200) :: param
      real(8) :: error,errorc
      real(8) :: rms_error
      real(8), dimension(200) :: derror
      real(8) :: dfx,dfy,dfz,dtx,dty,dtz,dot
      real(8) :: tfac
      integer nnn,nnnc

      error = 0.0d0 
      errorc = 0.0d0 
      derror = 0.0d0 

      nnn = nstepsfit * nmol_fm * 6
      nnnc = (nstep - nstepsfit) * nmol_fm * 6

      do n = 1,nstep
         call get_model_forces_and_torques(model_force,model_torque,n)

         do imol = 1,nmol_fm
            dfx = target_force(1,imol,n) - model_force(1,imol)
            dfy = target_force(2,imol,n) - model_force(2,imol)
            dfz = target_force(3,imol,n) - model_force(3,imol)

!            write(*,*) target_force(1,imol,n),target_force(2,imol,n),target_force(3,imol,n)
!            write(*,*) model_force(1,imol),model_force(2,imol),model_force(3,imol)
!            write(*,*)

            dtx = target_torque(1,imol,n) - model_torque(1,imol)
            dty = target_torque(2,imol,n) - model_torque(2,imol)
            dtz = target_torque(3,imol,n) - model_torque(3,imol)

            if(n.le.nstepsfit) then 
               error = error + dfx**2 + dfy**2 + dfz**2
               error = error + (dtx**2 + dty**2 + dtz**2) / (tfac**2)
            do iparam = 1,nparam
               dot = dfx * dmodel_force(1,iparam,imol,n) & 
     &             + dfy * dmodel_force(2,iparam,imol,n) &
     &             + dfz * dmodel_force(3,iparam,imol,n) &
     &             + (dtx * dmodel_torque(1,iparam,imol,n) &
     &             + dty * dmodel_torque(2,iparam,imol,n) &
     &             + dtz * dmodel_torque(3,iparam,imol,n) ) / (tfac**2)
               derror(iparam) = derror(iparam) - 2.0d0 * dot / dble(nnn)
            end do

            else if(n.gt.nstepsfit) then
               errorc = errorc + dfx**2 + dfy**2 + dfz**2
               errorc = errorc + (dtx**2 + dty**2 + dtz**2) / (tfac**2)
            endif


         end do 
      end do 


      error = error / dble(nnn)
      errorc = errorc / dble(nnnc)
      rms_error = dsqrt(error)
!      write(*,*) 'rms_error = ',rms_error

      end subroutine get_error


      subroutine get_model_forces_and_torques(model_force,model_torque,n)
      implicit none
      common/datafm/nmol_fm,nparam,nstep
      common/arrays/target_force,target_torque,model_force0,model_torque0,dmodel_force,dmodel_torque,param
      common/fit/errorc,nstepsfit
      real(8), dimension(3,200,20000) :: target_force,target_torque,model_force0,model_torque0
      real(8), dimension(3,200,50,20000) :: dmodel_force,dmodel_torque
      integer nmol_fm,nparam,nstep
      real(8), dimension(3,200) :: model_force,model_torque
      real(8), dimension(200) :: param
      real(8) :: errorc
      integer imol,n,iparam,nstepsfit
         do imol = 1,nmol_fm
            model_force(1,imol) = model_force0(1,imol,n)
            model_force(2,imol) = model_force0(2,imol,n)
            model_force(3,imol) = model_force0(3,imol,n)
            
            model_torque(1,imol) = model_torque0(1,imol,n)
            model_torque(2,imol) = model_torque0(2,imol,n)
            model_torque(3,imol) = model_torque0(3,imol,n)
            
            do iparam = 1,nparam
               model_force(1,imol) = model_force(1,imol) + dmodel_force(1,iparam,imol,n) * param(iparam)
               model_force(2,imol) = model_force(2,imol) + dmodel_force(2,iparam,imol,n) * param(iparam)
               model_force(3,imol) = model_force(3,imol) + dmodel_force(3,iparam,imol,n) * param(iparam)

               model_torque(1,imol) = model_torque(1,imol) + dmodel_torque(1,iparam,imol,n) * param(iparam)
               model_torque(2,imol) = model_torque(2,imol) + dmodel_torque(2,iparam,imol,n) * param(iparam)
               model_torque(3,imol) = model_torque(3,imol) + dmodel_torque(3,iparam,imol,n) * param(iparam)
            end do  

         end do 
         end subroutine get_model_forces_and_torques
      
      subroutine print
      implicit none
      common/datafm/nmol_fm,nparam,nstep
      common/arrays/target_force,target_torque,model_force0,model_torque0,dmodel_force,dmodel_torque,param
      real(8), dimension(3,200,20000) :: target_force,target_torque,model_force0,model_torque0
      real(8), dimension(3,200,50,20000) :: dmodel_force,dmodel_torque
      integer nmol_fm,nparam,nstep
      real(8), dimension(3,200) :: model_force,model_torque
      real(8), dimension(200) :: param
      integer n,imol

      open(63,file = 'fm.out')
      open(140,file = 'fm_force.out')
      open(141,file = 'fm_torque.out')

      do n = 1,nstep
         call get_model_forces_and_torques(model_force,model_torque,n)
         do imol =1,nmol_fm
            write(63,*) n,imol,model_force(1,imol),model_force(2,imol),model_force(3,imol)
            write(63,*) n,imol,target_force(1,imol,n),target_force(2,imol,n),target_force(3,imol,n)
            write(63,*) n,imol,model_torque(1,imol),model_torque(2,imol),model_torque(3,imol)
            write(63,*) n,imol,target_torque(1,imol,n),target_torque(2,imol,n),target_torque(3,imol,n)
            write(63,*) 
            write(140,*) model_force(1,imol),target_force(1,imol,n)
            write(140,*) model_force(2,imol),target_force(2,imol,n)
            write(140,*) model_force(3,imol),target_force(3,imol,n)

            write(141,*) model_torque(1,imol),target_torque(1,imol,n)
            write(141,*) model_torque(2,imol),target_torque(2,imol,n)
            write(141,*) model_torque(3,imol),target_torque(3,imol,n)
         end do 
      end do 
      end 


