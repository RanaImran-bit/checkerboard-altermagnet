!-----------------HoKinAnHop------------------!
!Purpose: For convenient input in other plotting
      SUBROUTINE cpOut
      use cpmc
      integer:: i,j,k
        open(unit=50,file='out.dat',status='unknown')
        write(50,*) 'n_site'
        write(50,*) nsites
        write(50,*) 'n_e'
        write(50,*) NELEC
        write(50,*) 'filling'
        write(50,*) NELEC/nsites
        write(50,*) 'n_up'
        write(50,*) NUP
        write(50,*) 'n_dn'
        write(50,*) NDN
        write(50,*) 'Ud'
        write(50,*) ud
        write(50,*) 't0'
        write(50,*) t0
        write(50,*) 't1'
        write(50,*) t1
        write(50,*) 't2'
        write(50,*) t2
        write(50,*) 'Vpd'
        write(50,*) vpd
        !write(50,*) 'alpha'
        !write(50,*) alpha
        write(50,*) 'tam'
        write(50,*) tam        
        write(50,*) 'alphat1'
        write(50,*) alphat1
        write(50,*) 'n_walker'
        write(50,*) NWLKRS
        write(50,*) 'dt'
        write(50,*) deltau
        write(50,*) 'bpStep'
        write(50,*) itvl_m
        write(50,*) 'ke'
        write(50,*) cor_top(NAVE-2,1),cor_top(NAVE-2,6)
        write(50,*) 'pe'
        write(50,*) cor_top(NAVE-1,1),cor_top(NAVE-1,6)
        write(50,*) 'totalEn'
        write(50,*) cor_top(NAVE,1),cor_top(NAVE,6)
        write(50,*) 'kSpace'
        do i=1,NSTATES
            write(50,'(2f15.9)') kSet(i,1),kSet(i,2)
        enddo
        write(50,*) 'Npair'
        k=(3*NSTATES+9)*NSTATES+5*NPAIR
        do i=1,NSTATES
            k=k+1
            write(50,'(2f15.9)')   cor_top(k,1),cor_top(k,6)
        enddo
        write(50,*) 'Dk'
        k=(3*NSTATES+9)*NSTATES+5*NPAIR
        do i=1,NSTATES
            k=k+1
            write(50,'(2f15.9)')   cor_top(k,1),cor_top(k,6)
        enddo
        close(50)
      endsubroutine
!-----------------HoKinAnHop------------------!
