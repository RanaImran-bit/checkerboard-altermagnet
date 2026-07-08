subroutine InitV
  use cpmc
  integer:: i,j,k,i1,j1,ixy,jxy
  integer::nls(2)
  ! interaction channel definition
  ! channel 1: same orbital Hubbard U
  spinlsi(1)=1;spinlsj(1)=2 ! up and down spin
  do i=1,nsites
      Vlist(i,1) = uxx
  enddo

  ! channel 2: interaction between orbital U'
  ! got 4 channels each for up and down spin coupling
  do i=1,nsites
      Vlist(i,2:5) = uxy
  enddo
  k=1
  do i=1,2
    do j=1,2
      k=k+1
      spinlsi(k)=i
      spinlsj(k)=j
      orblsi(k)=1
      orblsj(k)=2
    end do
  end do
  ! channel 3: interaction between nearest neighbor
  ! between x and y orbital, its negative
  ! between same orbital its positive
  do j1=1,2 ! neighbor index
    do ixy=1,2 ! xy orbital of i site
      do jxy=1,2 ! xy orbital of j site
        do i=1,2
          do j=1,2
            k=k+1
            spinlsi(k)=i
            spinlsj(k)=j
            orblsi(k)=ixy
            orblsj(k)=jxy
            do i1=1,nsites
              if (orblsi(k)==orblsj(k)) then
                Vlist(i1,k) = v          ! per-channel: +v intra-orbital
              else
                Vlist(i1,k) = -v         ! -v inter-orbital (was Vlist(i1,6:37), which
              endif                      ! overwrote ALL v-channels with the last sign)
            enddo
          end do
        end do
      end do
    end do
  end do
  nls(1)=NUP
  nls(2)=NDN
  do i=1,channels
    nlsi(i)=nls(spinlsi(i))
    nlsj(i)=nls(spinlsj(i))
  enddo
endsubroutine InitV

      ! subroutine VHubb(iw)
      !   use cpmc
      !   integer:: iw,i
      !   do i=1, nsites
      !     call Vee(iw, i, i, 1, 2, 1)
      !   enddo
      ! endsubroutine VHubb

      ! subroutine Vxy(iw)
      !   use cpmc
      !   integer:: iw,i
      !   ! between x and y orbital
      !   do i=1, lxy
      !     call Vee(iw, i, i+lxy, 1, 1, 2)
      !     call Vee(iw, i, i+lxy, 1, 2, 3)
      !     call Vee(iw, i, i+lxy, 2, 1, 4)
      !     call Vee(iw, i, i+lxy, 2, 2, 5)
      !   enddo
      ! endsubroutine Vxy

      ! subroutine Vznn(iw)
      !   use cpmc
      !   integer:: iw,i,j,j1,ix,iy,chan,ixy,jxy
      !   do i=1, lxy
      !     do j1=1,2
      !       j=idis(i,j1) ! neighbor of i
      !       do ixy=0,1 ! xy orbital of i site
      !         do jxy=0,1 ! xy orbital of j site
      !           isite = i + ixy*lxy
      !           jsite = j + jxy*lxy

      !           chan =  6 + ( (((j1-1)*2 + ixy)*2 + jxy) * 4 + 0 )
      !           call Vee(iw, isite, jsite, 1, 1, chan)

      !           chan = 6 + ( (((j1-1)*2 + ixy)*2 + jxy) * 4 + 1 )
      !           call Vee(iw, isite, jsite, 1, 2, chan)

      !           chan = 6 + ( (((j1-1)*2 + ixy)*2 + jxy) * 4 + 2 )
      !           call Vee(iw, isite, jsite, 2, 1, chan)

      !           chan = 6 + ( (((j1-1)*2 + ixy)*2 + jxy) * 4 + 3 )
      !           call Vee(iw, isite, jsite, 2, 2, chan)

      !         end do
      !       end do
      !     enddo
      !   enddo
      ! endsubroutine Vznn

subroutine VHubb(iw)
  use cpmc
  integer:: iw,i
  do i=1, nsites
    call Vee(iw, i, i, 1)
  enddo
endsubroutine VHubb

subroutine Vxy(iw)
  use cpmc
  integer:: iw,i
  ! between x and y orbital
  do i=1, lxy
    call Vee(iw, i, i+lxy, 2)
    call Vee(iw, i, i+lxy, 3)
    call Vee(iw, i, i+lxy, 4)
    call Vee(iw, i, i+lxy, 5)
  enddo
endsubroutine Vxy

subroutine Vznn(iw)
  use cpmc
  integer:: iw,i,j,j1,ix,iy,chan,ixy,jxy
  do i=1, lxy
    do j1=1,2
      j=idis(i,j1) ! neighbor of i
      do ixy=0,1 ! xy orbital of i site
        do jxy=0,1 ! xy orbital of j site
          isite = i + ixy*lxy
          jsite = j + jxy*lxy

          chan =  6 + ( (((j1-1)*2 + ixy)*2 + jxy) * 4 + 0 )
          call Vee(iw, isite, jsite, chan)

          chan = 6 + ( (((j1-1)*2 + ixy)*2 + jxy) * 4 + 1 )
          call Vee(iw, isite, jsite, chan)

          chan = 6 + ( (((j1-1)*2 + ixy)*2 + jxy) * 4 + 2 )
          call Vee(iw, isite, jsite, chan)

          chan = 6 + ( (((j1-1)*2 + ixy)*2 + jxy) * 4 + 3 )
          call Vee(iw, isite, jsite, chan)
        end do
      end do
    enddo
  enddo
endsubroutine Vznn


      SUBROUTINE MkExpV
      use cpmc
      use jiekou,only:invcosh
      use, intrinsic :: ieee_arithmetic
      logical :: has_nan
      integer::i,j,k,l,m,n,chan,ising
      real(sp)::sigma,alpha_u,tmpc1,tmpc2
      real(sp) :: invcosh_arg
      do chan=1, 37
        do i=1,nsites
          invcosh_arg = deltau * abs(Vlist(i,chan))
          call invcosh(invcosh_arg, alpha_u)

            if(Vlist(i,chan)>-0.001) then
             do ising=-1,1,2
               tmpc1= alpha_u*ising-hdeltau*Vlist(i,chan)
               tmpc2=-alpha_u*ising-hdeltau*Vlist(i,chan)

               expV(i,ising,1,chan)=exp( tmpc1 )
               expV(i,ising,2,chan)=exp( tmpc2 )

               DeltaV(i,ising,1,chan)=expV(i,ising,1,chan)-1.0_sp
               DeltaV(i,ising,2,chan)=expV(i,ising,2,chan)-1.0_sp
               coeffv(i,ising,chan)=1.0
             end do
            else

            ! call invcosh(-1.0*deltau*Vlist(i,chan),alpha_u)
             do ising=-1,1,2
               tmpc1=alpha_u*ising-hdeltau*Vlist(i,chan)
               tmpc2=alpha_u*ising-hdeltau*Vlist(i,chan)
               
               expV(i,ising,1,chan)=exp( tmpc1 )
               expV(i,ising,2,chan)=exp( tmpc2 )

               DeltaV(i,ising,1,chan)=expV(i,ising,1,chan)-1.0_sp
               DeltaV(i,ising,2,chan)=expV(i,ising,2,chan)-1.0_sp
               coeffv(i,ising,chan)=exp(-1.0*alpha_u*ising+hdeltau*Vlist(i,chan))
             end do
            end if
        end do
      enddo

      return
      end subroutine MkExpV

      !~~~~~Cosh(a) = Exp[ x/2 ]
      SUBROUTINE invcosh(x,a)
      use cpmc
      implicit none
      real(sp)::x,a,tmpxx
      integer::i,j,k,l,m,n
      tmpxx=exp(x/2.0)
      a=log( tmpxx+sqrt( tmpxx*tmpxx-1.0 ) )
      return
      end subroutine invcosh

      ! unify the propagation of electron electron interaction

      subroutine Vee(iw, isite, jsite, chan)
        use cpmc
        use jiekou, only: ranGen   ! NB: not bare `use jiekou` — jiekou declares
                                   ! Vee's own interface; importing it here makes
                                   ! `Vee` an ambiguous self-reference under gfortran.
        use, intrinsic :: ieee_arithmetic
        implicit none
        logical :: has_nan_in_detpbas 
      
        ! Input
        integer, intent(in) :: iw, isite, jsite, chan
      
        ! Local vars
        integer :: i, j, ibas, ising
        integer :: ni, nj
        real(sp) :: ovlpNEW, ovlpOLD, ptsum, ptest
        real(sp) :: rdet(NWFBAS, -1:1, 2), pt(-1:1), r_detbas(NWFBAS)
        real(sp) :: tmp1, tmp2
        logical :: same_spin
      
        ! Spin-dependent pointers into global arrays
        real(sp), allocatable :: phi_i(:,:), phi_j(:,:)
        real(sp), allocatable :: phiT_i(:,:,:), phiT_j(:,:,:)
        real(sp), allocatable :: g_i(:,:,:), g_j(:,:,:)
        ! Work arrays allocated after knowing ni/nj
        real(sp), allocatable :: gl_i(:,:,:), gr_i(:,:,:)
        real(sp), allocatable :: gl_j(:,:,:), gr_j(:,:,:)
        real(sp), allocatable :: gi_ii(:), gj_ii(:)
        real(sp), allocatable :: ggx1_i(:,:), ggx2_i(:,:), ggx1_i1(:,:), ggx2_i2(:,:)
        real(sp), allocatable :: ggx1_j(:,:), ggx2_j(:,:), ggx1_j1(:,:), ggx2_j2(:,:)
        ! Temporary G-matrix for sequential update
        real(sp), allocatable :: g_temp_i(:,:,:)
      
        i = isite
        j = jsite
        ovlpOLD = ovlpDET(iw)

        ! =========================== ALLOCATE ARRAYS ===========================
        if (use_spinor) then
          ni = NE
          nj = NE
          allocate(phi_i(nsites,NE), phi_j(nsites,NE))
          allocate(phiT_i(nsites,NE,NWFBAS), phiT_j(nsites,NE,NWFBAS))
          allocate(g_i(NE,NE,NWFBAS), g_j(NE,NE,NWFBAS))
        else
          ni = nlsi(chan)
          nj = nlsj(chan)
          allocate(phi_i(nsites,ni), phi_j(nsites,nj))
          allocate(phiT_i(nsites,ni,NWFBAS), phiT_j(nsites,nj,NWFBAS))
          allocate(g_i(ni,ni,NWFBAS), g_j(nj,nj,NWFBAS))
        end if
      
        allocate(gl_i(ni,nsites,NWFBAS), gr_i(nsites,ni,NWFBAS))
        allocate(gl_j(nj,nsites,NWFBAS), gr_j(nsites,nj,NWFBAS))
        allocate(gi_ii(NWFBAS), gj_ii(NWFBAS))
        allocate(ggx1_i(1,ni), ggx2_i(ni,1), ggx1_i1(1,ni), ggx2_i2(ni,1))
        allocate(ggx1_j(1,nj), ggx2_j(nj,1), ggx1_j1(1,nj), ggx2_j2(nj,1))
        allocate(g_temp_i(ni, ni, NWFBAS)) ! ni == nj if same_spin

        ! =========================== ASSIGN POINTERS ===========================
        if (use_spinor) then
          g_i    = g(:,:,:)
          g_j    = g(:,:,:)
          phi_i  = phi(:,:,iw)
          phi_j  = phi(:,:,iw)
          phiT_i = phiT(:,:,:)
          phiT_j = phiT(:,:,:)
        else
          if (spinlsi(chan) == 1) then
            g_i    = g_up(:,:,:)
            phi_i  = phi_up(:,:,iw)
            phiT_i = phiT_up(:,:,:)
          else
            g_i    = g_dn(:,:,:)
            phi_i  = phi_dn(:,:,iw)
            phiT_i = phiT_dn(:,:,:)
          end if
      
          if (spinlsj(chan) == 1) then
            g_j    = g_up(:,:,:)
            phi_j  = phi_up(:,:,iw)
            phiT_j = phiT_up(:,:,:)
          else
            g_j    = g_dn(:,:,:)
            phi_j  = phi_dn(:,:,iw)
            phiT_j = phiT_dn(:,:,:)
          end if
        end if
      
        ! ==================== BUILD LEFT/RIGHT OBJECTS (SITE i) ====================
        do ibas = 1, NWFBAS
          ! Site i
          ggx2_i(:,1) = phiT_i(i,:,ibas)
          call dgemm('N','N', ni, 1, ni, 1.0_8, g_i(:,:,ibas), ni, ggx2_i, ni, 0.0_8, ggx2_i2, ni)
          gl_i(:,i,ibas) = ggx2_i2(:,1)
      
          ggx1_i(1,:) = phi_i(i,:)
          call dgemm('N','N', 1, ni, ni, 1.0_8, ggx1_i, 1, g_i(:,:,ibas), ni, 0.0_8, ggx1_i1, 1)
          gr_i(i,:,ibas) = ggx1_i1(1,:)
        end do
      
        ! ==================== DIAGONAL TERMS (SITE i) ====================
        do ibas = 1, NWFBAS
          gi_ii(ibas) = dot_product(gr_i(i,:,ibas), phiT_i(i,:,ibas))
        end do

        ! === CHECK IF UPDATES ARE COUPLED (SAME SPIN) OR DECOUPLED (MIXED SPIN) ===
        same_spin = use_spinor .or. (spinlsi(chan) == spinlsj(chan))

        if (.not. same_spin) then
            ! Mixed-spin case: pre-calculate j-site objects from original G_j
            ! These are independent of the i-site update.
            do ibas = 1, NWFBAS
                ggx2_j(:,1) = phiT_j(j,:,ibas)
                call dgemm('N','N', nj, 1, nj, 1.0_8, g_j(:,:,ibas), nj, ggx2_j, nj, 0.0_8, ggx2_j2, nj)
                gl_j(:,j,ibas) = ggx2_j2(:,1)
            
                ggx1_j(1,:) = phi_j(j,:)
                call dgemm('N','N', 1, nj, nj, 1.0_8, ggx1_j, 1, g_j(:,:,ibas), nj, 0.0_8, ggx1_j1, 1)
                gr_j(j,:,ibas) = ggx1_j1(1,:)
            
                gj_ii(ibas) = dot_product(gr_j(j,:,ibas), phiT_j(j,:,ibas))
            end do
        end if
      
        ! ==================== IMPORTANCE SAMPLING OVER ising = ±1 ====================
        do ising = -1, 1, 2
          do ibas = 1, NWFBAS
            ! 1. Rdet for site i (always calculated first)
            rdet(ibas,ising,1) = 1.0_sp + DeltaV(i,ising,1,chan) * gi_ii(ibas)

            if (same_spin) then
                ! Same-spin case: j-site update depends on i-site update
                ! 1a. Calculate temporary G matrix after i-update
                tmp1 = DeltaV(i,ising,1,chan) / rdet(ibas,ising,1)
                ggx2_i(:,1) = gl_i(:,i,ibas)
                ggx1_i(1,:) = gr_i(i,:,ibas)
                g_temp_i(:,:,ibas) = g_i(:,:,ibas) - tmp1 * matmul(ggx2_i, ggx1_i)

                ! 1b. Calculate j-site objects using temporary G
                ggx2_j(:,1) = phiT_j(j,:,ibas)
                call dgemm('N','N', nj, 1, nj, 1.0_8, g_temp_i(:,:,ibas), nj, ggx2_j, nj, 0.0_8, ggx2_j2, nj)
                ! gl_j(:,j,ibas) = ggx2_j2(:,1) ! Not needed for rdet
                
                ggx1_j(1,:) = phi_j(j,:)
                call dgemm('N','N', 1, nj, nj, 1.0_8, ggx1_j, 1, g_temp_i(:,:,ibas), nj, 0.0_8, ggx1_j1, 1)
                ! gr_j(j,:,ibas) = ggx1_j1(1,:) ! Not needed for rdet
                
                ! 1c. Calculate j-site diagonal term from temp objects
                gj_ii(ibas) = dot_product(ggx1_j1(1,:), phiT_j(j,:,ibas))
                
                ! 2. Rdet for site j
                rdet(ibas,ising,2) = 1.0_sp + DeltaV(i,ising,2,chan) * gj_ii(ibas)
            else
                ! Mixed-spin case: independent, use pre-calculated gj_ii
                rdet(ibas,ising,2) = 1.0_sp + DeltaV(i,ising,2,chan) * gj_ii(ibas)
            end if

            ! 3. Total Rdet
            r_detbas(ibas) = rdet(ibas,ising,1) * rdet(ibas,ising,2) * detbas(ibas)
          end do
          ovlpNEW = sum(r_detbas)
          pt(ising) = coeffv(i,ising,chan) * ovlpNEW / ovlpOLD / two
          if (ovlpNEW * sgn(iw) <= zero) then
            wgtwlkr(iw) = wgtwlkr(iw) / (1.0_sp - pt(ising))
            pt(ising) = zero
          end if
        end do
      
        ! ==================== COLLAPSE THE FIELD ====================
        ptsum = pt(1) + pt(-1)
        if (ptsum == zero) then
          write(*,'(A,2I6,1X,2(ES12.4,1X),2(ES12.4,1X))') 'PTSUM=0 i,chan=', i, chan, pt(1), pt(-1), ovlpOLD, sgn(iw)
          wgtwlkr(iw) = zero
          return
        end if

        ptest = pt(1) / ptsum
        if (ptest > ranGen(ISEED)) then
          ising_v(i,chan) = 1
          ovlpNEW = pt(1) * ovlpOLD * two / coeffv(i,1,chan)
        else
          ising_v(i,chan) = -1
          ovlpNEW = pt(-1) * ovlpOLD * two / coeffv(i,-1,chan)
        end if
      
        ! ==================== UPDATE WAVEFUNCTION ====================
        if (use_spinor) then
          phi(i,:,iw) = expV(i, ising_v(i,chan), 1, chan) * phi(i,:,iw)
          phi(j,:,iw) = expV(i, ising_v(i,chan), 2, chan) * phi(j,:,iw) ! <-- FIXED: expV uses 'i'
        else
          if (spinlsi(chan) == 1) then
            phi_up(i,:,iw) = expV(i, ising_v(i,chan), 1, chan) * phi_up(i,:,iw)
          else
            phi_dn(i,:,iw) = expV(i, ising_v(i,chan), 1, chan) * phi_dn(i,:,iw)
          end if
      
          if (spinlsj(chan) == 1) then
            phi_up(j,:,iw) = expV(i, ising_v(i,chan), 2, chan) * phi_up(j,:,iw) ! <-- FIXED: expV uses 'i'
          else
            phi_dn(j,:,iw) = expV(i, ising_v(i,chan), 2, chan) * phi_dn(j,:,iw) ! <-- FIXED: expV uses 'i'
          end if
        end if
      
        ! ==================== UPDATE GREEN'S FUNCTION ====================
        do ibas = 1, NWFBAS
          tmp1 = DeltaV(i, ising_v(i,chan), 1, chan) / rdet(ibas, ising_v(i,chan), 1)
          tmp2 = DeltaV(i, ising_v(i,chan), 2, chan) / rdet(ibas, ising_v(i,chan), 2)
      
          ! Site i update (always applies to g_i)
          ggx2_i(:,1) = gl_i(:,i,ibas)
          ggx1_i(1,:) = gr_i(i,:,ibas)
          g_i(:,:,ibas) = g_i(:,:,ibas) - tmp1 * matmul(ggx2_i, ggx1_i)
      
          ! Site j update
          if (same_spin) then
            ! Must re-calculate gl_j/gr_j using the *updated* g_i
            ggx2_j(:,1) = phiT_j(j,:,ibas)
            call dgemm('N','N', nj, 1, nj, 1.0_8, g_i(:,:,ibas), nj, ggx2_j, nj, 0.0_8, ggx2_j2, nj) ! Use updated g_i
            
            ggx1_j(1,:) = phi_j(j,:)
            call dgemm('N','N', 1, nj, nj, 1.0_8, ggx1_j, 1, g_i(:,:,ibas), nj, 0.0_8, ggx1_j1, 1) ! Use updated g_i

            ! Apply update to g_i (which is also g_j)
            g_i(:,:,ibas) = g_i(:,:,ibas) - tmp2 * matmul(ggx2_j2, ggx1_j1)
          else
            ! Mixed-spin: update g_j using pre-calculated gl_j/gr_j
            ggx2_j(:,1) = gl_j(:,j,ibas)
            ggx1_j(1,:) = gr_j(j,:,ibas)
            g_j(:,:,ibas) = g_j(:,:,ibas) - tmp2 * matmul(ggx2_j, ggx1_j)
          end if
        end do
      
        ! ==================== FINALIZE ====================
        ovlpDET(iw) = ovlpNEW
        wgtwlkr(iw) = ptsum * wgtwlkr(iw)
      
        ! ==================== DEALLOCATE ====================
        deallocate(gl_i, gr_i, gl_j, gr_j, gi_ii, gj_ii)
        deallocate(ggx1_i, ggx2_i, ggx1_i1, ggx2_i2)
        deallocate(ggx1_j, ggx2_j, ggx1_j1, ggx2_j2)
        deallocate(phi_i, phi_j, phiT_i, phiT_j, g_i, g_j)
        deallocate(g_temp_i)
      
        return
      end subroutine Vee
      
      ! subroutine Vee(iw, isite, jsite, chan)
      !   use cpmc
      !   use jiekou
      !   implicit none
      !   integer, intent(in) :: iw, isite, jsite, chan
      !   integer :: i, j, ibas, ising
      !   integer :: ni, nj
      !   real(sp) :: ovlpNEW, ovlpOLD, ptsum, ptest
      !   real(sp) :: rdet(NWFBAS,-1:1,NSPIN), pt(-1:1), r_detbas(NWFBAS)
      !   real(sp) :: tmp1, tmp2

      !   ! Spin-dependent pointers into global arrays
      !   real(sp), allocatable :: phi_i(:,:), phi_j(:,:)
      !   real(sp), allocatable :: phiT_i(:,:,:), phiT_j(:,:,:)
      !   real(sp), allocatable :: g_i(:,:,:), g_j(:,:,:)

      !   ! Work arrays allocated after knowing ni/nj
      !   real(sp), allocatable :: gl_i(:,:,:), gr_i(:,:,:)
      !   real(sp), allocatable :: gl_j(:,:,:), gr_j(:,:,:)
      !   real(sp), allocatable :: gi_ii(:), gj_ii(:)
      !   real(sp), allocatable :: ggx1_i(:,:), ggx2_i(:,:), ggx1_i1(:,:), ggx2_i2(:,:)
      !   real(sp), allocatable :: ggx1_j(:,:), ggx2_j(:,:), ggx1_j1(:,:), ggx2_j2(:,:)

      !   if (use_spinor) then
      !     allocate(phi_i(nsites,NE), phi_j(nsites,NE))
      !     allocate(phiT_i(nsites,NE,NWFBAS), phiT_j(nsites,NE,NWFBAS))
      !     allocate(g_i(NE,NE,NWFBAS), g_j(NE,NE,NWFBAS))
      !     allocate(gl_i(NE,nsites,NWFBAS), gr_i(nsites,NE,NWFBAS))
      !     allocate(gl_j(NE,nsites,NWFBAS), gr_j(nsites,NE,NWFBAS))
      !     allocate(gi_ii(NWFBAS), gj_ii(NWFBAS))
      !     allocate(ggx1_i(1,NE), ggx2_i(NE,1), ggx1_i1(1,NE), ggx2_i2(NE,1))
      !     allocate(ggx1_j(1,NE), ggx2_j(NE,1), ggx1_j1(1,NE), ggx2_j2(NE,1))
      !     g_i    = g(:,:,:)
      !     phi_i  = phi(:,:,iw)
      !     phiT_i = phiT(:,:,:)
      !     g_j    = g(:,:,:)
      !     phi_j  = phi(:,:,iw)
      !     phiT_j = phiT(:,:,:)
      !   else
      !     allocate(phi_i(nsites,nlsi(chan)), phi_j(nsites,nlsj(chan)))
      !     allocate(phiT_i(nsites,nlsi(chan),NWFBAS), phiT_j(nsites,nlsj(chan),NWFBAS))
      !     allocate(g_i(nlsi(chan),nlsi(chan),NWFBAS), g_j(nlsj(chan),nlsj(chan),NWFBAS))
      !     allocate(gl_i(nlsi(chan),nsites,NWFBAS), gr_i(nsites,nlsi(chan),NWFBAS))
      !     allocate(gl_j(nlsj(chan),nsites,NWFBAS), gr_j(nsites,nlsj(chan),NWFBAS))
      !     allocate(gi_ii(NWFBAS), gj_ii(NWFBAS))
      !     allocate(ggx1_i(1,nlsi(chan)), ggx2_i(nlsi(chan),1), ggx1_i1(1,nlsi(chan)), ggx2_i2(nlsi(chan),1))
      !     allocate(ggx1_j(1,nlsj(chan)), ggx2_j(nlsj(chan),1), ggx1_j1(1,nlsj(chan)), ggx2_j2(nlsj(chan),1))
      !     if (spinlsi(chan)==1) then
      !       g_i    = g_up(:,:,:)
      !       phi_i  = phi_up(:,:,iw)
      !       phiT_i = phiT_up(:,:,:)
      !     else
      !       g_i    = g_dn(:,:,:)
      !       phi_i  = phi_dn(:,:,iw)
      !       phiT_i = phiT_dn(:,:,:)
      !     endif
      !     if (spinlsj(chan)==1) then
      !       g_j    = g_up(:,:,:)
      !       phi_j  = phi_up(:,:,iw)
      !       phiT_j = phiT_up(:,:,:)
      !     else
      !       g_j    = g_dn(:,:,:)
      !       phi_j  = phi_dn(:,:,iw)
      !       phiT_j = phiT_dn(:,:,:)
      !     endif
      !   endif

      !   if (use_spinor) then
      !     ni = NE
      !     nj = NE

      !     ! ===== Build left/right objects for site i and site j (spinor) =====
      !     do ibas = 1, NWFBAS
      !       ! ---- site i ----
      !       ggx2_i(:,1) = phiT_i(i,:,ibas)
      !       call dgemm('N','N', ni, 1, ni, 1.0_8, g_i(:,:,ibas), ni, ggx2_i, ni, 0.0_8, ggx2_i2, ni)
      !       ggx2_i = ggx2_i2
      !       gl_i(:,i,ibas) = ggx2_i(:,1)

      !       ggx1_i(1,:) = phi_i(i,:)
      !       call dgemm('N','N', 1, ni, ni, 1.0_8, ggx1_i, 1, g_i(:,:,ibas), ni, 0.0_8, ggx1_i1, 1)
      !       ggx1_i = ggx1_i1
      !       gr_i(i,:,ibas) = ggx1_i(1,:)

      !       ! ---- site j ----
      !       ggx2_j(:,1) = phiT_j(j,:,ibas)
      !       call dgemm('N','N', nj, 1, nj, 1.0_8, g_j(:,:,ibas), nj, ggx2_j, nj, 0.0_8, ggx2_j2, nj)
      !       ggx2_j = ggx2_j2
      !       gl_j(:,j,ibas) = ggx2_j(:,1)

      !       ggx1_j(1,:) = phi_j(j,:)
      !       call dgemm('N','N', 1, nj, nj, 1.0_8, ggx1_j, 1, g_j(:,:,ibas), nj, 0.0_8, ggx1_j1, 1)
      !       ggx1_j = ggx1_j1
      !       gr_j(j,:,ibas) = ggx1_j(1,:)
      !     end do

      !     ! Diagonal contractions
      !     do ibas = 1, NWFBAS
      !       gi_ii(ibas) = dot_product( gr_i(i,:,ibas), phiT_i(i,:,ibas) )
      !       gj_ii(ibas) = dot_product( gr_j(j,:,ibas), phiT_j(j,:,ibas) )
      !     end do

      !     ! ===== Importance sampling over ising = ±1 =====
      !     do ising = -1, 1, 2
      !       do ibas = 1, NWFBAS
      !         rdet(ibas,ising,1) = 1.0_sp + DeltaV(i,ising,1,chan) * gi_ii(ibas)
      !         rdet(ibas,ising,2) = 1.0_sp + DeltaV(i,ising,2,chan) * gj_ii(ibas)
      !         r_detbas(ibas) = rdet(ibas,ising,1) * rdet(ibas,ising,2) * detbas(ibas)
      !       end do
      !       ovlpNEW = sum(r_detbas)
      !       pt(ising) = coeffv(i,ising,chan) * ovlpNEW / ovlpOLD / two
      !       if (ovlpNEW*sgn(iw) <= zero) then
      !          wgtwlkr(iw) = wgtwlkr(iw) / ( 1.0_sp - pt(ising) )
      !          pt(ising) = zero
      !       end if
      !     end do

      !     ptsum = pt(1) + pt(-1)
      !     if (ptsum == zero) then
      !     write(*,'(A,2I6,1X,2(ES12.4,1X),2(ES12.4,1X))') &
      !     'PTSUM=0  i,chan=', i, chan, pt(1), pt(-1), ovlpOLD, sgn(iw)
      !        wgtwlkr(iw) = zero
      !        return
      !     end if

      !     ptest = pt(1)/ptsum
      !     if (ptest > ranGen(ISEED)) then
      !        ising_v(i,chan) =  1
      !        ovlpNEW = pt(1) * ovlpOLD * two / coeffv(i,1,chan)
      !     else
      !        ising_v(i,chan) = -1
      !        ovlpNEW = pt(-1) * ovlpOLD * two / coeffv(i,-1,chan)
      !     end if
      !     if (ising_v(i,chan) == 0) write(*,'(A,3I8)') 'SET ising_v(abs_i,chan,val)=', i, chan, ising_v(i,chan)

      !     ! ===== Update Slater rows (spinor) and Green's functions for the chosen field =====
      !     phi(i,:,iw) = expV(i, ising_v(i,chan), 1, chan) * phi(i,:,iw)
      !     phi(j,:,iw) = expV(j, ising_v(i,chan), 2, chan) * phi(j,:,iw)

      !     do ibas = 1, NWFBAS
      !       ! updates are per-basis
      !       tmp1 = DeltaV(i, ising_v(i,chan), 1,chan) / rdet(ibas, ising_v(i,chan), 1)
      !       tmp2 = DeltaV(i, ising_v(i,chan), 2,chan) / rdet(ibas, ising_v(i,chan), 2)

      !       ! Update G for site i
      !       ggx2_i(:,1) = gl_i(:,i,ibas)
      !       ggx1_i(1,:) = gr_i(i,:,ibas)
      !       g(:,:,ibas)  = g(:,:,ibas) - tmp1 * matmul(ggx2_i, ggx1_i)

      !       ! Update G for site j
      !       ggx2_j(:,1) = gl_j(:,j,ibas)
      !       ggx1_j(1,:) = gr_j(j,:,ibas)
      !       g(:,:,ibas)  = g(:,:,ibas) - tmp2 * matmul(ggx2_j, ggx1_j)
      !     end do

      !     ovlpDET(iw) = ovlpNEW
      !     wgtwlkr(iw) = ptsum * wgtwlkr(iw)
      !   else
      !   ovlpOLD = ovlpDET(iw)
      !   ! ===== Build left/right objects for site i (spin ispin) and site j (spin jspin) =====
      !   do ibas = 1, NWFBAS
      !     ! ---- spin at i ----
      !     ggx2_i(:,1) = phiT_i(i,:,ibas)
      !     call dgemm('N','N', ni, 1, ni, 1.0_8, g_i(:,:,ibas), ni, ggx2_i, ni, 0.0_8, ggx2_i2, ni)
      !     ggx2_i = ggx2_i2
      !     gl_i(:,i,ibas) = ggx2_i(:,1)

      !     ggx1_i(1,:) = phi_i(i,:)
      !     call dgemm('N','N', 1, ni, ni, 1.0_8, ggx1_i, 1, g_i(:,:,ibas), ni, 0.0_8, ggx1_i1, 1)
      !     ggx1_i = ggx1_i1
      !     gr_i(i,:,ibas) = ggx1_i(1,:)

      !     ! ---- spin at j ----
      !     ggx2_j(:,1) = phiT_j(j,:,ibas)
      !     call dgemm('N','N', nj, 1, nj, 1.0_8, g_j(:,:,ibas), nj, ggx2_j, nj, 0.0_8, ggx2_j2, nj)
      !     ggx2_j = ggx2_j2
      !     gl_j(:,j,ibas) = ggx2_j(:,1)

      !     ggx1_j(1,:) = phi_j(j,:)
      !     call dgemm('N','N', 1, nj, nj, 1.0_8, ggx1_j, 1, g_j(:,:,ibas), nj, 0.0_8, ggx1_j1, 1)
      !     ggx1_j = ggx1_j1
      !     gr_j(j,:,ibas) = ggx1_j(1,:)
      !   end do

      !   ! Diagonal contractions
      !   do ibas = 1, NWFBAS
      !     gi_ii(ibas) = dot_product( gr_i(i,:,ibas), phiT_i(i,:,ibas) )
      !     gj_ii(ibas) = dot_product( gr_j(j,:,ibas), phiT_j(j,:,ibas) )
      !   end do

      !   ! ===== Importance sampling over ising = ±1 =====
      !   do ising = -1, 1, 2
      !     do ibas = 1, NWFBAS
      !       rdet(ibas,ising,1) = 1.0_sp + DeltaV(i,ising,1,chan) * gi_ii(ibas)
      !       rdet(ibas,ising,2) = 1.0_sp + DeltaV(i,ising,2,chan) * gj_ii(ibas)
      !       r_detbas(ibas) = rdet(ibas,ising,1) * rdet(ibas,ising,2) * detbas(ibas)
      !     end do
      !     ovlpNEW = sum(r_detbas)
      !     pt(ising) = coeffv(i,ising,chan) * ovlpNEW / ovlpOLD / two
      !     if (ovlpNEW*sgn(iw) <= zero) then
      !        wgtwlkr(iw) = wgtwlkr(iw) / ( 1.0_sp - pt(ising) )
      !        pt(ising) = zero
      !     end if
      !   end do


      !   ptsum = pt(1) + pt(-1)
      !   if (ptsum == zero) then
      !   write(*,'(A,2I6,1X,2(ES12.4,1X),2(ES12.4,1X))') &
      !   'PTSUM=0  i,chan=', i, chan, pt(1), pt(-1), ovlpOLD, sgn(iw)
      !      wgtwlkr(iw) = zero
      !      return
      !   end if

      !   ptest = pt(1)/ptsum
      !   if (ptest > ranGen(ISEED)) then
      !      ising_v(i,chan) =  1
      !      ovlpNEW = pt(1) * ovlpOLD * two / coeffv(i,1,chan)
      !   else
      !      ising_v(i,chan) = -1
      !      ovlpNEW = pt(-1) * ovlpOLD * two / coeffv(i,-1,chan)
      !   end if

      !   ! ===== Update Slater rows and Green's functions for the chosen field =====
      !   if (ispin == 1) then
      !      phi_up(i,:,iw) = expV(i, ising_v(i,chan), 1, chan) * phi_up(i,:,iw)
      !   else
      !      phi_dn(i,:,iw) = expV(i, ising_v(i,chan), 1, chan) * phi_dn(i,:,iw)
      !   end if

      !   if (jspin == 1) then
      !      phi_up(j,:,iw) = expV(j, ising_v(i,chan), 2, chan) * phi_up(j,:,iw)
      !   else
      !      phi_dn(j,:,iw) = expV(j, ising_v(i,chan), 2, chan) * phi_dn(j,:,iw)
      !   end if


      !   do ibas = 1, NWFBAS
      !     ! updates are per-basis

      !     tmp1 = DeltaV(i, ising_v(i,chan), 1,chan) / rdet(ibas, ising_v(i,chan), 1)
      !     tmp2 = DeltaV(i, ising_v(i,chan), 2,chan) / rdet(ibas, ising_v(i,chan), 2)

      !     ! Update G for spin at i
      !     ggx2_i(:,1) = gl_i(:,i,ibas)
      !     ggx1_i(1,:) = gr_i(i,:,ibas)
      !     if (ispin == 1) then
      !       g_up(:,:,ibas) = g_up(:,:,ibas) - tmp1 * matmul(ggx2_i, ggx1_i)
      !     else
      !       g_dn(:,:,ibas) = g_dn(:,:,ibas) - tmp1 * matmul(ggx2_i, ggx1_i)
      !     end if

      !     ! Update G for spin at j
      !     ggx2_j(:,1) = gl_j(:,j,ibas)
      !     ggx1_j(1,:) = gr_j(j,:,ibas)
      !     if (jspin == 1) then
      !       g_up(:,:,ibas) = g_up(:,:,ibas) - tmp2 * matmul(ggx2_j, ggx1_j)
      !     else
      !       g_dn(:,:,ibas) = g_dn(:,:,ibas) - tmp2 * matmul(ggx2_j, ggx1_j)
      !     end if
      !   end do

      !   ovlpDET(iw) = ovlpNEW
      !   wgtwlkr(iw) = ptsum * wgtwlkr(iw)
      ! endif ! no spinor

      !   ! cleanup
      !   deallocate(gl_i, gr_i, gl_j, gr_j, gi_ii, gj_ii, ggx1_i, ggx2_i, ggx1_i1, ggx2_i2, ggx1_j, ggx2_j, ggx1_j1, ggx2_j2)

      !   return
      ! end subroutine Vee