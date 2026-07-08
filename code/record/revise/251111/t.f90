Subroutine Initlatt
  use cpmc
  !======================================================================!
  ! Two-orbital square lattice (d_xz, d_yz) with real-space hoppings for:
  !   eps_x (kx,ky)  = -2 t1 cos kx - 2 t2 cos ky - 4 t3 cos kx cos ky
  !   eps_y (kx,ky)  = -2 t2 cos kx - 2 t1 cos ky - 4 t3 cos kx cos ky
  !   eps_xy(kx,ky)  = -4 t4 sin kx sin ky
  !
  ! Mapping to real space:
  !   - NN intra-orbital: along x,y as above (sign = -t)
  !   - NNN (diagonals) intra-orbital: four diagonals with amplitude -t3
  !   - Inter-orbital (x <-> y) on diagonals with signs:
  !       +t4 on (±1,±1) with same sign dx=dy,
  !       -t4 on (±1,∓1) with opposite signs.
  !======================================================================!

  integer :: i, ix, iy, orb

  !! Build coordinates and orbital index
  i = 0
  do ix = 0, lx-1
    do iy = 0, ly-1
      ! orbital 1 ≡ d_xz
      i = i + 1
      ixv(i)     = ix
      iyv(i)     = iy
      sublatt(i) = 1
      iposit(ix,iy,1) = i
      iposit(ix,iy,2) = i
    end do
  end do
  do ix = 0, lx-1
      do iy = 0, ly-1
      ! orbital 2 ≡ d_yz
      i = i + 1
      ixv(i)     = ix
      iyv(i)     = iy
      sublatt(i) = 2
      iposit(ix,iy,2) = i
    end do
  end do

  !! Periodic extension of iposit so that (ix,iy) can be used with +/- shifts
  do orb = 1, 2
    do ix = -2*lx, 2*lx
      do iy = -2*ly, 2*ly
        iposit(ix,iy,orb) = iposit(mod(ix+lx,lx), mod(iy+ly,ly), orb)
      end do
    end do
  end do

  !! One-body hopping matrix (spin-independent)
  tk_ud = 0.0d0

  do i = 1, nsites
    ix  = ixv(i)
    iy  = iyv(i)
    orb = sublatt(i)

    ! -------- Intra-orbital nearest neighbours --------
    if (orb == 1) then
      ! d_xz: -t1 along x, -t2 along y
      tk_ud( i, iposit(ix+1,iy,1), 1 ) = tk_ud( i, iposit(ix+1,iy,1), 1 ) - t1
      tk_ud( i, iposit(ix-1,iy,1), 1 ) = tk_ud( i, iposit(ix-1,iy,1), 1 ) - t1
      tk_ud( i, iposit(ix,iy+1,1), 1 ) = tk_ud( i, iposit(ix,iy+1,1), 1 ) - t2
      tk_ud( i, iposit(ix,iy-1,1), 1 ) = tk_ud( i, iposit(ix,iy-1,1), 1 ) - t2
    else
      ! d_yz: -t2 along x, -t1 along y
      tk_ud( i, iposit(ix+1,iy,2), 1 ) = tk_ud( i, iposit(ix+1,iy,2), 1 ) - t2
      tk_ud( i, iposit(ix-1,iy,2), 1 ) = tk_ud( i, iposit(ix-1,iy,2), 1 ) - t2
      tk_ud( i, iposit(ix,iy+1,2), 1 ) = tk_ud( i, iposit(ix,iy+1,2), 1 ) - t1
      tk_ud( i, iposit(ix,iy-1,2), 1 ) = tk_ud( i, iposit(ix,iy-1,2), 1 ) - t1
    end if

    ! -------- Intra-orbital next-nearest neighbours (diagonals) --------
    tk_ud( i, iposit(ix+1,iy+1,orb), 1 ) = tk_ud( i, iposit(ix+1,iy+1,orb), 1 ) - t3
    tk_ud( i, iposit(ix-1,iy-1,orb), 1 ) = tk_ud( i, iposit(ix-1,iy-1,orb), 1 ) - t3
    tk_ud( i, iposit(ix+1,iy-1,orb), 1 ) = tk_ud( i, iposit(ix+1,iy-1,orb), 1 ) - t3
    tk_ud( i, iposit(ix-1,iy+1,orb), 1 ) = tk_ud( i, iposit(ix-1,iy+1,orb), 1 ) - t3

    ! -------- Inter-orbital hybridization on diagonals (d_xz <-> d_yz) --------
    if (orb == 1) then
      tk_ud( i, iposit(ix+1,iy+1,2), 1 ) = tk_ud( i, iposit(ix+1,iy+1,2), 1 ) + t4
      tk_ud( i, iposit(ix-1,iy-1,2), 1 ) = tk_ud( i, iposit(ix-1,iy-1,2), 1 ) + t4
      tk_ud( i, iposit(ix+1,iy-1,2), 1 ) = tk_ud( i, iposit(ix+1,iy-1,2), 1 ) - t4
      tk_ud( i, iposit(ix-1,iy+1,2), 1 ) = tk_ud( i, iposit(ix-1,iy+1,2), 1 ) - t4
    else
      tk_ud( i, iposit(ix+1,iy+1,1), 1 ) = tk_ud( i, iposit(ix+1,iy+1,1), 1 ) + t4
      tk_ud( i, iposit(ix-1,iy-1,1), 1 ) = tk_ud( i, iposit(ix-1,iy-1,1), 1 ) + t4
      tk_ud( i, iposit(ix+1,iy-1,1), 1 ) = tk_ud( i, iposit(ix+1,iy-1,1), 1 ) - t4
      tk_ud( i, iposit(ix-1,iy+1,1), 1 ) = tk_ud( i, iposit(ix-1,iy+1,1), 1 ) - t4
    end if
  end do

  tk_ud(:,:,2) = tk_ud(:,:,1)

  if (use_spinor) then
    tk = 0.0d0
    tk(1:nsites,1:nsites) = tk_ud(1:nsites,1:nsites,1)
    tk(nsites+1:NSO,nsites+1:NSO) = tk_ud(1:nsites,1:nsites,2)
    if (orb == 1) then
      ! d_xz: -t1 along x, -t2 along y
      tk( i, nsites+iposit(ix+1,iy,1) ) = tk( i, nsites+iposit(ix+1,iy,1) ) - SOC
      tk( i, nsites+iposit(ix-1,iy,1) ) = tk( i, nsites+iposit(ix-1,iy,1) ) - SOC
      tk( i, nsites+iposit(ix,iy+1,1) ) = tk( i, nsites+iposit(ix,iy+1,1) ) - SOC
      tk( i, nsites+iposit(ix,iy-1,1) ) = tk( i, nsites+iposit(ix,iy-1,1) ) - SOC
    else
      ! d_yz: -t2 along x, -t1 along y
      tk( i, nsites+iposit(ix+1,iy,2) ) = tk( i, nsites+iposit(ix+1,iy,2) ) - SOC
      tk( i, nsites+iposit(ix-1,iy,2) ) = tk( i, nsites+iposit(ix-1,iy,2) ) - SOC
      tk( i, nsites+iposit(ix,iy+1,2) ) = tk( i, nsites+iposit(ix,iy+1,2) ) - SOC
      tk( i, nsites+iposit(ix,iy-1,2) ) = tk( i, nsites+iposit(ix,iy-1,2) ) - SOC
    end if
  ! else
    ! tk(:,:,1) = 0.5d0 * ( tk(:,:,1) + transpose( tk(:,:,1) ) )
    ! tk(:,:,2) = 0.5d0 * ( tk(:,:,2) + transpose( tk(:,:,2) ) )
  endif
  
  ! Optional: neighbor list used elsewhere
  do i = 1, nsites
    idis(i,1) = iposit(ixv(i)+1, iyv(i), sublatt(i))
    idis(i,2) = iposit(ixv(i),   iyv(i)+1, sublatt(i))
  end do

end subroutine Initlatt

SUBROUTINE InitPhi
  use cpmc
  implicit none
  integer :: ibas, ip, is

  if (use_spinor) then    
    do ibas = 1, NWFBAS
      phiT(:, 1:NE, ibas) = z(:, 1:NE)
      phiB(:, 1:NE, ibas) = z(:, 1:NE)
      phiZ(:,:,ibas)      = transpose(phiT(:,:,ibas))
   end do
      
  else 
      if (iRead==1) then
        open(22,file='wfup.txt',status='old')
        read(22,*) phiT_up(:,:,1)
        close(22)
        
        open(22,file='wfdn.txt',status='old')
        read(22,*) phiT_dn(:,:,1)
        close(22)
    
      else
        do ip=1,NUP
          do is=1,nsites
             phiT_up(is,ip,1)=z_ud(is,ip,1)
          end do
        end do

       do ip=1,NDN
          do is=1,nsites
             phiT_dn(is,ip,1)=z_ud(is,ip,2)
          end do
       end do

      endif
      
      do ibas=1,NWFBAS
        phiB_up(:,:,ibas) = phiT_up(:,:,1)
        phiB_dn(:,:,ibas) = phiT_dn(:,:,1)
        phiZ_up(:,:,ibas) = transpose(phiT_up(:,:,1))
        phiZ_dn(:,:,ibas) = transpose(phiT_dn(:,:,1))
      end do

  end if
end subroutine InitPhi   

SUBROUTINE MkExpT(xdeltau)
  use cpmc
  use jiekou, only: tred2, tql2
  implicit none
  integer :: i, j, k, s, ierr, nsp
  real(sp) :: xdeltau, sum
  external :: dgemm
  
  if (use_spinor) then
    ! --- (use_spinor = .true.) ---
    
    if (any(isnan(T))) then
        if (myid == 0) print *, "FATAL ERROR in MkExpT: Input matrix T contains NaN!"
    end if

    call tred2(NSO, NSO, T, d, e, z)
    call tql2 (NSO, NSO, d, e, z, ierr)
    
    if (ierr /= 0) then
        write(*,*) 'FATAL ERROR in MkExpT: tql2 returned ierr =', ierr
    end if
    
    call InitPhi 
    
    ! expT = Z * diag(exp(-Δτ d)) * Z^T
    do i = 1, NSO
      do j = 1, NSO
        sum = 0.0_sp
        do k = 1, NSO
            sum = sum + z(i,k) * exp(-xdeltau * d(k)) * z(j,k)
        end do
        expT(i,j) = sum
      end do
    end do
    ! exp2T = expT * expT
    call dgemm('N','N', NSO, NSO, NSO, 1.0_8, expT, NSO, expT, NSO, 0.0_8, exp2T, NSO)

  else
    ! --- (use_spinor = .false.) ---
    
    do s = 1, 2

      if (any(isnan(tk_ud(:,:,s)))) then
          if (myid == 0) then
             print *, "FATAL ERROR in MkExpT: Input matrix tk_ud(:,:,", s, ") contains NaN!"
          end if
      end if
      
      call tred2(nsites, nsites, tk_ud(:,:,s), d_ud(:,s), e_ud(:,s), z_ud(:,:,s))
      call tql2 (nsites, nsites, d_ud(:,s), e_ud(:,s), z_ud(:,:,s), ierr)
      
      if (ierr /= 0) then
         write(*,*) 'FATAL ERROR in MkExpT: tql2 failed for spin', s, ' ierr=', ierr
      end if
      
      ! call InitPhi 
      
      do i = 1, nsites
        do j = 1, nsites
            sum = 0.0_sp
            do k = 1, nsites
              sum = sum + z_ud(i,k,s) * exp(-xdeltau * d_ud(k,s)) * z_ud(j,k,s)
            end do
            expT_ud(i,j,s) = sum
        end do
      end do
      call dgemm('N','N', nsites, nsites, nsites, 1.0_8, expT_ud(:,:,s), nsites, expT_ud(:,:,s), nsites, 0.0_8, exp2T_ud(:,:,s), nsites)
    end do
    
    call InitPhi
    
  end if

END subroutine MkExpT

! SUBROUTINE MkExpT(xdeltau)
!   use cpmc
!   use jiekou, only: tred2, tql2
!   implicit none
!   integer :: i, j, k, s, ierr, nsp
!   real(sp) :: xdeltau, sum
!   external :: dgemm
!   if (use_spinor) then
!     T = tk(:,:)
!     call tred2(NSO, NSO, T, d, e, z)
!     call tql2 (NSO, NSO, d, e, z, ierr)
!     if (ierr /= 0) then
!        write(*,*) 'MkExpT: tql2 returned ierr =', ierr
!     end if
!     call InitPhi
!     ! expT = Z * diag(exp(-Δτ d)) * Z^T
!     do i = 1, NSO
!        do j = 1, NSO
!           sum = 0.0_sp
!           do k = 1, NSO
!              sum = sum + z(i,k) * exp(-xdeltau * d(k)) * z(j,k)
!           end do
!           expT(i,j) = sum
!        end do
!     end do
!     ! exp2T = expT * expT
!     call dgemm('N','N', NSO, NSO, NSO, 1.0_8, expT, NSO, expT, NSO, 0.0_8, exp2T, NSO)

!   else
!     do s = 1, 2
!       call tred2(nsites, nsites, tk_ud(:,:,s), d_ud(:,s), e_ud(:,s), z_ud(:,:,s))
!       call tql2 (nsites, nsites, d_ud(:,s), e_ud(:,s), z_ud(:,:,s), ierr)
!       if (ierr /= 0) write(*,*) 'MkExpT tql2 ierr spin', s, ierr
!       call InitPhi
!       do i = 1, nsites
!         do j = 1, nsites
!             sum = 0.0_sp
!             do k = 1, nsites
!               sum = sum + z_ud(i,k,s) * exp(-xdeltau * d_ud(k,s)) * z_ud(j,k,s)
!             end do
!             expT_ud(i,j,s) = sum
!         end do
!       end do
!       call dgemm('N','N', nsites, nsites, nsites, 1.0_8, expT_ud(:,:,s), nsites, expT_ud(:,:,s), nsites, 0.0_8, exp2T_ud(:,:,s), nsites)
!     end do
!   end if

! END subroutine MkExpT

  !**************************************************!

  SUBROUTINE HalfT(iw,icase,istp)
  use cpmc
  use jiekou,only:Ovlps

  integer::iw,icase
  integer::i,j,k,l,m,n
  external :: dgemm

  ! real(sp)::tmp_up(nsites,NUP),tmp_dn(nsites,NDN)
  if (use_spinor) then
    !! Advance phi by exp(-deltau*T/2) for each spin
      if(icase==1) then
        call dgemm('N', 'N', NSO, NE, NSO, 1.0_8, expT(:,:), NSO, phi(:,:,iw), NSO, 0.0_8, tmp_ud, NSO) !cao
        phi(:,:,iw)=tmp_ud
      else
        call dgemm('N', 'N', NSO, NE, NSO, 1.0_8, exp2T(:,:), NSO, phi(:,:,iw), NSO, 0.0_8, tmp_ud, NSO) !cao
        phi(:,:,iw)=tmp_ud
      end if
    !! Compute determinant of the overlap integral, and Inverse
    if (maxval(abs(phi(:,:,iw))) <= 1.0d-10) then
      write(*,*) 'phi_up(:,:,iw) min/max = 0 =>', istp, minval(phi(:,:,iw)), maxval(phi(:,:,iw))
    end if
  else
  !! Advance phi by exp(-deltau*T/2) for each spin
    if(icase==1) then
      call dgemm('N', 'N', nsites, NUP, nsites, 1.0_8, expT_ud(:,:,1), nsites, phi_up(:,:,iw), nsites, 0.0_8, tmp_up, nsites) !cao
      phi_up(:,:,iw)=tmp_up
      call dgemm('N', 'N', nsites, NDN, nsites, 1.0_8, expT_ud(:,:,2), nsites, phi_dn(:,:,iw), nsites, 0.0_8, tmp_dn, nsites) !cao
      phi_dn(:,:,iw)=tmp_dn
      else
      call dgemm('N', 'N', nsites, NUP, nsites, 1.0_8, exp2T_ud(:,:,1), nsites, phi_up(:,:,iw), nsites, 0.0_8, tmp_up, nsites) !cao
      phi_up(:,:,iw)=tmp_up
      call dgemm('N', 'N', nsites, NDN, nsites, 1.0_8, exp2T_ud(:,:,2), nsites, phi_dn(:,:,iw), nsites, 0.0_8, tmp_dn, nsites) !cao
      phi_dn(:,:,iw)=tmp_dn
    end if
  !! Compute determinant of the overlap integral, and Inverse
  if (maxval(abs(phi_up(:,:,iw))) <= 1.0d-10) then
    write(*,*) 'phi_up(:,:,iw) min/max = 0 =>', istp, minval(phi_up(:,:,iw)), maxval(phi_up(:,:,iw))
  end if
  endif
  call Ovlps(iw)

  return
  END subroutine HalfT

  !**************************************************!
