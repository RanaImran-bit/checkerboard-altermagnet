
    SUBROUTINE Initph
      use cpmc
      use jiekou,only:ranGen
      integer::iw,i,j,k
      do iw=1,NWLKRS
        do i=1,nsites
        if(ranGen(ISEED)>0.5) then
          ising_ph(i,iw)=1
        else
          ising_ph(i,iw)=-1
        end if
        end do
      end do

    return
    end SUBROUTINE Initph

    SUBROUTINE MkExpVph
    use cpmc
    integer::i,j,k,l,m,n,ising
    real(sp)::tmpc1,tmpc2
    ! sigma,alpha_u,
    !------------------Zhongbing-phonon fields---------------------!
    do i=1,nsites
         do ising=-2,2,2
           tmpc1=-hdeltau*g_ph*ising
           tmpc2=-hdeltau*g_ph*ising
           eph(i,ising,1)=exp( tmpc1 )
           eph(i,ising,2)=exp( tmpc2 )
           Deltaph(i,ising,1)=eph(i,ising,1)-1.0_sp
           Deltaph(i,ising,2)=eph(i,ising,2)-1.0_sp
           if(ising==0) then
           coeffph(i,ising)=0.5*(exp(deltau*w_ph)-exp(-deltau*w_ph))
           else
           coeffph(i,ising)=0.5*(exp(deltau*w_ph)+exp(-deltau*w_ph))
           end if
         end do
    end do
    !---------------------------------------------------------------!
    end SUBROUTINE MkExpVph

      SUBROUTINE Fullph(iw)
      use cpmc
      use jiekou, only:ranGen
      integer::i,ising
      ! ,j,k,l,m,n
      integer::iw,itmph,itmp

      integer::ibas
      real(sp)::ovlpNEW,ovlpOLD
      real(sp)::tmpup,tmpdn
      real(sp)::rdet(NWFBAS,-1:1,NSPIN),pt(-1:1),ptsum,ptest
      real(sp)::gl_up(NUP,nsites,NWFBAS),gr_up(nsites,NUP,NWFBAS)
      real(sp)::gl_dn(NDN,nsites,NWFBAS),gr_dn(nsites,NDN,NWFBAS)
      real(sp)::gup_ii(NWFBAS),gdn_ii(NWFBAS),r_detbas(NWFBAS)
      real(sp)::ggx1_up(1,NUP),ggx2_up(NUP,1)
      real(sp)::ggx1_dn(1,NDN),ggx2_dn(NDN,1)
      real(sp)::ggx1_up1(1,NUP),ggx2_up2(NUP,1)
      real(sp)::ggx1_dn1(1,NDN),ggx2_dn2(NDN,1)
      external :: dgemm

      !!*****call MkExpV             ! Make e^V for each spin

      ovlpOLD=ovlpDET(iw)

      !Get Ising Variables, according to important sampling

   do 100 i=1,nsites

      !For fixed i; up & down; N^2
      !Make gl_s(k,i) = [ PhiT^+(ibas) Phi ]^-1(k,j) * PhiT^+(j,i); Ns x nsites
      !!!!!!gr_s(i,k) = Phi(i,j) * [ PhiT^+(ibas) Phi ]^-1(j,k);    nsites x Ns

      do 120 ibas=1,NWFBAS

          ggx2_up(:,1)=phiT_up(i,:,ibas)

          call dgemm('N', 'N', NUP, 1, NUP, 1.0_8, g_up(:,:,ibas), NUP, ggx2_up, NUP, 0.0_8, ggx2_up2, NUP) !cao
          ggx2_up=ggx2_up2

          !ggx2_up=matmul(g_up(:,:,ibas),ggx2_up)
          gl_up(:,i,ibas)=ggx2_up(:,1)

          ggx1_up(1,:)=phi_up(i,:,iw)

          call dgemm('N', 'N', 1, NUP, NUP, 1.0_8, ggx1_up, 1, g_up(:,:,ibas), NUP, 0.0_8, ggx1_up1, 1) !cao
          ggx1_up=ggx1_up1

          !ggx1_up=matmul(ggx1_up,g_up(:,:,ibas))
          gr_up(i,:,ibas)=ggx1_up(1,:)

          ggx2_dn(:,1)=phiT_dn(i,:,ibas)

          call dgemm('N', 'N', NDN, 1, NDN, 1.0_8, g_dn(:,:,ibas), NDN, ggx2_dn, NDN, 0.0_8, ggx2_dn2, NDN) !cao
          ggx2_dn=ggx2_dn2

          !ggx2_dn=matmul(g_dn(:,:,ibas),ggx2_dn)
          gl_dn(:,i,ibas)=ggx2_dn(:,1)

          ggx1_dn(1,:)=phi_dn(i,:,iw)

          call dgemm('N', 'N', 1, NDN, NDN, 1.0_8, ggx1_dn, 1, g_dn(:,:,ibas), NDN, 0.0_8, ggx1_dn1, 1) !cao
          ggx1_dn=ggx1_dn1

          !ggx1_dn=matmul(ggx1_dn,g_dn(:,:,ibas))
          gr_dn(i,:,ibas)=ggx1_dn(1,:)

      120 continue

      do 140 ibas=1,NWFBAS

          gup_ii(ibas)=dot_product(gr_up(i,:,ibas),phiT_up(i,:,ibas))
          gdn_ii(ibas)=dot_product(gr_dn(i,:,ibas),phiT_dn(i,:,ibas))

      140 continue

         !! check x(i), p( x(i) )

         do ising=-1,1,2

          do ibas=1,NWFBAS
           rdet(ibas,ising,1)=1.0+Deltaph(i,ising+ising_ph(i,iw),1)*gup_ii(ibas)
           rdet(ibas,ising,2)=1.0+Deltaph(i,ising+ising_ph(i,iw),2)*gdn_ii(ibas)

           r_detbas(ibas)=rdet(ibas,ising,1)*rdet(ibas,ising,2)*detbas(ibas)
           !! cwfbas(ibas) is included(ovlps)

          end do
          ovlpNEW=sum(r_detbas)

          pt(ising)=coeffph(i,ising+ising_ph(i,iw))*ovlpNEW/ovlpOLD
          if(ovlpNEW*sgn(iw) <= ZERO) then
           wgtwlkr(iw)=wgtwlkr(iw)/(1.0-pt(ising))
           pt(ising)=zero
          end if

         end do

         !! calculate p(1)/( p(1)+p(-1) )

         ptsum=pt(1)+pt(-1)

         if(ptsum == zero) then
          wgtwlkr(iw)=ZERO     ! walker is terminated, nothing left
          return
         end if

         ptest=pt(1)/ptsum
         if( ptest > ranGen(ISEED) ) then
          ising_sum(i)=1+ising_ph(i,iw)
          ising_ph(i,iw)=1
          ovlpNEW=pt(1)*ovlpOLD/coeffph(i,ising_sum(i)) !see the above line
         else
          ising_sum(i)=-1+ising_ph(i,iw)
          ising_ph(i,iw)=-1
          ovlpNEW=pt(-1)*ovlpOLD/coeffph(i,ising_sum(i))!see the above line
         end if

         itmph=ising_ph(i,iw)
         itmp=ising_sum(i)                     ! picked Ising_sum(i)

         !! New Wave function, only i_th row are changed
         !! advance phi by exp(-deltau*V) for picked Ising variables
         tmpup=eph(i,itmp,1)
         tmpdn=eph(i,itmp,2)

         phi_up(i,:,iw)=tmpup*phi_up(i,:,iw)
         phi_dn(i,:,iw)=tmpdn*phi_dn(i,:,iw)

         !! Updating Green's(?) functions due to change of ising_u(i)
      do 200 ibas=1,NWFBAS

         tmpup=Deltaph(i,itmp,1)/rdet(ibas,itmph,1)

         ggx2_up(:,1)=gl_up(:,i,ibas)
         ggx1_up(1,:)=gr_up(i,:,ibas)
         g_up(:,:,ibas)=g_up(:,:,ibas)-tmpup*matmul(ggx2_up,ggx1_up)

         tmpdn=Deltaph(i,itmp,2)/rdet(ibas,itmph,2)

         ggx2_dn(:,1)=gl_dn(:,i,ibas)
         ggx1_dn(1,:)=gr_dn(i,:,ibas)
         g_dn(:,:,ibas)=g_dn(:,:,ibas)-tmpdn*matmul(ggx2_dn,ggx1_dn)

      200 continue

         !! compute determinant of the overlap integral
         do ibas=1,NWFBAS
          detbas(ibas)=rdet(ibas,itmph,1)*rdet(ibas,itmph,2)*detbas(ibas)
         end do

         ovlpOLD=ovlpNEW                ! ready for next move
         wgtwlkr(iw)=ptsum*wgtwlkr(iw)  ! renormalized walker weight

  100   continue

        ovlpDET(iw)=ovlpOLD

        return
        END subroutine Fullph
      !--------------------------------------------------!
