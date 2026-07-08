!==================SECOND PART=======================!

function ranGen(idum)
use cpmc
implicit none

real(sp)::ranGen
integer,intent(inout)::idum

integer,parameter::ia=16807,im=2147483647,iq=127773,ir=2836
real(sp),save::am
integer,save::ix=-1,iy=-1,kp

if(idum<=0 .or. iy<0) then
am=nearest(1.0,-1.0)/im
iy=ior(ieor(888889999,abs(idum)),1)
ix=ieor(777755555,abs(idum))
idum=abs(idum)+1
end if

ix=ieor(ix,ishft(ix,13))
ix=ieor(ix,ishft(ix,-17))
ix=ieor(ix,ishft(ix,5))
kp=iy/iq
iy=ia*(iy-kp*iq)-ir*kp
if(iy<0) iy=iy+im
ranGen=am*ior(iand(im,ieor(ix,iy)),1)
return
end function ranGen


subroutine stblz
use cpmc
use jiekou,only:modgs
implicit none

integer::iw
real(sp)::rescale
!-----------------HoKinPara-------------------!
do iw=iwStart(myID),iwEnd(myID)
!-----------------HoKinPara-------------------!
   !!skip walkers with zero weight.

   if(wgtwlkr(iw)/=0.0_sp) then

     call modgs(phi_up(:,:,iw),NSTATES,NUP,rescale)
     ovlpDET(iw)=rescale*ovlpDET(iw)

     call modgs(phi_dn(:,:,iw),NSTATES,NDN,rescale)
     ovlpDET(iw)=rescale*ovlpDET(iw)

   end if

end do

return
end subroutine stblz

subroutine modgs(phi,nss,npp,rescale)
use cpmc
use jiekou,only:snrm2
implicit none

real(sp)::rescale,tmp,hld,sdot
integer::nss,npp,ip,jp,n1

real(sp),dimension(:,:)::phi

rescale=1.0_sp
n1=1

do ip=1,npp

   tmp=1.0_sp/snrm2(nss,phi(:,ip),n1)
   rescale=rescale*tmp
   phi(:,ip)=tmp*phi(:,ip)

   do jp=ip+1,npp
   hld=dot_product(phi(:,ip),phi(:,jp))
   phi(:,jp)=phi(:,jp)-hld*phi(:,ip)
   end do
end do

return
end subroutine modgs

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

      SUBROUTINE SGEFA(A,N,N1,IPVT,INFO)
      use cpmc
      implicit none

      integer::N,N1,INFO
      integer,dimension(:)::IPVT
      REAL(sp),dimension(:,:)::A
      REAL(sp)::T,smax,xmag
      integer::isamax,J,K,KP1,L,NM1,I
      INFO = 0
      NM1 = N - 1
      IF (NM1 < 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1

         L=maxval(maxloc(abs(A(K:N,K))))+K-1

         IPVT(K) = L
         IF (A(L,K) == 0.0E0_sp) GO TO 40
            IF (L == K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
            T = -1.0E0_sp/A(K,K)

            A(K+1:N,K)=T*A(K+1:N,K)

            DO 30 J = KP1, N
               T = A(L,J)
               IF (L == K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE


               A(K+1:N,J)=A(K+1:N,J)+T*A(K+1:N,K)


   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) == 0.0E0_sp) INFO = N
      RETURN
      END SUBROUTINE SGEFA


      SUBROUTINE SGEDI(A,N,N1,IPVT,DET,WORK,JOB)
      use cpmc
      implicit none

      integer::N,N1,JOB
      integer,dimension(:)::IPVT
      REAL(sp),dimension(:)::DET,WORK
      real(sp),dimension(N)::tmp1
      real(sp),dimension(:,:)::A
      REAL(sp)::T
      integer::I,J,K,KB,KP1,L,NM1
      IF (JOB/10 == 0) GO TO 70
         DET(1) = 1.0E0_sp
         DET(2) = 0.0E0_sp
         DO 50 I = 1, N
            IF (IPVT(I) /= I) DET(1) = -DET(1)
	            DET(1) = A(I,I)*DET(1)
            IF (DET(1) == 0.0E0_sp) EXIT

            do
               IF (ABS(DET(1)) >= 1.0E0_sp) EXIT
               DET(1) = TEN*DET(1)
               DET(2) = DET(2) - 1.0E0_sp
            end do

            do
               IF (ABS(DET(1)) < TEN) EXIT
               DET(1) = DET(1)/TEN
               DET(2) = DET(2) + 1.0E0_sp
            end do

   50    CONTINUE

   70  CONTINUE
      IF (MOD(JOB,10) == 0) GO TO 150
         DO 100 K = 1, N
            A(K,K) = 1.0E0_sp/A(K,K)
            T = -A(K,K)

            A(1:K-1,K)=T*A(1:K-1,K)

            KP1 = K + 1
            IF (N < KP1) CYCLE
            do J = KP1, N
               T = A(K,J)
               A(K,J) = 0.0E0_sp


               A(1:K,J)=A(1:K,J)+T*A(1:K,K)


            end do

  100    CONTINUE
         NM1 = N - 1
         IF (NM1 < 1) GO TO 140
         DO 130 KB = 1, NM1
            K = N - KB
            KP1 = K + 1

            do I = KP1, N
               WORK(I) = A(I,K)
               A(I,K) = 0.0E0_sp
            end do

            do J = KP1, N

               T = WORK(J)
               A(:,K)=A(:,K)+T*A(:,J)

            end do

            L = IPVT(K)

            IF (L /= K) then
            tmp1=A(:,K);A(:,K)=A(:,L);A(:,L)=tmp1
            end if

  130    CONTINUE
  140    CONTINUE
  150 CONTINUE
      RETURN
      END SUBROUTINE SGEDI

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

      FUNCTION SNRM2(N,SX,INCX)
      use cpmc
      implicit none

      real(sp)::SNRM2
      real(sp),dimension(:)::SX
      integer::N,NN,I,J,INCX,NEXT
      REAL(sp)::CUTLO, CUTHI, HITEST, SUM, XMAX, ch

      CUTLO=4.441E-50_sp
      CUTHI=1.304E+50_sp
      ch=1.0_sp*N

      IF(N > 0) GO TO 10
         SNRM2  = ZERO
         GO TO 300
   10 NEXT=30
      SUM = ZERO
      NN = N * INCX
      I = 1
   20 select case(NEXT)
        case (30); goto 30
        case (50); goto 50
        case (70); goto 70
        case (110); goto 110
      end select

   30 IF( ABS(SX(I)) > CUTLO) GO TO 85
      NEXT=50
      XMAX = ZERO
   50 IF( SX(I) == ZERO) GO TO 200
      IF( ABS(SX(I)) > CUTLO) GO TO 85
      NEXT=70
      GO TO 105
  100 I = J
      NEXT=110
      SUM = (SUM / SX(I)) / SX(I)
  105 XMAX = ABS(SX(I))
      GO TO 115
   70 IF( ABS(SX(I)) > CUTLO ) GO TO 75
  110 IF( ABS(SX(I)) <= XMAX ) GO TO 115
         SUM = ONE + SUM * (XMAX / SX(I))**2
         XMAX = ABS(SX(I))
         GO TO 200
  115 SUM = SUM + (SX(I)/XMAX)**2
      GO TO 200
   75 SUM = (SUM * XMAX) * XMAX

   85 HITEST = CUTHI/ch
      DO 95 J =I,NN,INCX
      IF(ABS(SX(J)) >= HITEST) GO TO 100
   95    SUM = SUM + SX(J)**2
      SNRM2 = SQRT( SUM )
      GO TO 300
  200 CONTINUE
      I = I + INCX
      IF ( I <= NN ) GO TO 20
      SNRM2 = XMAX * SQRT(SUM)
  300 CONTINUE
      RETURN
      END FUNCTION SNRM2

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

      subroutine tred2(nm,n,a,d,e,z)
      use cpmc
      implicit none

      integer::i,j,k,l,n,ii,nm,jp1
      real(sp),dimension(:,:)::a,z
      real(sp),dimension(:)::d,e
      real(sp)::f,g,h,hh,scale

       z = a

      if (n == 1) go to 320
      do 300 ii = 2, n
         i = n + 2 - ii
         l = i - 1
         h = 0.0e0_sp
         scale = 0.0e0_sp
         if (l < 2) go to 130
         do 120 k = 1, l
  120    scale = scale + abs(z(i,k))
         if (scale /= 0.0e0_sp) go to 140
  130    e(i) = z(i,l)
         go to 290
  140    do 150 k = 1, l
            z(i,k) = z(i,k) / scale
            h = h + z(i,k) * z(i,k)
  150    continue
         f = z(i,l)
         g = -sign(sqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         z(i,l) = f - g
         f = 0.0e0_sp
         do 240 j = 1, l
            z(j,i) = z(i,j) / h
            g = 0.0e0_sp
            do 180 k = 1, j
  180       g = g + z(j,k) * z(i,k)
            jp1 = j + 1
            if (l < jp1) go to 220
            do 200 k = jp1, l
  200       g = g + z(k,j) * z(i,k)
  220       e(j) = g / h
            f = f + e(j) * z(i,j)
  240    continue
         hh = f / (h + h)
         do 260 j = 1, l
            f = z(i,j)
            g = e(j) - hh * f
            e(j) = g
            do 260 k = 1, j
               z(j,k) = z(j,k) - f * e(k) - g * z(i,k)
  260    continue
  290    d(i) = h
  300 continue
  320 d(1) = 0.0e0_sp
      e(1) = 0.0e0_sp
      do 500 i = 1, n
         l = i - 1
         if (d(i) == 0.0e0_sp) go to 380
         do 360 j = 1, l
            g = 0.0e0_sp
            do 340 k = 1, l
  340       g = g + z(i,k) * z(k,j)
            do 360 k = 1, l
               z(k,j) = z(k,j) - g * z(k,i)
  360    continue
  380    d(i) = z(i,i)
         z(i,i) = 1.0e0_sp
         if (l < 1) go to 500
         do 400 j = 1, l
            z(i,j) = 0.0e0_sp
            z(j,i) = 0.0e0_sp
  400    continue
  500 continue
      return
      end subroutine tred2

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)
      use cpmc
      implicit none

      integer::I,J,K,L,M,N,II,L1,L2,NM,MML,IERR
      REAL(sp),dimension(:)::D,E
      real(sp),dimension(:,:)::Z
      REAL(sp)::B,C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2
      REAL(sp)::PYTHAG

      IERR = 0
      IF (N == 1) GO TO 1001
      DO 100 I = 2, N
  100 E(I-1) = E(I)
      F = 0.0E0_sp
      B = 0.0E0_sp
      E(N) = 0.0E0_sp
      DO 240 L = 1, N
         J = 0
         H = ABS(D(L)) + ABS(E(L))
         IF (B < H) B = H
         DO 110 M = L, N
            IF (B + ABS(E(M)) == B) GO TO 120
  110    CONTINUE
  120    IF (M == L) GO TO 220
  130    IF (J == 30) GO TO 1000
         J = J + 1
         L1 = L + 1
         L2 = L1 + 1
         G = D(L)
         P = (D(L1) - G) / (2.0E0_sp * E(L))
         R = PYTHAG(P,1.0E0_sp)
         D(L) = E(L) / (P + SIGN(R,P))
         D(L1) = E(L) * (P + SIGN(R,P))
         DL1 = D(L1)
         H = G - D(L)
         IF (L2 > N) GO TO 145
         DO 140 I = L2, N
  140    D(I) = D(I) - H
  145    F = F + H
         P = D(M)
         C = 1.0E0_sp
         C2 = C
         EL1 = E(L1)
         S = 0.0E0_sp
         MML = M - L
         DO 200 II = 1, MML
            C3 = C2
            C2 = C
            S2 = S
            I = M - II
            G = C * E(I)
            H = C * P
            IF (ABS(P) < ABS(E(I))) GO TO 150
            C = E(I) / P
            R = SQRT(C*C+1.0E0_sp)
            E(I+1) = S * P * R
            S = C / R
            C = 1.0E0_sp / R
            GO TO 160
  150       C = P / E(I)
            R = SQRT(C*C+1.0E0_sp)
            E(I+1) = S * E(I) * R
            S = 1.0E0_sp / R
            C = C * S
  160       P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE
  200    CONTINUE
         P = -S * S2 * C3 * EL1 * E(L) / DL1
         E(L) = S * P
         D(L) = C * P
         IF (B + ABS(E(L)) > B) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
         DO 260 J = II, N
            IF (D(J) >= P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
         IF (K == I) GO TO 300
         D(K) = D(I)
         D(I) = P
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
  300 CONTINUE
      GO TO 1001
 1000 IERR = L
 1001 RETURN
      END subroutine tql2

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

      function pythag(a,b)
      use cpmc
      implicit none

      real(sp)::pythag
      real(sp)::a,b
      real(sp)::absa,absb
      absa=abs(a)
      absb=abs(b)

      if(absa>absb) then
        pythag=absa*sqrt(1.0_sp+(absb/absa)**2)
      else
        if(absb==0.0_sp) then
          pythag=0.0_sp
        else
          pythag=absb*sqrt(1.0_sp+(absa/absb)**2)
        end if
      end if

      return
      end function pythag

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

      SUBROUTINE MkInitOvlps
      use cpmc
      use jiekou,only:sgefa,sgedi

      integer::i,j,k,l,m,n
      integer::ipvt(NSTATES),info,ibas,jbas
      ! real(sp)::tmp_up(NUP,NUP),tmp_dn(NDN,NDN)
      real(sp)::det_up(NWFBAS),det_dn(NWFBAS)
      real(sp)::s(NWFBAS)
      real(sp)::det(2),work(NSTATES)

      do jbas=1,NIW
         do ibas=1,NWFBAS
            call dgemm('N', 'N', NUP, NUP, NSTATES, 1.0_8, phiZ_up(:,:,ibas), NUP, phiB_up(:,:,jbas), NSTATES, 0.0_8, tmp_up, NUP) !li
             
            !tmp_up=matmul(phiZ_up(:,:,ibas),phiB_up(:,:,jbas))
            !write(*,*)tmp_up

            call sgefa(tmp_up,NUP,NUP,ipvt,info)
	    if(info /= 0) stop 'Problem in sgefa routine'
            call sgedi(tmp_up,NUP,NUP,ipvt,det,work,10_i4b)
	    det_up(ibas)=det(1)*TEN**det(2)
            
            call dgemm('N', 'N', NDN, NDN, NSTATES, 1.0_8, phiZ_dn(:,:,ibas), NDN, phiB_dn(:,:,jbas), NSTATES, 0.0_8, tmp_dn, NDN) !li

            !tmp_dn=matmul(phiZ_dn(:,:,ibas),phiB_dn(:,:,jbas))

            call sgefa(tmp_dn,NDN,NDN,ipvt,info)
	    if(info /= 0) stop 'Problem in sgefa routine'
            call sgedi(tmp_dn,NDN,NDN,ipvt,det,work,10_i4b)
	    det_dn(ibas)=det(1)*TEN**det(2)

            s(ibas)=cwfbas(ibas)*det_up(ibas)*det_dn(ibas)

         end do
         ovlpINIT(jbas)=sum(s)
      end do

      do ibas=1,NIW
         sgnINIT(ibas)=sign(ONE,ovlpINIT(ibas))
         ovlpINIT(ibas)=abs(ovlpINIT(ibas))
      end do

      return
      END subroutine MkInitOvlps

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

      SUBROUTINE InitPop
      use cpmc
      integer::i,j,k,l,m,n
      integer::ibas,istart,iend,istop,iw
      real(sp)::cnorm,exch(NSTATES,1)

      cnorm=0.0
      do i=1,NIW
      cnorm=cnorm+abs(cwibas(i))
      end do

      istart=1
      do ibas=1,NIW
         istop=istart+nint(abs(cwibas(ibas)/cnorm)*NWLKRS)
         if(ibas == NIW) istop=NWLKRS
         iend=min(istop,NWLKRS)
         do iw=istart,iend

            do l=1,NSTATES

               do k=1,NUP
                  phi_up(l,k,iw)=phiB_up(l,k,ibas)
               end do

               do k=1,NDN
                  phi_dn(l,k,iw)=phiB_dn(l,k,ibas)
               end do

            end do

            if(sgnINIT(ibas)<0.0) then
            exch(:,1)=phi_up(:,1,iw)
            phi_up(:,1,iw)=phi_up(:,2,iw)
            phi_up(:,2,iw)=exch(:,1)
            end if

            ovlpDET(iw)=ONE
            wgtwlkr(iw)=ONE
            sgn(iw)=ONE

         end do

         istart=iend+1
      end do

      return
      END subroutine InitPop

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

      SUBROUTINE MkExpV
      use cpmc
      use jiekou,only:coupling
      integer::i,j,k,l,m,n
      real(sp)::hdeltau,sigma,alpha_u,tmpc1,tmpc2

      hdeltau=half*deltau
      do i=1,nsites
          if(hub_u(i)>-0.001) then
           call coupling(deltau*hub_u(i),alpha_u)
           do ising=-1,1,2
             tmpc1= alpha_u*ising-hdeltau*hub_u(i)
             tmpc2=-alpha_u*ising-hdeltau*hub_u(i)
             expV(i,ising,1)=exp( tmpc1 )
             expV(i,ising,2)=exp( tmpc2 )
             DeltaV(i,ising,1)=expV(i,ising,1)-1.0_sp
             DeltaV(i,ising,2)=expV(i,ising,2)-1.0_sp
             coeff(i,ising)=1.0
           end do
          else
           call coupling(-1.0*deltau*hub_u(i),alpha_u)
           do ising=-1,1,2
             tmpc1=alpha_u*ising-hdeltau*hub_u(i)
             tmpc2=alpha_u*ising-hdeltau*hub_u(i)
             expV(i,ising,1)=exp( tmpc1 )
             expV(i,ising,2)=exp( tmpc2 )
             DeltaV(i,ising,1)=expV(i,ising,1)-1.0_sp
             DeltaV(i,ising,2)=expV(i,ising,2)-1.0_sp
             coeff(i,ising)=exp(-1.0*alpha_u*ising+hdeltau*hub_u(i))
           end do
          end if
      end do

      do i=1,nsites
         if(vpd>-0.001) then
          call coupling(deltau*vpd,alpha_u)
          do ising=-1,1,2
            tmpc1= alpha_u*ising-hdeltau*vpd
            tmpc2=-alpha_u*ising-hdeltau*vpd
            eVpd(i,ising,1)=exp( tmpc1 )
            eVpd(i,ising,2)=exp( tmpc2 )
            DeVpd(i,ising,1)=eVpd(i,ising,1)-1.0_sp
            DeVpd(i,ising,2)=eVpd(i,ising,2)-1.0_sp
            coeffv(i,ising)=1.0
          end do
         else
          call coupling(-deltau*vpd,alpha_u)
          do ising=-1,1,2
            tmpc1=alpha_u*ising-hdeltau*vpd
            tmpc2=alpha_u*ising-hdeltau*vpd
            eVpd(i,ising,1)=exp( tmpc1 )
            eVpd(i,ising,2)=exp( tmpc2 )
            DeVpd(i,ising,1)=eVpd(i,ising,1)-1.0_sp
            DeVpd(i,ising,2)=eVpd(i,ising,2)-1.0_sp
            coeffv(i,ising)=exp(-1.0*alpha_u*ising+hdeltau*vpd)
          end do
         end if
      end do

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

      return
      end subroutine MkExpV

      !~~~~~Cosh(a) = Exp[ x/2 ]
      SUBROUTINE Coupling(x,a)
      use cpmc
      implicit none
      real(sp)::x,a,tmpxx
      integer::i,j,k,l,m,n

      tmpxx=exp(x/2.0)
      a=log( tmpxx+sqrt( tmpxx*tmpxx-1.0 ) )

      return
      end subroutine Coupling

      subroutine savephi
      use cpmc
      integer::iw,is,ip

        do 100 iw=1,NWLKRS
           phi_cup(:,:,iw) = phi_up(:,:,iw)
           phi_cdn(:,:,iw) = phi_dn(:,:,iw)
100     continue

      return
      end subroutine savephi

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

      SUBROUTINE Comb
      use cpmc
      use jiekou,only:PopCopy,ranGen
      integer::i,j,k,l,m,n
      real(sp)::wgtTOT,wgtAVG,offset
      real(sp)::psum(0:NWLKRS),wcomb(NWLKRS)
      integer::iw
      integer::imult(NWLKRS),iempty(NWLKRS)

     !! form the partial sum over the weights
      psum(0)=ZERO
      psum(1)=wgtwlkr(1)
      do iw=2,NWLKRS
         psum(iw)=psum(iw-1)+wgtwlkr(iw)
      end do

      wgtTOT=psum(NWLKRS)
      wgtAVG=wgtTOT/real(NWLKRS)

     !! make the comb and reset weights


      do iw=1,NWLKRS
         offset=ranGen(ISEED)*wgtAVG
         wgtwlkr(iw)=wgtAVG
         wcomb(iw)=offset+real(iw-1)*wgtAVG
      end do

      !! make list marking walkers to be replicated or eliminated
      k=1
      do iw=1,NWLKRS
         imult(iw)=0
 11      if(k <= NWLKRS) then
           if(psum(iw-1) < wcomb(k)) then
               if(wcomb(k) <= psum(iw)) then
                  k=k+1
                  imult(iw)=imult(iw)+1
                  goto 11
               end if
            end if
         end if
      end do

      !! create a list marking memory gaps
      k=0
      do iw=1,NWLKRS
         if(imult(iw) == 0) then
            k=k+1
            iempty(k)=iw
         end if
      end do

      !! replicate walkers and fill in memory
      k=0
      do j=1,NWLKRS
         if(imult(j) > 1) then
            do l=2,imult(j)
               k=k+1
               call PopCopy(iempty(k),j)
            end do
         end if
      end do

      return
      END subroutine Comb

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

      SUBROUTINE PopCopy(iwempty,iw)
      use cpmc
      integer::i,j,k,l,m,n
      integer::iwempty,iw,my

      phi_up(:,:,iwempty)=phi_up(:,:,iw)
      phi_dn(:,:,iwempty)=phi_dn(:,:,iw)

      ovlpDET(iwempty)=ovlpDET(iw)
      sgn(iwempty)=sgn(iw)

      if(measl==1.and.mstep>0) then

      phi_cup(:,:,iwempty)=phi_cup(:,:,iw)
      phi_cdn(:,:,iwempty)=phi_cdn(:,:,iw)


         do my=1,mstep       ! parent(mstep,iwempty)=parent(mstep,iw)
         do i=1,nsites
            kexpV_s(i,my,iwempty)=kexpV_s(i,my,iw)
         end do
         end do

         if(g_ph*w_ph>0.001) then
         do my=1,mstep       ! parent(mstep,iwempty)=parent(mstep,iw)
         do i=1,nsites
	    keph_s(i,my,iwempty)=keph_s(i,my,iw)
         end do
         end do
         end if

         if(abs(vpd)>0.001) then
         do my=1,mstep
         do i=1,lxy
         do m=1,2
         do n=1,4

           keVpd_s(i,my,iwempty,m,n)=keVpd_s(i,my,iw,m,n)

         end do
         end do
         end do
         end do
         end if

      end if

      return
      END subroutine PopCopy

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

      SUBROUTINE InitEnergy
      use cpmc
      use jiekou,only:sgefa,sgedi
      integer::i,j,k,l,m,n,ibas,jbas
      real(sp)::energy

      real(sp)::ekin,ecoul,detu,detd
      ! real(sp)::ovlpINV_up(NUP,NUP),tmp_up(NSTATES,NUP)
      ! real(sp)::ovlpINV_dn(NDN,NDN),tmp_dn(NSTATES,NDN)
      real(sp)::pnew_up(NUP,NUP),pnew_dn(NDN,NDN)
      real(sp)::pg_up(NUP,NUP),pg_dn(NDN,NDN)
      real(sp)::detpbas(NIW)
      ! real(sp)::gx_up(NSTATES,NSTATES),gx_dn(NSTATES,NSTATES)
      real(sp)::ewlkr(NIW),ebas(NIW),detwlkr(NIW)
      real(sp)::det(2),work(NSTATES)
      integer::info,ipvt(NSTATES)

      !g = <a^{\dagger}a> = I-G = [R(LR)^{-1}L]/det(LR) /=0!
      !g = computing adjoin matrix det(LR) =0 or /=0!

      estptop=ZERO
      estpbtm=ZERO
      ewlkr=0.0
      detwlkr=0.0
      do 600 ibas=1,NIW
        do 700 jbas=1,NIW

         ebas(jbas)=ZERO

         !----UP----compute overlap matrix LR, etc. for spin up!

         call dgemm('T', 'N', NUP, NUP, NSTATES, 1.0_8, phiB_up(:,:,ibas), NSTATES, phiB_up(:,:,jbas), NSTATES, 0.0_8, ovlpINV_up, NUP) !li
         !write(*,*) ovlpINV_up
         !write(*,*) '858'

 
         !ovlpINV_up=matmul(transpose(phiB_up(:,:,ibas)),phiB_up(:,:,jbas))
         !write(*,*) ovlpINV_up
         !write(*,*) '865'

         pnew_up=ovlpINV_up

         !compute determinant of overlap matrix LR!
         call sgefa(pnew_up,NUP,NUP,ipvt,info)
         if(info /= 0) then
         detp_up(jbas)=0.0
         else
         call sgedi(pnew_up,NUP,NUP,ipvt,det,work,10_i4b)
         detp_up(jbas)=det(1)*TEN**det(2)
         end if

         !compute g(R[(LR)^{-1}]L)!
         pg_up=0.0
         do i=1,NUP
           do j=1,NUP
           pnew_up=ovlpINV_up
           pnew_up(j,:)=0.0; pnew_up(:,i)=0.0; pnew_up(j,i)=1.0
           call sgefa(pnew_up,NUP,NUP,ipvt,info)
           if(info /= 0) then
            detu=0.0
           else
            call sgedi(pnew_up,NUP,NUP,ipvt,det,work,10_i4b)
            detu=det(1)*TEN**det(2)
           end if
           pg_up(i,j)=detu
           end do
         end do

         call dgemm('N', 'N', NSTATES, NUP, NUP, 1.0_8, phiB_up(:,:,jbas), NSTATES, pg_up, NUP, 0.0_8, tmp_up, NSTATES) !li

         !tmp_up=matmul(phiB_up(:,:,jbas),pg_up)

         call dgemm('N', 'T', NSTATES, NSTATES, NUP, -1.0_8, tmp_up, NSTATES, phiB_up(:,:,ibas), NSTATES, 0.0_8, gx_up(:,:), NSTATES)!li
         
         !gx_up(:,:)=-matmul(tmp_up,transpose(phiB_up(:,:,ibas)))
         !write(*,*)'904'
         !write(*,*)gx_up(:,:)

         do i=1,NSTATES
         gx_up(i,i)=detp_up(jbas)+gx_up(i,i)
         end do

         !---DOWN---compute overlap matrix LR, etc. for spin down!

         call dgemm('T', 'N', NDN, NDN, NSTATES, 1.0_8, phiB_dn(:,:,ibas), NSTATES, phiB_dn(:,:,jbas), NSTATES, 0.0_8, ovlpINV_dn, NDN)!li check

         !ovlpINV_dn=matmul(transpose(phiB_dn(:,:,ibas)),phiB_dn(:,:,jbas))
         pnew_dn=ovlpINV_dn

         !compute inverse and determinant of overlap matrix LR!
         call sgefa(pnew_dn,NDN,NDN,ipvt,info)
         if(info /= 0) then
          detp_dn(jbas)=0.0
         else
          call sgedi(pnew_dn,NDN,NDN,ipvt,det,work,10_i4b)
          detp_dn(jbas)=det(1)*TEN**det(2)
         end if

         !compute g (R[(LR)^1]L)!
         pg_dn=0.0
         do i=1,NDN
           do j=1,NDN
           pnew_dn=ovlpINV_dn
           pnew_dn(j,:)=0.0; pnew_dn(:,i)=0.0; pnew_dn(j,i)=1.0
           call sgefa(pnew_dn,NDN,NDN,ipvt,info)
           if(info /= 0) then
            detd=0.0
           else
            call sgedi(pnew_dn,NDN,NDN,ipvt,det,work,10_i4b)
            detd=det(1)*TEN**det(2)
           end if
           pg_dn(i,j)=detd
           end do
         end do

         call dgemm('N', 'N', NSTATES, NDN, NDN, 1.0_8, phiB_dn(:,:,jbas), NSTATES, pg_dn, NDN, 0.0_8, tmp_dn, NSTATES) !li
          
         !tmp_dn=matmul(phiB_dn(:,:,jbas),pg_dn)
    
         call dgemm('N', 'T', NSTATES, NSTATES, NDN, -1.0_8, tmp_dn, NSTATES, phiB_dn(:,:,ibas), NSTATES, 0.0_8, gx_dn(:,:), NSTATES) !li
          
         !gx_dn(:,:)=-matmul(tmp_dn,transpose(phiB_dn(:,:,ibas))) 
         
         

         do i=1,NSTATES
          gx_dn(i,i)=detp_dn(jbas)+gx_dn(i,i)
         end do

         detpbas(jbas)=cwibas(jbas)*detp_up(jbas)*detp_dn(jbas)

         !**compute one-body contribution to energy**!
         ekin=ZERO
         do j=1,NSTATES
            do i=1,NSTATES
!---------------HKAnHop-------------------!
            ekin=ekin-tk(j,i,1)*gx_up(i,j)*detp_dn(jbas)-&
                    & tk(j,i,2)*gx_dn(i,j)*detp_up(jbas)
!---------------HKAnHop-------------------!
            end do
         end do
         ebas(jbas)=ebas(jbas)+ekin

         !**compute two-body-contributions to energy**!
         ecoul=ZERO
         do i=1,NSTATES
            ecoul=ecoul+hub_u(i)*gx_up(i,i)*gx_dn(i,i)-&
                    &hub_u(i)*detp_up(jbas)*gx_dn(i,i)-&
                    &hub_u(i)*detp_dn(jbas)*gx_up(i,i)+&
                    &hub_u(i)*detp_up(jbas)*detp_dn(jbas)
         end do
         ebas(jbas)=ebas(jbas)+ecoul

         ebas(jbas)=ebas(jbas)*cwibas(jbas)

700      continue
         detwlkr(ibas)=cwibas(ibas)*sum(detpbas)
         ewlkr(ibas)=cwibas(ibas)*sum(ebas)
600   continue

      estptop=sum(ewlkr)
      estpbtm=sum(detwlkr)

      etrial=estptop/(estpbtm)
      e_var=etrial
      write(*,*)estptop,estpbtm,etrial

      write(*,*)'Initial Energy=',etrial
      write(*,*)'Initial Energy Per N=',etrial/lxy
      write(*,*)'###########################################'
      write(*,*)

      return
      END subroutine InitEnergy

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

      SUBROUTINE EstEtrial
      use cpmc
      use jiekou,only:Step,Stblz,comb,Step1,Stblz1
      integer::i,j,k,l,m,n
      integer::istp,iw
      real(sp)::egrowth,w_up,w_down,wgtsum

      write(*,*) ' '
      write(*,*) 'Growth Estimation:'
      write(*,'((a),t11,(a),t28,(a),t38,(a))')&
            & ' Iblk','Egrowth', 'W_Up', 'W_Down'

      egrowth=ZERO
      do i=1,nblkgr*8
         wgtwlkr=(real(NWLKRS)/sum(wgtwlkr))*wgtwlkr
         w_up=sum(wgtwlkr)

         do istp=1,itvlorth
            call Step1(istp)
         end do

         w_down=sum(wgtwlkr)
         
         call Stblz1
         
         call comb

         egrowth=egrowth+w_up/w_down
         write(*,'(i5,3e12.4)') i,egrowth,w_up,w_down

      end do

      etrial=etrial+log(egrowth/real(nblkgr*8))/deltau/real(itvlorth)

      write(*,*)
      write(*,*)'Growth estimate of Etrial: ',etrial
      write(*,*)'Growth estimate of Etrial Per N: ',etrial/lxy
      write(*,*)

      return
      end subroutine EstEtrial

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

      SUBROUTINE StpMeas
      use cpmc
      use jiekou,only:sgefa,sgedi,ranGen
      integer::i,j,k,l,m,n
      integer::iw,ibas

      real(sp)::ekin,ecoul,ovlpNEW
      ! real(sp)::ovlpINV_up(NUP,NUP),tmp_up(NSTATES,NUP),det_up(NWFBAS)
      ! real(sp)::ovlpINV_dn(NDN,NDN),tmp_dn(NSTATES,NDN),det_dn(NWFBAS)
      real(sp)::detpbas(NWFBAS)
      ! real(sp)::gx_up(NSTATES,NSTATES),gx_dn(NSTATES,NSTATES)
      real(sp)::ewlkr(NWLKRS),ebas(NWFBAS),detwlkr(NWLKRS)
      real(sp)::det(2),work(NSTATES)
      integer::kp
      integer::info,ipvt(NSTATES)

      !! g = <a^{\dagger}a> = I-G = [R(LR)^{-1}L]/det(LR)

      ewlkr=0.0
      detwlkr=0.0
      do 600 iw=1,NWLKRS

        if(wgtwlkr(iw) /= ZERO) then

        do 700 ibas=1,NWFBAS
         ebas(ibas)=ZERO

         !!--UP--compute overlap matrix LR, etc. for spin up


         call dgemm('N', 'N', NUP, NUP, NSTATES, 1.0_8, phiZ_up(:,:,ibas), NUP, phi_up(:,:,iw), NSTATES, 0.0_8, ovlpINV_up, NUP) !li 11.2

         !ovlpINV_up=matmul(phiZ_up(:,:,ibas),phi_up(:,:,iw))

         !compute inverse and determinant of overlap matrix LR
         call sgefa(ovlpINV_up,NUP,NUP,ipvt,info)
         if(info /= 0) stop 'Problem in sgefa routine'
         call sgedi(ovlpINV_up,NUP,NUP,ipvt,det,work,11_i4b)
         detp_up(ibas)=det(1)*TEN**det(2)

         !! compute g(R[(LR)^{-1}]L)



         call dgemm('N', 'N', NSTATES, NUP, NUP, 1.0_8, phi_up(:,:,iw), NSTATES, ovlpINV_up, NUP, 0.0_8, tmp_up, NSTATES) !li

         !tmp_up=matmul(phi_up(:,:,iw),ovlpINV_up)


         call dgemm('N', 'N', NSTATES, NSTATES, NUP, 1.0_8, tmp_up, NSTATES, phiZ_up(:,:,ibas), NUP, 0.0_8, gx_up(:,:), NSTATES) !cao

         !gx_up(:,:)=matmul(tmp_up,phiZ_up(:,:,ibas))

         !!--DOWN--compute overlap matrix LR, etc. for spin down

         call dgemm('N', 'N', NDN, NDN, NSTATES, 1.0_8, phiZ_dn(:,:,ibas), NDN, phi_dn(:,:,iw), NSTATES, 0.0_8, ovlpINV_dn, NDN) !cao

         !ovlpINV_dn=matmul(phiZ_dn(:,:,ibas),phi_dn(:,:,iw))

         !compute inverse and determinant of overlap matrix LR
         call sgefa(ovlpINV_dn,NDN,NDN,ipvt,info)
         if(info /= 0) stop 'Problem in sgefa routine'
         call sgedi(ovlpINV_dn,NDN,NDN,ipvt,det,work,11_i4b)
         detp_dn(ibas)=det(1)*TEN**det(2)

         !! compute g (R[(LR)^1]L)

         call dgemm('N', 'N', NSTATES, NDN, NDN, 1.0_8, phi_dn(:,:,iw), NSTATES, ovlpINV_dn, NDN, 0.0_8, tmp_dn, NSTATES) !li

         !tmp_dn=matmul(phi_dn(:,:,iw),ovlpINV_dn)

         call dgemm('N', 'N', NSTATES, NSTATES, NDN, 1.0_8, tmp_dn, NSTATES, phiZ_dn(:,:,ibas), NDN, 0.0_8, gx_dn, NSTATES) !cao

         !gx_dn(:,:)=matmul(tmp_dn,phiZ_dn(:,:,ibas))

         detpbas(ibas)=cwfbas(ibas)*detp_up(ibas)*detp_dn(ibas)

         !! compute one-body contribution to energy
         ekin=ZERO
         do j=1,NSTATES
           do i=1,NSTATES
           ekin=ekin+tk(j,i,1)*gx_up(i,j)+tk(j,i,2)*gx_dn(i,j)
           end do
         end do
         ebas(ibas)=ebas(ibas)+ekin


         !! compute two-body-contributions to energy
         ecoul=ZERO
         do i=1,NSTATES
           ecoul=ecoul+hub_u(i)*gx_up(i,i)*gx_dn(i,i)
         end do
         ebas(ibas)=ebas(ibas)+ecoul

         if(abs(vpd)>0.001) then
         ecoul=ZERO
         do i=1,nsites
         do j1=1,2
            j=idis(i,j1)
            ecoul=ecoul+vpd*(gx_up(i,i)*gx_up(j,j)-gx_up(i,j)*gx_up(j,i))&
                     & +vpd*(gx_dn(i,i)*gx_dn(j,j)-gx_dn(i,j)*gx_dn(j,i))&
                     & +vpd*(gx_up(i,i)*gx_dn(j,j)+gx_up(j,j)*gx_dn(i,i))
         end do
         end do
         ebas(ibas)=ebas(ibas)+ecoul
         end if

         ebas(ibas)=ebas(ibas)*detpbas(ibas)

700      continue

         detwlkr(iw)=wgtwlkr(iw)*sum(detpbas)/ovlpDET(iw)
         ewlkr(iw)=wgtwlkr(iw)*sum(ebas)/ovlpDET(iw)

         end if

600   continue

      estptop=estptop+sum(ewlkr)
      estpbtm=estpbtm+sum(detwlkr)
      nstp=nstp+1

      return
      END subroutine StpMeas

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    SUBROUTINE preprog
      use cpmc
      use jiekou,only:ranGen
      integer::iw,i,j,k

    do iw=1,NWLKRS
      do i=1,NSTATES
      if(ranGen(ISEED)>0.5) then
        ising_ph(i,iw)=1
      else
        ising_ph(i,iw)=-1
      end if
      end do
    end do

    return
    end SUBROUTINE preprog

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    !**************************************************!

      SUBROUTINE Step(istp)
      use cpmc
      use jiekou,only:HalfT,FullV
      integer::i,j,k,l,m,n,iselect
      integer::iw,istp,mx
      real(sp)::hdeltau,alpha_u

      iselect=mod(istp,itvlorth)

!------------------HoKinPara----------------!
!      do 100 iw=1,NWLKRS
      do 100 iw=iwStart(myID),iwEnd(myID)   
!------------------HoKinPara----------------!

       !! SKIP WALKERS WITH wgt=0
         if(wgtwlkr(iw) <= ZERO) go to 100

       !! APPLY growth control factor
         wgtwlkr(iw)=exp(deltau*etrial)*wgtwlkr(iw)

       !! APPLY EXP(-0.5*DelTau*T):

         if(iselect==1) then
         call HalfT(iw,1)
         else
         call halfT(iw,2)
         end if

         if(wgtwlkr(iw) <= ZERO) go to 100

       !! APPLY exp(-DelTau*V):
       !-------Zhongbing-phonon fields------!
         if(g_ph*w_ph>0.001) then
         call Fullph(iw)
         if(wgtwlkr(iw) <= ZERO) go to 100
         end if
       !------------------------------------!
         call FullV(iw)
         if(wgtwlkr(iw) <= ZERO) go to 100

         if(abs(vpd)>0.001) then
         call V_upup(iw)
         if(wgtwlkr(iw) <=ZERO) goto 100
         call V_updn(iw)
         if(wgtwlkr(iw) <=ZERO) goto 100
         call V_dnup(iw)
         if(wgtwlkr(iw) <=ZERO) goto 100
         call V_dndn(iw)
         if(wgtwlkr(iw) <=ZERO) goto 100
         end if

         mx=mod(istp,itvl_m)
         if(mx==0) mx=itvl_m

         if(measl==1) then

          hdeltau=half*deltau
          do i=1,nsites             ! Saved V(istp,iw)
            kexpV_s(i,mx,iw)=ising_u(i)
          end do

          if(g_ph*w_ph>0.001) then
          do i=1,nsites             ! Saved phonon fields
            keph_s(i,mx,iw)=ising_sum(i)
          end do
          end if

          if(abs(vpd)>0.001) then
          do i=1,lxy               ! Saved Vpd(istp,iw)
          do m=1,2
          do n=1,4

            keVpd_s(i,mx,iw,m,n)=ising_v(i,m,n)

          end do
          end do
          end do
          end if

         end if

       !! APPLY exp(-0.5*DelTau*T):

         if(iselect==0) call HalfT(iw,1)

100    continue

      return
      END subroutine Step

      !**************************************************!

      SUBROUTINE HalfT(iw,icase)
      use cpmc
      use jiekou,only:Ovlps

      integer::iw,icase
      integer::i,j,k,l,m,n
      ! real(sp)::tmp_up(NSTATES,NUP),tmp_dn(NSTATES,NDN)

      !! Advance phi by exp(-deltau*T/2) for each spin

      if(icase==1) then

      call dgemm('N', 'N', NSTATES, NUP, NSTATES, 1.0_8, expT(:,:,1), NSTATES, phi_up(:,:,iw), NSTATES, 0.0_8, tmp_up, NSTATES) !cao 

      !tmp_up=matmul(expT(:,:,1),phi_up(:,:,iw))
      phi_up(:,:,iw)=tmp_up

      call dgemm('N', 'N', NSTATES, NDN, NSTATES, 1.0_8, expT(:,:,2), NSTATES, phi_dn(:,:,iw), NSTATES, 0.0_8, tmp_dn, NSTATES) !cao

      !tmp_dn=matmul(expT(:,:,2),phi_dn(:,:,iw))
      phi_dn(:,:,iw)=tmp_dn

      else
      
      call dgemm('N', 'N', NSTATES, NUP, NSTATES, 1.0_8, exp2T(:,:,1), NSTATES, phi_up(:,:,iw), NSTATES, 0.0_8, tmp_up, NSTATES) !cao 

      !tmp_up=matmul(exp2T(:,:,1),phi_up(:,:,iw))
      phi_up(:,:,iw)=tmp_up

      call dgemm('N', 'N', NSTATES, NDN, NSTATES, 1.0_8, exp2T(:,:,2), NSTATES, phi_dn(:,:,iw), NSTATES, 0.0_8, tmp_dn, NSTATES) !cao 

      !tmp_dn=matmul(exp2T(:,:,2),phi_dn(:,:,iw))
      phi_dn(:,:,iw)=tmp_dn

      end if

      !! Compute determinant of the overlap integral, and Inverse

      call Ovlps(iw)

      return
      END subroutine HalfT

      !**************************************************!

      SUBROUTINE Ovlps(iw)
      use cpmc
      use jiekou,only:sgefa,sgedi
      integer::i,j,k,l,m,n
      integer::iw,ibas,info

      real(sp)::ovlpNEW
      real(sp)::det(2),work(NSTATES)
      ! real(sp)::ovlpINV_up(NUP,NUP),ovlpINV_dn(NDN,NDN)
      real(sp)::det_up(NWFBAS),det_dn(NWFBAS)
      integer::ipvt(NSTATES)

      !! Compute determinant of the overlap integral, and Inverse
      do ibas=1,NWFBAS

         call dgemm('N', 'N', NUP, NUP, NSTATES, 1.0_8, phiZ_up(:,:,ibas), NUP, phi_up(:,:,iw), NSTATES, 0.0_8, ovlpINV_up, NUP) !cao  

         !ovlpINV_up=matmul(phiZ_up(:,:,ibas),phi_up(:,:,iw))

         call sgefa(ovlpINV_up,NUP,NUP,ipvt,info)
         if(info /= 0) stop 'Problem in sgefa routine'
         call sgedi(ovlpINV_up,NUP,NUP,ipvt,det,work,11_i4b)
         det_up(ibas)=det(1)*TEN**det(2)

       !! Rename G to be the inverse,
          g_up(:,:,ibas)=ovlpINV_up

          call dgemm('N', 'N', NDN, NDN, NSTATES, 1.0_8, phiZ_dn(:,:,ibas), NDN, phi_dn(:,:,iw), NSTATES, 0.0_8, ovlpINV_dn, NDN) !cao  

         !ovlpINV_dn=matmul(phiZ_dn(:,:,ibas),phi_dn(:,:,iw))

         call sgefa(ovlpINV_dn,NDN,NDN,ipvt,info)
         if(info /= 0) stop 'Problem in sgefa routine'
         call sgedi(ovlpINV_dn,NDN,NDN,ipvt,det,work,11_i4b)
         det_dn(ibas)=det(1)*TEN**det(2)

       !! Rename G to be the inverse,
          g_dn(:,:,ibas)=ovlpINV_dn

          detbas(ibas)=cwfbas(ibas)*det_up(ibas)*det_dn(ibas)

      end do

      ovlpNEW=sum(detbas)

      if(ovlpNEW*sgn(iw) <= ZERO) then
         wgtwlkr(iw)=ZERO
      else
         wgtwlkr(iw)=(ovlpNEW/ovlpDET(iw))*wgtwlkr(iw)
         ovlpDET(iw)=ovlpNEW
      end if

      return
      END subroutine Ovlps

      !**************************************************!
      !--------------Zhongbing-phonon fields-------------!
      SUBROUTINE Fullph(iw)
      use cpmc
      use jiekou, only:ranGen
      integer::i,j,k,l,m,n
      integer::iw,itmph

      integer::ibas
      real(sp)::ovlpNEW,ovlpOLD
      real(sp)::tmpup,tmpdn
      real(sp)::rdet(NWFBAS,-1:1,NSPIN),pt(-1:1),ptsum,ptest
      real(sp)::gl_up(NUP,NSTATES,NWFBAS),gr_up(NSTATES,NUP,NWFBAS)
      real(sp)::gl_dn(NDN,NSTATES,NWFBAS),gr_dn(NSTATES,NDN,NWFBAS)
      real(sp)::gup_ii(NWFBAS),gdn_ii(NWFBAS),r_detbas(NWFBAS)
      real(sp)::ggx1_up(1,NUP),ggx2_up(NUP,1)
      real(sp)::ggx1_dn(1,NDN),ggx2_dn(NDN,1)

      real(sp)::ggx1_up1(1,NUP),ggx2_up2(NUP,1)
      real(sp)::ggx1_dn1(1,NDN),ggx2_dn2(NDN,1)

      !!*****call MkExpV             ! Make e^V for each spin

      ovlpOLD=ovlpDET(iw)

      !Get Ising Variables, according to important sampling

   do 100 i=1,nsites

      !For fixed i; up & down; N^2
      !Make gl_s(k,i) = [ PhiT^+(ibas) Phi ]^-1(k,j) * PhiT^+(j,i); Ns x NSTATES
      !!!!!!gr_s(i,k) = Phi(i,j) * [ PhiT^+(ibas) Phi ]^-1(j,k);    NSTATES x Ns

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

      !**************************************************!
      SUBROUTINE FullV(iw)
      use cpmc
      use jiekou,only:ranGen
      integer::i,j,k,l,m,n
      integer::iw

      integer::ibas
      real(sp)::ovlpNEW,ovlpOLD
      real(sp)::tmpup,tmpdn
      !!real(sp)::s_detbas(NWFBAS),s_wgtwlkr,s_ovlpDET         ! Later
      real(sp)::rdet(NWFBAS,-1:1,NSPIN),pt(-1:1),ptsum,ptest
      real(sp)::gl_up(NUP,NSTATES,NWFBAS),gr_up(NSTATES,NUP,NWFBAS)
      real(sp)::gl_dn(NDN,NSTATES,NWFBAS),gr_dn(NSTATES,NDN,NWFBAS)
      real(sp)::gup_ii(NWFBAS),gdn_ii(NWFBAS),r_detbas(NWFBAS)
      real(sp)::ggx1_up(1,NUP),ggx2_up(NUP,1)
      real(sp)::ggx1_dn(1,NDN),ggx2_dn(NDN,1)

      real(sp)::ggx1_up1(1,NUP),ggx2_up2(NUP,1)
      real(sp)::ggx1_dn1(1,NDN),ggx2_dn2(NDN,1)

      !!*****call MkExpV             ! Make e^V for each spin

      ovlpOLD=ovlpDET(iw)

      !Get Ising Variables, according to important sampling

   do 100 i=1,nsites

      !For fixed i; up & down; N^2
      !Make gl_s(k,i) = [ PhiT^+(ibas) Phi ]^-1(k,j) * PhiT^+(j,i); Ns x NSTATES
      !!!!!!gr_s(i,k) = Phi(i,j) * [ PhiT^+(ibas) Phi ]^-1(j,k);    NSTATES x Ns

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
           rdet(ibas,ising,1)=1.0+ DeltaV(i,ising,1)*gup_ii(ibas)
           rdet(ibas,ising,2)=1.0+ DeltaV(i,ising,2)*gdn_ii(ibas)

           r_detbas(ibas)=rdet(ibas,ising,1)*rdet(ibas,ising,2)*detbas(ibas)
           !! cwfbas(ibas) is included(ovlps)

          end do
          ovlpNEW=sum(r_detbas)

          pt(ising)=coeff(i,ising)*ovlpNEW/ovlpOLD/two !original probability
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
          ising_u(i)= 1
          ovlpNEW=pt(1)*ovlpOLD*two/coeff(i,ising_u(i)) !see the above line
         else
          ising_u(i)=-1
          ovlpNEW=pt(-1)*ovlpOLD*two/coeff(i,ising_u(i))!see the above line
         end if

         itmp=ising_u(i)                    ! picked Ising_u(i)

         !! New Wave function, only i_th row are changed
         !! advance phi by exp(-deltau*V) for picked Ising variables
         tmpup=expV(i,itmp,1)
         tmpdn=expV(i,itmp,2)

         phi_up(i,:,iw)=tmpup*phi_up(i,:,iw)
         phi_dn(i,:,iw)=tmpdn*phi_dn(i,:,iw)

         !! Updating Green's(?) functions due to change of ising_u(i)
      do 200 ibas=1,NWFBAS

         tmpup=DeltaV( i,itmp,1) /rdet( ibas,itmp,1 )

         ggx2_up(:,1)=gl_up(:,i,ibas)
         ggx1_up(1,:)=gr_up(i,:,ibas)
         g_up(:,:,ibas)=g_up(:,:,ibas)-tmpup*matmul(ggx2_up,ggx1_up)

         tmpdn=DeltaV( i,itmp,2 )/rdet( ibas,itmp,2 )

         ggx2_dn(:,1)=gl_dn(:,i,ibas)
         ggx1_dn(1,:)=gr_dn(i,:,ibas)
         g_dn(:,:,ibas)=g_dn(:,:,ibas)-tmpdn*matmul(ggx2_dn,ggx1_dn)

      200 continue

         !! compute determinant of the overlap integral
         do ibas=1,NWFBAS
          detbas(ibas)=rdet(ibas,itmp,1)*rdet(ibas,itmp,2)*detbas(ibas)
         end do

         ovlpOLD=ovlpNEW                ! ready for next move
         wgtwlkr(iw)=ptsum*wgtwlkr(iw)  ! renormalized walker weight

  100   continue

        ovlpDET(iw)=ovlpOLD

        return
        END subroutine FullV

      !*********************************************************!

      SUBROUTINE V_upup(iw)
      use cpmc
      use jiekou
      integer::i,j,k,l,m,n,j1
      integer::iw

      integer::ibas
      real(sp)::ovlpNEW,ovlpOLD
      real(sp)::tmpup,tmpup1
      !!real(sp)::s_detbas(NWFBAS),s_wgtwlkr,s_ovlpDET         ! Later
      real(sp)::rdet(NWFBAS,-1:1,NSPIN),pt(-1:1),ptsum,ptest
      real(sp)::gl_up(NUP,NSTATES,NWFBAS),gr_up(NSTATES,NUP,NWFBAS)
      real(sp)::gl1_up(NUP,NSTATES,NWFBAS),gr1_up(NSTATES,NUP,NWFBAS)
      real(sp)::gup_ii(NWFBAS),r_detbas(NWFBAS)
      real(sp)::gup_jj(NWFBAS)
      real(sp)::ggx1_up(1,NUP),ggx2_up(NUP,1)
      real(sp)::gt_up(NUP,NUP,NWFBAS,-1:1)

      real(sp)::ggx1_up1(1,NUP),ggx2_up2(NUP,1)

      ovlpOLD=ovlpDET(iw)

      !!Get Ising Variables, according to important sampling

   do 100 i=1,nsites

   do 100 j1=1,2

         j=idis(i,j1)

      !For fixed i; up & down; N^2
      !Make gl_s(k,i) = [ PhiT^+(ibas) Phi ]^-1(k,j) * PhiT^+(j,i); Ns x NSTATES
      !gr_s(i,k) = Phi(i,j) * [ PhiT^+(ibas) Phi ]^-1(j,k);    NSTATES x Ns

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

      120 continue

      do 140 ibas=1,NWFBAS

         gup_ii(ibas)=dot_product(gr_up(i,:,ibas),phiT_up(i,:,ibas))

      140 continue

      !!check x(i), p( x(i) )

      do ising=-1,1,2

         do ibas=1,NWFBAS

            rdet(ibas,ising,1)=1.0+ DeVpd(i,ising,1)*gup_ii(ibas)

            tmpup=DeVpd(i,ising,1)/rdet(ibas,ising,1)
            ggx1_up(1,:)=gr_up(i,:,ibas)
            ggx2_up(:,1)=gl_up(:,i,ibas)

            gt_up(:,:,ibas,ising)=g_up(:,:,ibas)-tmpup*matmul(ggx2_up,ggx1_up)

            ggx2_up(:,1)=phiT_up(j,:,ibas)

            call dgemm('N', 'N', NUP, 1, NUP, 1.0_8, gt_up(:,:,ibas,ising), NUP, ggx2_up, NUP, 0.0_8, ggx2_up2, NUP) !cao
            ggx2_up=ggx2_up2

            !ggx2_up=matmul(gt_up(:,:,ibas,ising),ggx2_up)
            gl1_up(:,j,ibas)=ggx2_up(:,1)

            ggx1_up(1,:)=phi_up(j,:,iw)

            call dgemm('N', 'N', 1, NUP, NUP, 1.0_8, ggx1_up, 1, gt_up(:,:,ibas,ising), NUP, 0.0_8, ggx1_up1, 1) !cao
            ggx1_up=ggx1_up1

            !ggx1_up=matmul(ggx1_up,gt_up(:,:,ibas,ising))
            gr1_up(j,:,ibas)=ggx1_up(1,:)

            gup_jj(ibas)=dot_product(gr1_up(j,:,ibas),phiT_up(j,:,ibas))

            rdet(ibas,ising,2)=1.0+DeVpd(i,ising,2)*gup_jj(ibas)

         end do

         do ibas=1,NWFBAS
            r_detbas(ibas)=rdet(ibas,ising,1)*rdet(ibas,ising,2)*detbas(ibas)
         end do

         ovlpNEW=sum(r_detbas)

         pt(ising)=coeffv(i,ising)*ovlpNEW/ovlpOLD/two  !original probability
         if(ovlpNEW*sgn(iw) <= ZERO) then
           wgtwlkr(iw)=wgtwlkr(iw)/( 1.0-pt(ising) )
           pt(ising)=zero
         end if

      end do

      !!calculate p(1)/( p(1)+p(-1) )

      ptsum=pt(1)+pt(-1)

      if(ptsum == zero) then
        wgtwlkr(iw)=ZERO     ! walker is terminated, nothing left
        return
      end if

      ptest=pt(1)/ptsum
      if( ptest > ranGen(ISEED) ) then
        ising_v(i,j1,1)= 1
        ovlpNEW=pt(1)*ovlpOLD*two/coeffv(i,1)  !see the above line
      else
        ising_v(i,j1,1)=-1
        ovlpNEW=pt(-1)*ovlpOLD*two/coeffv(i,-1) !see the above line
      end if

      itmp=ising_v(i,j1,1)           ! picked Ising_u(i)

      !!New Wave function, only i_th row are changed
      !!advance phi by exp(-deltau*V) for picked Ising variables

      tmpup=eVpd(i,itmp,1)
      tmpup1=eVpd(i,itmp,2)

      phi_up(i,:,iw)=tmpup*phi_up(i,:,iw)
      phi_up(j,:,iw)=tmpup1*phi_up(j,:,iw)

      !!Updating Green's(?) functions due to change of ising_v(i)
      do 200 ibas=1,NWFBAS

         tmpup=DeVpd( i,itmp,2) /rdet( ibas,itmp,2 )

         ggx2_up(:,1)=gl1_up(:,j,ibas)
         ggx1_up(1,:)=gr1_up(j,:,ibas)
         g_up(:,:,ibas)=gt_up(:,:,ibas,itmp)-tmpup*matmul(ggx2_up,ggx1_up)

      200 continue

      !!compute determinant of the overlap integral
      do ibas=1,NWFBAS
         detbas(ibas)=rdet(ibas,itmp,1)*rdet(ibas,itmp,2)*detbas(ibas)
      end do

      ovlpOLD=ovlpNEW                ! ready for next move
      wgtwlkr(iw)=ptsum*wgtwlkr(iw)  ! renormalized walker weight

  100 continue

      ovlpDET(iw)=ovlpOLD

      return
      END subroutine V_upup

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

      SUBROUTINE V_updn(iw)
      use cpmc
      use jiekou
      integer::i,j,k,l,m,n,j1
      integer::iw

      integer::ibas
      real(sp)::ovlpNEW,ovlpOLD
      real(sp)::tmpup,tmpdn
      !!real(sp)::s_detbas(NWFBAS),s_wgtwlkr,s_ovlpDET         ! Later
      real(sp)::rdet(NWFBAS,-1:1,NSPIN),pt(-1:1),ptsum,ptest
      real(sp)::gl_up(NUP,NSTATES,NWFBAS),gr_up(NSTATES,NUP,NWFBAS)
      real(sp)::gl_dn(NDN,NSTATES,NWFBAS),gr_dn(NSTATES,NDN,NWFBAS)
      real(sp)::gup_ii(NWFBAS),gdn_ii(NWFBAS),r_detbas(NWFBAS)
      real(sp)::ggx1_up(1,NUP),ggx2_up(NUP,1)
      real(sp)::ggx1_dn(1,NDN),ggx2_dn(NDN,1)

      real(sp)::ggx1_up1(1,NUP),ggx2_up2(NUP,1)
      real(sp)::ggx1_dn1(1,NDN),ggx2_dn2(NDN,1)

      ovlpOLD=ovlpDET(iw)

      !!Get Ising Variables, according to important sampling

   do 100 i=1,nsites

   do 100 j1=1,2

      j=idis(i,j1)

      !!For fixed i; up & down; N^2
      !!Make gl_s(k,i) = [ PhiT^+(ibas) Phi ]^-1(k,j) * PhiT^+(j,i); Ns x NSTATES
      !!gr_s(i,k) = Phi(i,j) * [ PhiT^+(ibas) Phi ]^-1(j,k);    NSTATES x Ns

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

         ggx2_dn(:,1)=phiT_dn(j,:,ibas)

         call dgemm('N', 'N', NDN, 1, NDN, 1.0_8, g_dn(:,:,ibas), NDN, ggx2_dn, NDN, 0.0_8, ggx2_dn2, NDN) !cao
         ggx2_dn=ggx2_dn2

         !ggx2_dn=matmul(g_dn(:,:,ibas),ggx2_dn)
         gl_dn(:,j,ibas)=ggx2_dn(:,1)

         ggx1_dn(1,:)=phi_dn(j,:,iw)

         call dgemm('N', 'N', 1, NDN, NDN, 1.0_8, ggx1_dn, 1, g_dn(:,:,ibas), NDN, 0.0_8, ggx1_dn1, 1) !cao
         ggx1_dn=ggx1_dn1

         !ggx1_dn=matmul(ggx1_dn,g_dn(:,:,ibas))
         gr_dn(j,:,ibas)=ggx1_dn(1,:)

      120 continue

      do 140 ibas=1,NWFBAS

             gup_ii(ibas)=dot_product(gr_up(i,:,ibas),phiT_up(i,:,ibas))
             gdn_ii(ibas)=dot_product(gr_dn(j,:,ibas),phiT_dn(j,:,ibas))

      140 continue
      !!check x(i), p( x(i) )

      do ising=-1,1,2

         do ibas=1,NWFBAS

            rdet(ibas,ising,1)=1.0+ DeVpd(i,ising,1)*gup_ii(ibas)
            rdet(ibas,ising,2)=1.0+ DeVpd(i,ising,2)*gdn_ii(ibas)

            r_detbas(ibas)=rdet(ibas,ising,1)*rdet(ibas,ising,2)*detbas(ibas)
            ! cwfbas(ibas) is included(ovlps)

         end do

         ovlpNEW=sum(r_detbas)

         pt(ising)=coeffv(i,ising)*ovlpNEW/ovlpOLD/two !original probability
         if(ovlpNEW*sgn(iw) <= ZERO) then
           wgtwlkr(iw)=wgtwlkr(iw)/( 1.0-pt(ising) )
           pt(ising)=zero
         end if

      end do

      !!calculate p(1)/( p(1)+p(-1) )

      ptsum=pt(1)+pt(-1)

      if(ptsum == zero) then
        wgtwlkr(iw)=ZERO     ! walker is terminated, nothing left
        return
      end if

      ptest=pt(1)/ptsum
      if( ptest > ranGen(ISEED) ) then
        ising_v(i,j1,2)= 1
        ovlpNEW=pt(1)*ovlpOLD*two/coeffv(i,1) !see the above line
      else
        ising_v(i,j1,2)=-1
        ovlpNEW=pt(-1)*ovlpOLD*two/coeffv(i,-1)!see the above line
      end if

      itmp=ising_v(i,j1,2)           ! picked Ising_u(i)

      !!New Wave function, only i_th row are changed
      !!advance phi by exp(-deltau*V) for picked Ising variables

      tmpup=eVpd(i,itmp,1)
      tmpdn=eVpd(i,itmp,2)

      phi_up(i,:,iw)=tmpup*phi_up(i,:,iw)
      phi_dn(j,:,iw)=tmpdn*phi_dn(j,:,iw)

      !!Updating Green's(?) functions due to change of ising_u(i)
      do 200 ibas=1,NWFBAS

         tmpup=DeVpd(i,itmp,1) /rdet( ibas,itmp,1 )

         ggx2_up(:,1)=gl_up(:,i,ibas)
         ggx1_up(1,:)=gr_up(i,:,ibas)
         g_up(:,:,ibas)=g_up(:,:,ibas)-tmpup*matmul(ggx2_up,ggx1_up)

         tmpdn=DeVpd(i,itmp,2)/rdet( ibas,itmp,2 )

         ggx2_dn(:,1)=gl_dn(:,j,ibas)
         ggx1_dn(1,:)=gr_dn(j,:,ibas)
         g_dn(:,:,ibas)=g_dn(:,:,ibas)-tmpdn*matmul(ggx2_dn,ggx1_dn)

      200 continue

      !!compute determinant of the overlap integral
      do ibas=1,NWFBAS
         detbas(ibas)=rdet(ibas,itmp,1)*rdet(ibas,itmp,2)*detbas(ibas)
      end do

      ovlpOLD=ovlpNEW                ! ready for next move
      wgtwlkr(iw)=ptsum*wgtwlkr(iw)  ! renormalized walker weight

  100 continue

      ovlpDET(iw)=ovlpOLD

      return
      END subroutine V_updn

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

      SUBROUTINE V_dnup(iw)
      use cpmc
      use jiekou
      integer::i,j,k,l,m,n,j1
      integer::iw

      integer::ibas
      real(sp)::ovlpNEW,ovlpOLD
      real(sp)::tmpup,tmpdn
      !!real(sp)::s_detbas(NWFBAS),s_wgtwlkr,s_ovlpDET         ! Later
      real(sp)::rdet(NWFBAS,-1:1,NSPIN),pt(-1:1),ptsum,ptest
      real(sp)::gl_up(NUP,NSTATES,NWFBAS),gr_up(NSTATES,NUP,NWFBAS)
      real(sp)::gl_dn(NDN,NSTATES,NWFBAS),gr_dn(NSTATES,NDN,NWFBAS)
      real(sp)::gup_ii(NWFBAS),gdn_ii(NWFBAS),r_detbas(NWFBAS)
      real(sp)::ggx1_up(1,NUP),ggx2_up(NUP,1)
      real(sp)::ggx1_dn(1,NDN),ggx2_dn(NDN,1)

      real(sp)::ggx1_up1(1,NUP),ggx2_up2(NUP,1)
      real(sp)::ggx1_dn1(1,NDN),ggx2_dn2(NDN,1)

      ovlpOLD=ovlpDET(iw)

      !!Get Ising Variables, according to important sampling

   do 100 i=1,nsites

   do 100 j1=1,2

      j=idis(i,j1)

      !!For fixed i; up & down; N^2
      !!Make gl_s(k,i) = [ PhiT^+(ibas) Phi ]^-1(k,j) * PhiT^+(j,i); Ns x NSTATES
      !!gr_s(i,k) = Phi(i,j) * [ PhiT^+(ibas) Phi ]^-1(j,k);    NSTATES x Ns

      do 120 ibas=1,NWFBAS

         ggx2_up(:,1)=phiT_up(j,:,ibas)

         call dgemm('N', 'N', NUP, 1, NUP, 1.0_8, g_up(:,:,ibas), NUP, ggx2_up, NUP, 0.0_8, ggx2_up2, NUP) !cao
         ggx2_up=ggx2_up2

         !ggx2_up=matmul(g_up(:,:,ibas),ggx2_up)
         gl_up(:,j,ibas)=ggx2_up(:,1)

         ggx1_up(1,:)=phi_up(j,:,iw)

         call dgemm('N', 'N', 1, NUP, NUP, 1.0_8, ggx1_up, 1, g_up(:,:,ibas), NUP, 0.0_8, ggx1_up1, 1) !cao
         ggx1_up=ggx1_up1

         !ggx1_up=matmul(ggx1_up,g_up(:,:,ibas))
         gr_up(j,:,ibas)=ggx1_up(1,:)

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

         gup_ii(ibas)=dot_product(gr_up(j,:,ibas),phiT_up(j,:,ibas))
         gdn_ii(ibas)=dot_product(gr_dn(i,:,ibas),phiT_dn(i,:,ibas))

      140 continue

      !!check x(i), p( x(i) )

      do ising=-1,1,2

         do ibas=1,NWFBAS

           rdet(ibas,ising,1)=1.0+ DeVpd(i,ising,1)*gup_ii(ibas)
           rdet(ibas,ising,2)=1.0+ DeVpd(i,ising,2)*gdn_ii(ibas)

           r_detbas(ibas)=rdet(ibas,ising,1)*rdet(ibas,ising,2)*detbas(ibas)
           ! cwfbas(ibas) is included(ovlps)

         end do

         ovlpNEW=sum(r_detbas)

         pt(ising)=coeffv(i,ising)*ovlpNEW/ovlpOLD/two !original probability
         if(ovlpNEW*sgn(iw) <= ZERO) then
           wgtwlkr(iw)=wgtwlkr(iw)/( 1.0-pt(ising) )
           pt(ising)=zero
         end if

      end do

      !!calculate p(1)/( p(1)+p(-1) )

      ptsum=pt(1)+pt(-1)

      if(ptsum == zero) then
        wgtwlkr(iw)=ZERO     ! walker is terminated, nothing left
        return
      end if

      ptest=pt(1)/ptsum
      if( ptest > ranGen(ISEED) ) then
        ising_v(i,j1,3)= 1
        ovlpNEW=pt(1)*ovlpOLD*two/coeffv(i,1) ! see the above line
      else
        ising_v(i,j1,3)=-1
        ovlpNEW=pt(-1)*ovlpOLD*two/coeffv(i,-1)! see the above line
      end if

      itmp=ising_v(i,j1,3)           ! picked Ising_u(i)

      !!New Wave function, only i_th row are changed
      !!advance phi by exp(-deltau*V) for picked Ising variables

      tmpup=eVpd(i,itmp,1)
      tmpdn=eVpd(i,itmp,2)

      phi_up(j,:,iw)=tmpup*phi_up(j,:,iw)
      phi_dn(i,:,iw)=tmpdn*phi_dn(i,:,iw)

      !!Updating Green's(?) functions due to change of ising_u(i)
      do 200 ibas=1,NWFBAS

         tmpup=DeVpd(i,itmp,1) /rdet( ibas,itmp,1 )

         ggx2_up(:,1)=gl_up(:,j,ibas)
         ggx1_up(1,:)=gr_up(j,:,ibas)
         g_up(:,:,ibas)=g_up(:,:,ibas)-tmpup*matmul(ggx2_up,ggx1_up)

         tmpdn=DeVpd(i,itmp,2)/rdet( ibas,itmp,2 )

         ggx2_dn(:,1)=gl_dn(:,i,ibas)
         ggx1_dn(1,:)=gr_dn(i,:,ibas)
         g_dn(:,:,ibas)=g_dn(:,:,ibas)-tmpdn*matmul(ggx2_dn,ggx1_dn)

      200 continue

      !!compute determinant of the overlap integral
      do ibas=1,NWFBAS
         detbas(ibas)=rdet(ibas,itmp,1)*rdet(ibas,itmp,2)*detbas(ibas)
      end do

      ovlpOLD=ovlpNEW                ! ready for next move
      wgtwlkr(iw)=ptsum*wgtwlkr(iw)  ! renormalized walker weight

  100 continue

      ovlpDET(iw)=ovlpOLD

      return
      END subroutine V_dnup

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

      SUBROUTINE V_dndn(iw)
      use cpmc
      use jiekou
      integer::i,j,k,l,m,n,j1
      integer::iw

      integer::ibas
      real(sp)::ovlpNEW,ovlpOLD
      real(sp)::tmpdn,tmpdn1
      !!real(sp)::s_detbas(NWFBAS),s_wgtwlkr,s_ovlpDET         ! Later
      real(sp)::rdet(NWFBAS,-1:1,NSPIN),pt(-1:1),ptsum,ptest
      real(sp)::gl_dn(NDN,NSTATES,NWFBAS),gr_dn(NSTATES,NDN,NWFBAS)
      real(sp)::gl1_dn(NDN,NSTATES,NWFBAS),gr1_dn(NSTATES,NDN,NWFBAS)
      real(sp)::gdn_ii(NWFBAS),r_detbas(NWFBAS)
      real(sp)::gdn_jj(NWFBAS)
      real(sp)::ggx1_dn(1,NDN),ggx2_dn(NDN,1)
      real(sp)::gt_dn(NDN,NDN,NWFBAS,-1:1)

      real(sp)::ggx1_dn1(1,NDN),ggx2_dn2(NDN,1)

      ovlpOLD=ovlpDET(iw)

      !!Get Ising Variables, according to important sampling

   do 100 i=1,nsites

   do 100 j1=1,2

      j=idis(i,j1)

      !!For fixed i; up & down; N^2
      !!Make gl_s(k,i) = [ PhiT^+(ibas) Phi ]^-1(k,j) * PhiT^+(j,i); Ns x NSTATES
      !!gr_s(i,k) = Phi(i,j) * [ PhiT^+(ibas) Phi ]^-1(j,k);    NSTATES x Ns

      do 120 ibas=1,NWFBAS

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

         gdn_ii(ibas)=dot_product(gr_dn(i,:,ibas),phiT_dn(i,:,ibas))

      140 continue

      !!check x(i), p( x(i) )

      do ising=-1,1,2

         do ibas=1,NWFBAS

           rdet(ibas,ising,1)=1.0+ DeVpd(i,ising,1)*gdn_ii(ibas)

           tmpdn=DeVpd(i,ising,1)/rdet(ibas,ising,1)
           ggx1_dn(1,:)=gr_dn(i,:,ibas)
           ggx2_dn(:,1)=gl_dn(:,i,ibas)

           gt_dn(:,:,ibas,ising)=g_dn(:,:,ibas)-tmpdn*matmul(ggx2_dn,ggx1_dn)

           ggx2_dn(:,1)=phiT_dn(j,:,ibas)

           call dgemm('N', 'N', NDN, 1, NDN, 1.0_8, gt_dn(:,:,ibas,ising), NDN, ggx2_dn, NDN, 0.0_8, ggx2_dn2, NDN) !cao
           ggx2_dn=ggx2_dn2

           !ggx2_dn=matmul(gt_dn(:,:,ibas,ising),ggx2_dn)
           gl1_dn(:,j,ibas)=ggx2_dn(:,1)

           ggx1_dn(1,:)=phi_dn(j,:,iw)

           call dgemm('N', 'N', 1, NDN, NDN, 1.0_8, ggx1_dn, 1, gt_dn(:,:,ibas,ising), NDN, 0.0_8, ggx1_dn1, 1) !cao
           ggx1_dn=ggx1_dn1

           !ggx1_dn=matmul(ggx1_dn,gt_dn(:,:,ibas,ising))
           gr1_dn(j,:,ibas)=ggx1_dn(1,:)

           gdn_jj(ibas)=dot_product(gr1_dn(j,:,ibas),phiT_dn(j,:,ibas))

           rdet(ibas,ising,2)=1.0+DeVpd(i,ising,2)*gdn_jj(ibas)

         end do

         do ibas=1,NWFBAS
           r_detbas(ibas)=rdet(ibas,ising,1)*rdet(ibas,ising,2)*detbas(ibas)
         end do

         ovlpNEW=sum(r_detbas)

         pt(ising)=coeffv(i,ising)*ovlpNEW/ovlpOLD/two !original probability
         if(ovlpNEW*sgn(iw) <= ZERO) then
           wgtwlkr(iw)=wgtwlkr(iw)/( 1.0-pt(ising) )
           pt(ising)=zero
         end if

      end do

      !!calculate p(1)/( p(1)+p(-1) )

      ptsum=pt(1)+pt(-1)

      if(ptsum == zero) then
        wgtwlkr(iw)=ZERO     ! walker is terminated, nothing left
        return
      end if

      ptest=pt(1)/ptsum
      if( ptest > ranGen(ISEED) ) then
        ising_v(i,j1,4)= 1
        ovlpNEW=pt(1)*ovlpOLD*two/coeffv(i,1)  !see the above line
      else
        ising_v(i,j1,4)=-1
        ovlpNEW=pt(-1)*ovlpOLD*two/coeffv(i,-1) !see the above line
      end if

      itmp=ising_v(i,j1,4)           ! picked Ising_u(i)

      !!New Wave function, only i_th row are changed
      !!advance phi by exp(-deltau*V) for picked Ising variables

      tmpdn=eVpd(i,itmp,1)
      tmpdn1=eVpd(i,itmp,2)

      phi_dn(i,:,iw)=tmpdn*phi_dn(i,:,iw)
      phi_dn(j,:,iw)=tmpdn1*phi_dn(j,:,iw)

      !!Updating Green's(?) functions due to change of ising_u(i)
      do 200 ibas=1,NWFBAS

         tmpdn=DeVpd( i,itmp,2) /rdet( ibas,itmp,2 )

         ggx2_dn(:,1)=gl1_dn(:,j,ibas)
         ggx1_dn(1,:)=gr1_dn(j,:,ibas)
         g_dn(:,:,ibas)=gt_dn(:,:,ibas,itmp)-tmpdn*matmul(ggx2_dn,ggx1_dn)

      200      continue

      !!compute determinant of the overlap integral
      do ibas=1,NWFBAS
         detbas(ibas)=rdet(ibas,itmp,1)*rdet(ibas,itmp,2)*detbas(ibas)
      end do

      ovlpOLD=ovlpNEW                ! ready for next move
      wgtwlkr(iw)=ptsum*wgtwlkr(iw)  ! renormalized walker weight

  100 continue

      ovlpDET(iw)=ovlpOLD

      return
      END subroutine V_dndn


      !**************************************************!

	subroutine backphi(kp,L,phip_up,phip_dn)
	use cpmc
	use jiekou,only:stblzbk
   integer::i,j,k,L,m,n,count,istep
   integer::n1,n2,n3,itest
   integer::kp,mx,ip,pk,lsave
   real(sp)::phip_up(NUP,NSTATES),phip_dn(NDN,NSTATES)
	real(sp)::up_tmp(NUP,NSTATES,2),dn_tmp(NDN,NSTATES,2)

     !----------------------------------------------------------------!
     ! < Phi(k',L) | = < Phi_T(L) | B(N+m,kp)B(N+m-1,kp) ... B(N+1,kp)!
     !               = < Phi_T(L) | Exp(-K/2) V(N+m,kp)
     !                              Exp(-K)   V(N+m-1,kp) ...
     !                              Exp(-K)   V(N-1,kp)   Exp(-K/2)!
     !-------------------------------------------------------------!

     !-----------------------------------------------------------------!
     ! For spin up
     ! < Phi_T(L) | Exp(-K/2); (Np x Nstates) * (Nstates x Nstates) ...!
     !-----------------------------------------------------------------!

      count=0

      up_tmp(:,:,1)=phiZ_up(:,:,L)
      dn_tmp(:,:,1)=phiZ_dn(:,:,L)

      n1=1
      n2=2

      do mx=itvl_m,1,-1

	count=count+1

         if(mod(count,itvlorth)==1) then

         call dgemm('N', 'N', NUP, NSTATES, NSTATES, 1.0_8, up_tmp(:,:,n1), NUP, expT(:,:,1), NSTATES, 0.0_8, up_tmp(:,:,n2), NUP) !cao
         call dgemm('N', 'N', NDN, NSTATES, NSTATES, 1.0_8, dn_tmp(:,:,n1), NDN, expT(:,:,2), NSTATES, 0.0_8, dn_tmp(:,:,n2), NDN) !cao

          !up_tmp(:,:,n2)=matmul(up_tmp(:,:,n1),expT(:,:,1))
          !dn_tmp(:,:,n2)=matmul(dn_tmp(:,:,n1),expT(:,:,2))
         else

         call dgemm('N', 'N', NUP, NSTATES, NSTATES, 1.0_8, up_tmp(:,:,n1), NUP, exp2T(:,:,1), NSTATES, 0.0_8, up_tmp(:,:,n2), NUP) !cao
         call dgemm('N', 'N', NDN, NSTATES, NSTATES, 1.0_8, dn_tmp(:,:,n1), NDN, exp2T(:,:,2), NSTATES, 0.0_8, dn_tmp(:,:,n2), NDN) !cao

          !up_tmp(:,:,n2)=matmul(up_tmp(:,:,n1),exp2T(:,:,1))
          !dn_tmp(:,:,n2)=matmul(dn_tmp(:,:,n1),exp2T(:,:,2))
         end if

     !----------------------------------------------------------!
     ! V(N+mx,kp); (Np x Nstates) * (Nstates x Nstates) Diagonal!
     !----------------------------------------------------------!

	   do 100 ip=1,NUP
	   do 100 j=1,NSTATES

	      lsave=kexpV_s(j,mx,kp)
	      up_tmp(ip,j,n2)=up_tmp(ip,j,n2)*expV(j,lsave,1)

100        continue

            do 200 ip=1,NDN
            do 200 j=1,NSTATES

               lsave=kexpV_s(j,mx,kp)
	      dn_tmp(ip,j,n2)=dn_tmp(ip,j,n2)*expV(j,lsave,2)

200        continue

           !--------------Zhongbing-phonon fields-----------!
           if(g_ph*w_ph>0.001) then
	   do 101 ip=1,NUP
	   do 101 j=1,NSTATES

	      lsave=keph_s(j,mx,kp)
	      up_tmp(ip,j,n2)=up_tmp(ip,j,n2)*eph(j,lsave,1)

101        continue

            do 202 ip=1,NDN
            do 202 j=1,NSTATES

               lsave=keph_s(j,mx,kp)
	      dn_tmp(ip,j,n2)=dn_tmp(ip,j,n2)*eph(j,lsave,2)

202        continue
           end if
           !-------------------------------------------------!

            if(abs(vpd)>0.001) then

            do 601 i=1,lxy
            do 601 j=1,2

               lsave=keVpd_s(i,mx,kp,j,1)
               up_tmp(:,i,n2)=up_tmp(:,i,n2)*eVpd(i,lsave,1)
               up_tmp(:,idis(i,j),n2)=up_tmp(:,idis(i,j),n2)*eVpd(i,lsave,2)

601         continue

            do 602 i=1,lxy
            do 602 j=1,2

               lsave=keVpd_s(i,mx,kp,j,2)
               up_tmp(:,i,n2)=up_tmp(:,i,n2)*eVpd(i,lsave,1)
               dn_tmp(:,idis(i,j),n2)=dn_tmp(:,idis(i,j),n2)*eVpd(i,lsave,2)

602         continue

            do 603 i=1,lxy
            do 603 j=1,2

               lsave=keVpd_s(i,mx,kp,j,3)
               dn_tmp(:,i,n2)=dn_tmp(:,i,n2)*eVpd(i,lsave,2)
               up_tmp(:,idis(i,j),n2)=up_tmp(:,idis(i,j),n2)*eVpd(i,lsave,1)

603         continue
            do 604 i=1,lxy
            do 604 j=1,2

               lsave=keVpd_s(i,mx,kp,j,4)
               dn_tmp(:,i,n2)=dn_tmp(:,i,n2)*eVpd(i,lsave,1)
               dn_tmp(:,idis(i,j),n2)=dn_tmp(:,idis(i,j),n2)*eVpd(i,lsave,2)

604         continue

            end if

            n3=n2
            n2=n1
            n1=n3

            if(mod(count,itvlorth)==0) then

            call dgemm('N', 'N', NUP, NSTATES, NSTATES, 1.0_8, up_tmp(:,:,n1), NUP, expT(:,:,1), NSTATES, 0.0_8, up_tmp(:,:,n2), NUP) !cao
            call dgemm('N', 'N', NDN, NSTATES, NSTATES, 1.0_8, dn_tmp(:,:,n1), NDN, expT(:,:,2), NSTATES, 0.0_8, dn_tmp(:,:,n2), NDN) !cao

            !up_tmp(:,:,n2)=matmul(up_tmp(:,:,n1),expT(:,:,1))
            !dn_tmp(:,:,n2)=matmul(dn_tmp(:,:,n1),expT(:,:,2))   !! Exp(-K/2)!!

            n3=n2
            n2=n1
            n1=n3

            call stblzbk(kp,up_tmp(:,:,n1),dn_tmp(:,:,n1))

            end if

          end do

         phip_up=up_tmp(:,:,n1)
         phip_dn=dn_tmp(:,:,n1)

     !======================!
     ! Testing; delete later!
     !----------------------!
	itest=0
	if(itest == 1) then
	   do 350 ip=1,NUP
	   do 350 j=1,NSTATES
	      if( phip_up(ip,j) /= phiT_up(j,ip,L) ) then
		 write(*,*) phip_up(ip,j),phiT_up(j,ip,L)
		 stop
	      end if
350        continue
	   do 450 ip=1,NDN
	   do 450 j=1,NSTATES
	      if( phip_dn(ip,j) /= phiT_dn(j,ip,L) ) then
		 write(*,*) phip_dn(ip,j),phiT_dn(j,ip,L)
		 stop
	      end if
450        continue
	end if

	itest=0
	if(itest == 1) then         ! Use Trial Wave Function for Testing
	   do 300 ip=1,NUP
	   do 300 j=1,NSTATES
	      phip_up(ip,j)=phiT_up(j,ip,L)
300        continue
	   do 400 ip=1,NDN
	   do 400 j=1,NSTATES
	      phip_dn(ip,j)=phiT_dn(j,ip,L)
400        continue
	end if

	return
	end subroutine backphi


	subroutine stblzbk(kp,phi1,phi2)
	use cpmc
	use jiekou,only:modgs
	implicit none

         integer::kp
         real(sp),dimension(NUP,NSTATES)::phi1
         real(sp),dimension(NDN,NSTATES)::phi2
         real(sp)::rescale,tmpup(NSTATES,NUP),tmpdn(NSTATES,NDN)

           tmpup=transpose(phi1)
           tmpdn=transpose(phi2)

           call modgs(tmpup,NSTATES,NUP,rescale)

           call modgs(tmpdn,NSTATES,NDN,rescale)

           phi1=transpose(tmpup)
           phi2=transpose(tmpdn)

	return
	end subroutine stblzbk


     !**************************************************!

      SUBROUTINE Step1(istp)
      use cpmc
      use jiekou,only:HalfT,FullV
      integer(i4b)::i,j,k,l,m,n,iselect
      integer(i4b)::iw,istp,mx
      real(sp)::hdeltau,alpha_u

      iselect=mod(istp,itvlorth)

      do 100 iw=1,NWLKRS

       !! SKIP WALKERS WITH wgt=0
         if(wgtwlkr(iw) <= ZERO) go to 100

       !! APPLY growth control factor
         wgtwlkr(iw)=exp(deltau*etrial)*wgtwlkr(iw)

       !! APPLY EXP(-0.5*DelTau*T):

         if(iselect==1) then
         call HalfT(iw,1)
         else
         call halfT(iw,2)
         end if

         if(wgtwlkr(iw) <= ZERO) go to 100

       !! APPLY exp(-DelTau*V):
       !-------Zhongbing-phonon fields------!
         if(g_ph*w_ph>0.001) then
         call Fullph(iw)
         if(wgtwlkr(iw) <= ZERO) go to 100
         end if
       !------------------------------------!
         call FullV(iw)
         if(wgtwlkr(iw) <= ZERO) go to 100

         if(abs(vpd)>0.001) then
         call V_upup(iw)
         if(wgtwlkr(iw) <=ZERO) goto 100
         call V_updn(iw)
         if(wgtwlkr(iw) <=ZERO) goto 100
         call V_dnup(iw)
         if(wgtwlkr(iw) <=ZERO) goto 100
         call V_dndn(iw)
         if(wgtwlkr(iw) <=ZERO) goto 100
         end if

         mx=mod(istp,itvl_m)
         if(mx==0) mx=itvl_m

         if(measl==1) then

          hdeltau=half*deltau
          do i=1,nsites             ! Saved V(istp,iw)
            kexpV_s(i,mx,iw)=ising_u(i)
          end do

          if(g_ph*w_ph>0.001) then
          do i=1,nsites             ! Saved phonon fields
            keph_s(i,mx,iw)=ising_sum(i)
          end do
          end if

          if(abs(vpd)>0.001) then
          do i=1,lxy               ! Saved Vpd(istp,iw)
          do m=1,2
          do n=1,4

            keVpd_s(i,mx,iw,m,n)=ising_v(i,m,n)

          end do
          end do
          end do
          end if

         end if

       !! APPLY exp(-0.5*DelTau*T):

         if(iselect==0) call HalfT(iw,1)

100    continue

      return
      END subroutine Step1

!**************************************************!

subroutine stblz1
use cpmc
use jiekou,only:modgs
implicit none

integer::iw
real(sp)::rescale
!-----------------HoKinPara-------------------!
do iw=1,NWLKRS
!-----------------HoKinPara-------------------!
   !!skip walkers with zero weight.

   if(wgtwlkr(iw)/=0.0_sp) then

     call modgs(phi_up(:,:,iw),NSTATES,NUP,rescale)
     ovlpDET(iw)=rescale*ovlpDET(iw)

     call modgs(phi_dn(:,:,iw),NSTATES,NDN,rescale)
     ovlpDET(iw)=rescale*ovlpDET(iw)

   end if

end do

return
end subroutine stblz1