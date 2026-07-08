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


      subroutine tred2(nm,n,a,d,e,z)
      use cpmc, only: sp, use_spinor, nsites, NSO

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
      use cpmc, only: sp, use_spinor, nsites, NSO
      implicit none

      integer::I,J,K,L,M,N,II,L1,L2,NM,MML,IERR
      REAL(sp),dimension(:)::D,E
      real(sp),dimension(:,:)::Z
      REAL(sp)::B,C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2
      REAL(sp), EXTERNAL::PYTHAG

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

  subroutine dgefa(a,lda,n,ipvt,info)
   integer lda,n,ipvt(1),info
   double precision a(lda,1)
!     dgefa factors a double precision matrix by gaussian elimination.
!     dgefa is usually called by dgeco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
!     on entry
!        a       double precision(lda, n)
!                the matrix to be factored.
!        lda     integer
!                the leading dimension of the array  a .
!        n       integer
!                the order of the matrix  a .
!     on return
!        a       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!        info    integer
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that dgesl or dgedi will divide by zero
!                     if called.  use  rcond  in dgeco for a reliable
!                     indication of singularity.
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!     subroutines and functions
!     blas daxpy,dscal,idamax
!     internal variables
   double precision t
   integer j,k,kp1,l,nm1
   integer, external :: idamax      ! 函数，返回 integer
   external :: dscal, daxpy         ! 子程序   
!     gaussian elimination with partial pivoting
   info = 0
   nm1 = n - 1
   if (nm1 .lt. 1) go to 70
   do 60 k = 1, nm1
      kp1 = k + 1

!        find l = pivot index
      l = idamax(n-k+1,a(k,k),1) + k - 1
      ipvt(k) = l

!        zero pivot implies this column already triangularized
      if (a(l,k) .eq. 0.0d0) go to 40
!           interchange if necessary
         if (l .eq. k) go to 10
            t = a(l,k)
            a(l,k) = a(k,k)
            a(k,k) = t
10       continue
!           compute multipliers
         t = -1.0d0/a(k,k)
         call dscal(n-k,t,a(k+1,k),1)
!           row elimination with column indexing
         do 30 j = kp1, n
            t = a(l,j)
            if (l .eq. k) go to 20
               a(l,j) = a(k,j)
               a(k,j) = t
20          continue
            call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
30       continue
      go to 50
40    continue
         info = k
50    continue
60 continue
70 continue
   ipvt(n) = n
   if (a(n,n) .eq. 0.0d0) info = n
   return
   end

   subroutine dgedi(a,lda,n,ipvt,det,work,job)
      integer lda,n,ipvt(1),job
      double precision a(lda,1),det(2),work(1)
!     dgedi computes the determinant and inverse of a matrix
!     using the factors computed by dgeco or dgefa.
!     on entry
!        a       double precision(lda, n)
!                the output from dgeco or dgefa.
!        lda     integer
!                the leading dimension of the array  a .
!        n       integer
!                the order of the matrix  a .
!        ipvt    integer(n)
!                the pivot vector from dgeco or dgefa.
!        work    double precision(n)
!                work vector.  contents destroyed.
!        job     integer
!                = 11   both determinant and inverse.
!                = 01   inverse only.
!                = 10   determinant only.
!     on return
!        a       inverse of original matrix if requested.
!                otherwise unchanged.
!        det     double precision(2)
!                determinant of original matrix if requested.
!                otherwise not referenced.
!                determinant = det(1) * 10.0**det(2)
!                with  1.0 .le. dabs(det(1)) .lt. 10.0
!                or  det(1) .eq. 0.0 .
!     error condition
!        a division by zero will occur if the input factor contains
!        a zero on the diagonal and the inverse is requested.
!        it will not occur if the subroutines are called correctly
!        and if dgeco has set rcond .gt. 0.0 or dgefa has set
!        info .eq. 0 .
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!     subroutines and functions
!     blas 
      external :: daxpy,dscal,dswap
!     fortran dabs,mod
!     internal variables
      double precision t
      double precision ten
      integer i,j,k,kb,kp1,l,nm1
!     compute determinant
      if (job/10 .eq. 0) go to 70
         det(1) = 1.0d0
         det(2) = 0.0d0
         ten = 10.0d0
         do 50 i = 1, n
            if (ipvt(i) .ne. i) det(1) = -det(1)
            det(1) = a(i,i)*det(1)
!        ...exit
            if (det(1) .eq. 0.0d0) go to 60
   10       if (dabs(det(1)) .ge. 1.0d0) go to 20
               det(1) = ten*det(1)
               det(2) = det(2) - 1.0d0
            go to 10
   20       continue
   30       if (dabs(det(1)) .lt. ten) go to 40
               det(1) = det(1)/ten
               det(2) = det(2) + 1.0d0
            go to 30
   40       continue
   50    continue
   60    continue
   70 continue
!     compute inverse(u)
      if (mod(job,10) .eq. 0) go to 150
         do 100 k = 1, n
            a(k,k) = 1.0d0/a(k,k)
            t = -a(k,k)
            call dscal(k-1,t,a(1,k),1)
            kp1 = k + 1
            if (n .lt. kp1) go to 90
            do 80 j = kp1, n
               t = a(k,j)
               a(k,j) = 0.0d0
               call daxpy(k,t,a(1,k),1,a(1,j),1)
   80       continue
   90       continue
  100    continue
!        form inverse(u)*inverse(l)
         nm1 = n - 1
         if (nm1 .lt. 1) go to 140
         do 130 kb = 1, nm1
            k = n - kb
            kp1 = k + 1
            do 110 i = kp1, n
               work(i) = a(i,k)
               a(i,k) = 0.0d0
  110       continue
            do 120 j = kp1, n
               t = work(j)
               call daxpy(n,t,a(1,j),1,a(1,k),1)
  120       continue
            l = ipvt(k)
            if (l .ne. k) call dswap(n,a(1,k),1,a(1,l),1)
  130    continue
  140    continue
  150 continue
      return
      end   