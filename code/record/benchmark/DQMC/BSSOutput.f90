subroutine analysisOutput
    use mbss
    implicit none
    integer:: i,j,q,ix,iy
    character(50):: fileNameForOutput
    character(len=256) :: k_filename, r_filename
    character(len=32) :: kNames(20), rNames(5)
    !for analysis using Matlab
    
    !---------------variable location-------------!
    !       overall property
    !        1	sn      <n>
    !        2	sm      <m>
    !        3	sn2     <n2>
    !        4	sm2     <m2>
    !        5	s2oc    <double occ>
    !        6	s0oc    <vancacy>
    !        7	eu      <PE>
    !        8	eh      <KE>
    !        9	emu     <CE>
    !        10	elan    <LE>
    !        11	totalE
    !        12	0.0d0;
    !        13 onpair  =>on-site pair corr
    !        14 pairlen => pair correlation length
    !        15 afmssf
    !        21 spsus
    !        22 sepsus
    !        23 dpsus
    !        24 ppsus
    !        25 puppsus
    !        26 spcor
    !        27 sepcor
    !        28 dpcor
    !        29 ppcor
    !        30 puppcor
    !        31 spsussp
    !        32 sepsussp
    !        33 dpsussp
    !        34 ppsussp
    !        35 puppsussp
    !---------------site property------------------!
    !        Meas+0*Nxy	sns
    !        Meas+1*Nxy	sms
    !        Meas+2*Nxy	s2ocs
    !        Meas+3*Nxy	sds
    !        Meas+4*Nxy	sn2s
    !        Meas+5*Nxy	sm2s
    !------------------k-space propery---------------------!
    !        MEAS13+0*NKK	    csf(q)
    !        MEAS13+1*NKK	    ssfz(q)
    !        MEAS13+2*NKK	    ssfxy(q)
    !        MEAS13+3*NKK
    !        MEAS13+4*NKK	    csup(q)
    !        MEAS13+5*NKK	    ssupz(q)
    !        MEAS13+6*NKK	    ssupxy(q)
    !        MEAS13+7*NKK
    !        MEAS13+8*NKK	    ssupz(q) => single particle contribution
    !        MEAS13+9*NKK	    ssupxy(q) => single particle contribution
    !        MEAS13+10*NKK
    !        MEAS13+11*NKK        fdf(q)
    !        MEAS13+12*NKK        fdfUp(q)
    !        MEAS13+13*NKK       fdfDn(q)
    !        MEAS13+14*NKK       pmdf(q)
    !        MEAS13+15*NKK       dsf(q)
    !---------------site-site correlation------------------!
    !        Meas1+0*Nxy*Nxy	csn(i,j), <ninj> for i>j, position q=(i-1)*NXY+j ,<mimj>z for j>i, position q=(j-1)*NXY+i
    !        Meas1+1*Nxy*Nxy	sdwxy(i,j) <mimj>xy for i>j, position q=(i-1)*NXY+j
    !        Meas1+2*Nxy*Nxy	0.d0
    !        Meas1+3*Nxy*Nxy	cdw1(i,j)
    !        Meas1+4*Nxy*Nxy	sdwz1(i,j)
    !        Meas1+5*Nxy*Nxy	sdwxy1(i,j)
    !        Meas1+6*Nxy*Nxy	pmdf(i,j)
    !        Meas1+7*Nxy*Nxy	dsf(i,j)
    !        Meas1+8*Nxy*Nxy    spl
    !        Meas1+9*Nxy*Nxy    sepl
    !        Meas1+10*Nxy*Nxy    dpl
    !        Meas1+11*Nxy*Nxy    ppl
    !        Meas1+12*Nxy*Nxy    puppl
    !        Meas1+13*Nxy*Nxy    sp
    !        Meas1+14*Nxy*Nxy    sep
    !        Meas1+15*Nxy*Nxy    dp
    !        Meas1+16*Nxy*Nxy    pp
    !        Meas1+17*Nxy*Nxy    pupp
    !---------------current-current correlation--------------!
    !        Meas21+0*NT*Nxy*Nxy	curcor(i,j,tau)
    !---------------current-current correlation--------------!
    !        Meas22+0*NWM*Nxy	    curqw(i,tau)
    !---------------pair correlation------------------------!
    
    !-------------------------------------------------------!
    !        Meas4+1 nsign
    !        Meas4+2 acc
    !        Meas4+3 negs
    
    !---------------complex no------------------------------!
    !        1+0*NXYS*NXYS      gt_vt(i,j), gvt(1)=sum(1,NT) g(t,1)
    !        1+1*NXYS*NXYS      gt2(i,j,0) <c^+ c>equal time
    !        1+2*NXYS*NXYS	    gt1(i,j,0) <c c^+>equal time
    !---------------k space G-----------------------------!
    !        MEAC1+0*NK         Gkupup <c^+ c>
    !        MEAC1+1*NK         Gkdndn <c^+ c>
    !        MEAC1+2*NK         Gkupdn <c^+dn cup>
    !        MEAC1+3*NK         Gkdnup <c^+up cdn>
    !        MEAC1+4*NK         Gkupup <c c^+>
    !        MEAC1+5*NK         Gkdndn <c c^+>
    !        MEAC1+6*NK         Gkupdn <c c^+>
    !        MEAC1+7*NK         Gkdnup <c c^+>
    !---------------------Green function------------------!
    !        MEAC2+0*(NT+1)*Nxys*Nxys	gt2(i,j,time)cplc0 <c^+ c>
    !        MEAC2+1*(NT+1)*Nxys*Nxys	gt1(i,j,time)clcp0 <c c^+>
    name1(:)=' '
    name1(1)='n'
    name1(2)='m'
    name1(3)='n2'
    name1(4)='m2'
    name1(5)='doc'
    name1(6)='vac'
    name1(7)='PE'
    name1(8)='KE'
    name1(9)='CE'
    name1(10)='LE'
    name1(11)='SumE'
    name1(14)='pairlen'
    name1(15)='afmssf'
    name1(21)='spsus'
    name1(22)='sepsus'
    name1(23)='dpsus'
    name1(24)='ppsus'
    name1(25)='puppsus'
    name1(26)='spcor'
    name1(27)='sepcor'
    name1(28)='dpcor'
    name1(29)='ppcor'
    name1(30)='puppcor'
    name1(31)='spsussp'
    name1(32)='sepsussp'
    name1(33)='dpsussp'
    name1(34)='ppsussp'
    name1(35)='puppsussp'
    name1(36)='ChiralSDW'
    name1(37)='absChiralSDW'
    do i=1,NXY
        name1(MEAS+i)='nsite'
        name1(MEAS+NXY+i)='msite'
        name1(MEAS+2*NXY+i)='s2site'
        name1(MEAS+3*NXY+i)='dnsite'
        name1(MEAS+4*NXY+i)='n2site'
        name1(MEAS+5*NXY+i)='m2site'
        name1(MEAS+6*NXY+i)='ChiralSDWsite'
        name1(MEAS+7*NXY+i)='absChiralSDWsite'
    enddo
    do i=1,NK
        name1(MEAS13+0*NKK+i)='csfAA'
        name1(MEAS13+1*NKK+i)='ssfzAA'
        name1(MEAS13+2*NKK+i)='ssfxyAA'
        name1(MEAS13+4*NKK+i)='csusAA'
        name1(MEAS13+5*NKK+i)='ssuszAA'
        name1(MEAS13+6*NKK+i)='ssusxyAA'
        name1(MEAS13+8*NKK+i)='ssupzspAA'
        name1(MEAS13+9*NKK+i)='ssupxyspAA'
        name1(MEAS13+11*NKK+i)='fdfAA'
        name1(MEAS13+12*NKK+i)='fdfUpAA'
        name1(MEAS13+13*NKK+i)='fdfDnAA'
        name1(MEAS13+14*NKK+i)='pmdfAA'
        name1(MEAS13+15*NKK+i)='dsfAA'
    enddo
    
    name1(MEAS1+1:MEAS1+NXY*NXY)='csn'
    name1(MEAS1+NXY*NXY+1:MEAS1+2*NXY*NXY)='sdwxy'
    name1(MEAS1+3*NXY*NXY+1:MEAS1+4*NXY*NXY)='cdw1'
    name1(MEAS1+4*NXY*NXY+1:MEAS1+5*NXY*NXY)='sdwz1'
    name1(MEAS1+5*NXY*NXY+1:MEAS1+6*NXY*NXY)='sdwxy1'
    name1(MEAS1+6*NXY*NXY+1:MEAS1+7*NXY*NXY)='pmdf'
    name1(MEAS1+7*NXY*NXY+1:MEAS1+8*NXY*NXY)='dsf'
    name1(MEAS1+8*NXY*NXY+1:MEAS1+9*NXY*NXY)='spl'
    name1(MEAS1+9*NXY*NXY+1:MEAS1+10*NXY*NXY)='sepl'
    name1(MEAS1+10*NXY*NXY+1:MEAS1+11*NXY*NXY)='dpl'
    name1(MEAS1+11*NXY*NXY+1:MEAS1+12*NXY*NXY)='ppl'
    name1(MEAS1+12*NXY*NXY+1:MEAS1+13*NXY*NXY)='puppl'
    name1(MEAS1+13*NXY*NXY+1:MEAS1+14*NXY*NXY)='sp'
    name1(MEAS1+14*NXY*NXY+1:MEAS1+15*NXY*NXY)='sep'
    name1(MEAS1+15*NXY*NXY+1:MEAS1+16*NXY*NXY)='dp'
    name1(MEAS1+16*NXY*NXY+1:MEAS1+17*NXY*NXY)='pp'
    name1(MEAS1+17*NXY*NXY+1:MEAS1+18*NXY*NXY)='pupp'

    ! Vertex
    name1(MEAS_VTX_R + 0*NXY*NXY + 1 : MEAS_VTX_R + 1*NXY*NXY) = 'vertex_s'
    name1(MEAS_VTX_R + 1*NXY*NXY + 1 : MEAS_VTX_R + 2*NXY*NXY) = 'vertex_se'
    name1(MEAS_VTX_R + 2*NXY*NXY + 1 : MEAS_VTX_R + 3*NXY*NXY) = 'vertex_d'
    name1(MEAS_VTX_R + 3*NXY*NXY + 1 : MEAS_VTX_R + 4*NXY*NXY) = 'vertex_p'
    name1(MEAS_VTX_R + 4*NXY*NXY + 1 : MEAS_VTX_R + 5*NXY*NXY) = 'vertex_pup'

    ! (k-space) Vertex
    do i=1, NK
        name1(MEAS_VTX_K + 0*NK + i) = 'vertex_s_k'
        name1(MEAS_VTX_K + 1*NK + i) = 'vertex_se_k'
        name1(MEAS_VTX_K + 2*NK + i) = 'vertex_d_k'
        name1(MEAS_VTX_K + 3*NK + i) = 'vertex_p_k'
        name1(MEAS_VTX_K + 4*NK + i) = 'vertex_pup_k'
    enddo

    name1(Meas21+1:Meas21+NT*Nxy*Nxy)='curcor'
    name1(Meas22+1:Meas22+NWM*Nxy)='curqw'
    name1(Meas4+1)='nsign'
    name1(Meas4+2)='acc'
    name1(Meas4+3)='negs'
    
    name2(1:NXYS*NXYS)='gtvt'
    name2(1*NXYS*NXYS+1:2*NXYS*NXYS)='g2'
    name2(2*NXYS*NXYS+1:3*NXYS*NXYS)='g1'
    name2(MEAC1+0*NK+1:MEAC1+1*NK)='g2kuu'
    name2(MEAC1+1*NK+1:MEAC1+2*NK)='g2kdd'
    name2(MEAC1+2*NK+1:MEAC1+3*NK)='g2kud'
    name2(MEAC1+3*NK+1:MEAC1+4*NK)='g2kdu'
    name2(MEAC1+4*NK+1:MEAC1+5*NK)='g1kuu'
    name2(MEAC1+5*NK+1:MEAC1+6*NK)='g1kdd'
    name2(MEAC1+6*NK+1:MEAC1+7*NK)='g1kud'
    name2(MEAC1+7*NK+1:MEAC1+8*NK)='g1kdu'
    name2(MEAC2+0*(NT+1)*Nxys*Nxys+1:MEAC2+1*(NT+1)*Nxys*Nxys)='cplc0'
    name2(MEAC2+1*(NT+1)*Nxys*Nxys+1:MEAC2+2*(NT+1)*Nxys*Nxys)='clcp0'
    name3(1)='NT'
    name3(2)='dt'
    name3(3)='beta'
    name3(4)='U'
    name3(5)='mu'
    name3(6)='warms'
    name3(7)='runs'
    name3(8)='sweeps'
    name3(9)='landa'
    name3(10)='NX'
    name3(11)='NY'
    name3(12)='NLA'
    name3(13)='h'
    name3(14)='tam'
    name3(15)='NWM'
    do i=1,MEAS21
        write(800+ncpu,'(1 a,2 f25.15)') name1(i),ram(i,1),ram(i,2)
    enddo
    do i=Meas22+1,Meas2
        write(800+ncpu,'(1 a,2 f25.15)') name1(i),ram(i,1),ram(i,2)
    enddo
    !do i=1,MEAS2
    !    write(800+ncpu,'(1 a,2 f20.15)') name1(i),ram(i,1),ram(i,2)
    !enddo
    !do i=MEAS3,MEAS4
    !    write(800+ncpu,'(1 a,2 f20.15)') name1(i),ram(i,1),ram(i,2)
    !enddo
    !do i=1,MEAC2
    !    write(800+ncpu,'(1 a,4 f20.15)') name2(i),real(ramc(i,1)),aimag(ramc(i,1)),real(ramc(i,2)),aimag(ramc(i,2))
    !enddo
    !do i=1,MEAC3
    !    write(800+ncpu,'(1 a,4 f20.15)') name2(i),real(ramc(i,1)),aimag(ramc(i,1)),real(ramc(i,2)),aimag(ramc(i,2))
    !enddo
    write(800+ncpu,'(1 a,2 f25.15)') name1(MEAS4+1),rsign(1),rsign(2)
    do i=MEAS4+2,MEAS5
        write(800+ncpu,'(1 a,f25.15)') name1(i),mah(i)
    enddo
    
    write(800+ncpu,'(1 a,i5)') name3(1),NT
    write(800+ncpu,'(1 a,f20.15)') name3(2),dt
    write(800+ncpu,'(1 a,f20.15)') name3(3),beta
    write(800+ncpu,'(1 a,f20.15)') name3(4),ue
    write(800+ncpu,'(1 a,f20.15)') name3(5),mu
    write(800+ncpu,'(1 a,i5)') name3(6),warms
    write(800+ncpu,'(1 a,i5)') name3(7),runs
    write(800+ncpu,'(1 a,i5)') name3(8),sweeps
    write(800+ncpu,'(1 a,2 f20.15)') name3(9),lamda
    write(800+ncpu,'(1 a,i5)') name3(10),NX
    write(800+ncpu,'(1 a,i5)') name3(11),NY
    write(800+ncpu,'(1 a,i5)') name3(12),NLA
    write(800+ncpu,'(1 a,f20.15)') name3(13),h
    write(800+ncpu,'(1 a,f20.15)') name3(14),tam
    write(800+ncpu,'(1 a,i5)') name3(15),NWM
    
    write(350,NXYSstr) dreal(ramc(MEAC2+1:MEAC2+NXYS*NXYS,1))
    
    ! k-space output
    kNames = [character(len=32) :: 'spsus', 'sepsus', 'dpsus', 'ppsus', 'puppsus', &
    'spcor', 'sepcor', 'dpcor', 'ppcor', 'puppcor', &
    'spsussp', 'sepsussp', 'dpsussp', 'ppsussp', 'puppsussp', &
    'vertex_s', 'vertex_se', 'vertex_d', 'vertex_p', 'vertex_pup']

    do j = 1, 15
        write(k_filename, '("PairKSpace/", A, "_CPU", I2.2, ".dat")') trim(kNames(j)), ncpu
        open(unit=200+j, file=trim(k_filename), status='unknown')

        write(200+j, '(A15, A15, A20, A20)') 'kx', 'ky', 'value', 'error'
        do q = 1, NK
            write(200+j, '(2F15.9, 2E20.10)') kx(q), ky(q), ram(MEAS14+(j-1)*NKK+q, 1), ram(MEAS14+(j-1)*NKK+q, 2)
        enddo
        close(200+j)
    enddo

    do j = 1, 5
        write(k_filename, '("PairKSpace/", A, "_CPU", I2.2, ".dat")') trim(kNames(15+j)), ncpu
        open(unit=215+j, file=trim(k_filename), status='unknown')

        write(215+j, '(A15, A15, A20, A20)') 'kx', 'ky', 'value', 'error'
        do q = 1, NK
            write(215+j, '(2F15.9, 2E20.10)') kx(q), ky(q), ram(MEAS_VTX_K+(j-1)*NK+q, 1), ram(MEAS_VTX_K+(j-1)*NK+q, 2)
        enddo
        close(215+j)
    enddo

    rNames = [character(len=32) :: 'vertex_s', 'vertex_se', 'vertex_d', 'vertex_p', 'vertex_pup']

    do j = 1, 5
        write(r_filename, '("PairRSpace/", A, "_CPU", I2.2, ".dat")') trim(rNames(j)), ncpu
        open(unit=300+j, file=trim(r_filename), status='unknown')

        write(300+j, '(A8, A8, A20, A20)') 'ix', 'iy', 'value', 'error'
        do ix = 1, NXY
            do iy = 1, NXY
                q = (ix-1)*NXY + iy
                write(300+j, '(2I8, 2E20.10)') ix, iy, ram(MEAS_VTX_R+(j-1)*NXY*NXY+q, 1), ram(MEAS_VTX_R+(j-1)*NXY*NXY+q, 2)
            enddo
        enddo
        close(300+j)
    enddo

    endsubroutine analysisOutput
    
    !if having parallel computing cpuLabel>0
    subroutine binOutput(cpuLabel,intBin)
    use mbss
    implicit none
    integer:: i,j,q
    integer:: intBin
    integer:: cpuLabel
    character(50):: fileNameForOutput
    if (cpuLabel==0) then
        write (fileNameForOutput,'(I2.2,A)') intBin,'.txt'
        fileNameForOutput=trim(fileNameForOutput)
    else
        write (fileNameForOutput,'(I2.2,A,I2.2,A)') intBin,'CPU',cpuLabel,'.txt'
        fileNameForOutput=trim(fileNameForOutput)
    endif
    
    open(unit=8,file='DataBin'//fileNameForOutput)
    
    !---------------variable location-------------!
    !       overall property
    !        1	sn      <n>
    !        2	sm      <m>
    !        3	sn2     <n2>
    !        4	sm2     <m2>
    !        5	s2oc    <double occ>
    !        6	s0oc    <vancacy>
    !        7	eu      <PE>
    !        8	eh      <KE>
    !        9	emu     <CE>
    !        10	elan    <LE>
    !        11	totalE
    !        12	0.0d0;
    !        13 onpair  =>on-site pair corr
    !        14 pairlen => pair correlation length
    !        15 afmssf
    !        21 spsus
    !        22 sepsus
    !        23 dpsus
    !        24 ppsus
    !        25 puppsus
    !        26 spcor
    !        27 sepcor
    !        28 dpcor
    !        29 ppcor
    !        30 puppcor
    !        31 spsussp
    !        32 sepsussp
    !        33 dpsussp
    !        34 ppsussp
    !        35 puppsussp
    !---------------site property------------------!
    !        Meas+0*Nxy	sns
    !        Meas+1*Nxy	sms
    !        Meas+2*Nxy	s2ocs
    !        Meas+3*Nxy	sds
    !        Meas+4*Nxy	sn2s
    !        Meas+5*Nxy	sm2s
    !------------------k-space propery---------------------!
    !        MEAS13+0*NKK	    csf(q)
    !        MEAS13+1*NKK	    ssfz(q)
    !        MEAS13+2*NKK	    ssfxy(q)
    !        MEAS13+3*NKK
    !        MEAS13+4*NKK	    csup(q)
    !        MEAS13+5*NKK	    ssupz(q)
    !        MEAS13+6*NKK	    ssupxy(q)
    !        MEAS13+7*NKK
    !        MEAS13+8*NKK	    ssupz(q) => single particle contribution
    !        MEAS13+9*NKK	    ssupxy(q) => single particle contribution
    !        MEAS13+10*NKK
    !        MEAS13+11*NKK        fdf(q)
    !        MEAS13+12*NKK        fdfUp(q)
    !        MEAS13+13*NKK       fdfDn(q)
    !        MEAS13+14*NKK       pmdf(q)
    !        MEAS13+15*NKK       dsf(q)
    !---------------site-site correlation------------------!
    !        Meas1+0*Nxy*Nxy	csn(i,j), <ninj> for i>j, position q=(i-1)*NXY+j ,<mimj>z for j>i, position q=(j-1)*NXY+i
    !        Meas1+1*Nxy*Nxy	sdwxy(i,j) <mimj>xy for i>j, position q=(i-1)*NXY+j
    !        Meas1+2*Nxy*Nxy	0.d0
    !        Meas1+3*Nxy*Nxy	cdw1(i,j)
    !        Meas1+4*Nxy*Nxy	sdwz1(i,j)
    !        Meas1+5*Nxy*Nxy	sdwxy1(i,j)
    !        Meas1+6*Nxy*Nxy	pmdf(i,j)
    !        Meas1+7*Nxy*Nxy	dsf(i,j)
    !        Meas1+8*Nxy*Nxy    spl
    !        Meas1+9*Nxy*Nxy    sepl
    !        Meas1+10*Nxy*Nxy    dpl
    !        Meas1+11*Nxy*Nxy    ppl
    !        Meas1+12*Nxy*Nxy    puppl
    !        Meas1+13*Nxy*Nxy    sp
    !        Meas1+14*Nxy*Nxy    sep
    !        Meas1+15*Nxy*Nxy    dp
    !        Meas1+16*Nxy*Nxy    pp
    !        Meas1+17*Nxy*Nxy    pupp
    !---------------current-current correlation--------------!
    !        Meas21+0*NT*Nxy*Nxy	curcor(i,j,tau)
    !---------------current-current correlation--------------!
    !        Meas22+0*NWM*Nxy	    curqw(i,tau)
    !---------------pair correlation------------------------!
    
    !-------------------------------------------------------!
    !        Meas4+1 nsign
    !        Meas4+2 acc
    !        Meas4+3 negs
    
    !---------------complex no------------------------------!
    !        1+0*NXYS*NXYS      gt_vt(i,j), gvt(1)=sum(1,NT) g(t,1)
    !        1+1*NXYS*NXYS      gt2(i,j,0) <c^+ c>equal time
    !        1+2*NXYS*NXYS	    gt1(i,j,0) <c c^+>equal time
    !---------------k space G-----------------------------!
    !        MEAC1+0*NK         Gkupup <c^+ c>
    !        MEAC1+1*NK         Gkdndn <c^+ c>
    !        MEAC1+2*NK         Gkupdn <c^+dn cup>
    !        MEAC1+3*NK         Gkdnup <c^+up cdn>
    !        MEAC1+4*NK         Gkupup <c c^+>
    !        MEAC1+5*NK         Gkdndn <c c^+>
    !        MEAC1+6*NK         Gkupdn <c c^+>
    !        MEAC1+7*NK         Gkdnup <c c^+>
    !---------------------Green function------------------!
    !        MEAC2+0*(NT+1)*Nxys*Nxys	gt2(i,j,time)cplc0 <c^+ c>
    !        MEAC2+1*(NT+1)*Nxys*Nxys	gt1(i,j,time)clcp0 <c c^+>
    name1(:)=' '
    name1(1)='n'
    name1(2)='m'
    name1(3)='n2'
    name1(4)='m2'
    name1(5)='doc'
    name1(6)='vac'
    name1(7)='PE'
    name1(8)='KE'
    name1(9)='CE'
    name1(10)='LE'
    name1(11)='SumE'
    name1(13)='onpair'
    name1(14)='pairlen'
    name1(15)='afmssf'
    name1(21)='spsus'
    name1(22)='sepsus'
    name1(23)='dpsus'
    name1(24)='ppsus'
    name1(25)='puppsus'
    name1(26)='spcor'
    name1(27)='sepcor'
    name1(28)='dpcor'
    name1(29)='ppcor'
    name1(30)='puppcor'
    name1(31)='spsussp'
    name1(32)='sepsussp'
    name1(33)='dpsussp'
    name1(34)='ppsussp'
    name1(35)='puppsussp'
    name1(36)='ChiralSDW'
    name1(37)='absChiralSDW'
    do i=1,NXY
        name1(MEAS+i)='nsite'
        name1(MEAS+NXY+i)='msite'
        name1(MEAS+2*NXY+i)='s2site'
        name1(MEAS+3*NXY+i)='dnsite'
        name1(MEAS+4*NXY+i)='n2site'
        name1(MEAS+5*NXY+i)='m2site'
        name1(MEAS+6*NXY+i)='ChiralSDWsite'
        name1(MEAS+7*NXY+i)='absChiralSDWsite'
    enddo
    do i=1,NK
        name1(MEAS13+0*NKK+i)='csfAA'
        name1(MEAS13+1*NKK+i)='ssfzAA'
        name1(MEAS13+2*NKK+i)='ssfxyAA'
        name1(MEAS13+4*NKK+i)='csusAA'
        name1(MEAS13+5*NKK+i)='ssuszAA'
        name1(MEAS13+6*NKK+i)='ssusxyAA'
        name1(MEAS13+8*NKK+i)='ssupzspAA'
        name1(MEAS13+9*NKK+i)='ssupxyspAA'
        name1(MEAS13+11*NKK+i)='fdfAA'
        name1(MEAS13+12*NKK+i)='fdfUpAA'
        name1(MEAS13+13*NKK+i)='fdfDnAA'
        name1(MEAS13+14*NKK+i)='pmdfAA'
        name1(MEAS13+15*NKK+i)='dsfAA'
    enddo
    
    name1(MEAS1+1:MEAS1+NXY*NXY)='csn'
    name1(MEAS1+NXY*NXY+1:MEAS1+2*NXY*NXY)='sdwxy'
    name1(MEAS1+3*NXY*NXY+1:MEAS1+4*NXY*NXY)='cdw1'
    name1(MEAS1+4*NXY*NXY+1:MEAS1+5*NXY*NXY)='sdwz1'
    name1(MEAS1+5*NXY*NXY+1:MEAS1+6*NXY*NXY)='sdwxy1'
    name1(MEAS1+6*NXY*NXY+1:MEAS1+7*NXY*NXY)='pmdf'
    name1(MEAS1+7*NXY*NXY+1:MEAS1+8*NXY*NXY)='dsf'
    name1(MEAS1+8*NXY*NXY+1:MEAS1+9*NXY*NXY)='spl'
    name1(MEAS1+9*NXY*NXY+1:MEAS1+10*NXY*NXY)='sepl'
    name1(MEAS1+10*NXY*NXY+1:MEAS1+11*NXY*NXY)='dpl'
    name1(MEAS1+11*NXY*NXY+1:MEAS1+12*NXY*NXY)='ppl'
    name1(MEAS1+12*NXY*NXY+1:MEAS1+13*NXY*NXY)='puppl'
    name1(MEAS1+13*NXY*NXY+1:MEAS1+14*NXY*NXY)='sp'
    name1(MEAS1+14*NXY*NXY+1:MEAS1+15*NXY*NXY)='sep'
    name1(MEAS1+15*NXY*NXY+1:MEAS1+16*NXY*NXY)='dp'
    name1(MEAS1+16*NXY*NXY+1:MEAS1+17*NXY*NXY)='pp'
    name1(MEAS1+17*NXY*NXY+1:MEAS1+18*NXY*NXY)='pupp'

    ! Vertex
    name1(MEAS_VTX_R + 0*NXY*NXY + 1 : MEAS_VTX_R + 1*NXY*NXY) = 'vertex_s'
    name1(MEAS_VTX_R + 1*NXY*NXY + 1 : MEAS_VTX_R + 2*NXY*NXY) = 'vertex_se'
    name1(MEAS_VTX_R + 2*NXY*NXY + 1 : MEAS_VTX_R + 3*NXY*NXY) = 'vertex_d'
    name1(MEAS_VTX_R + 3*NXY*NXY + 1 : MEAS_VTX_R + 4*NXY*NXY) = 'vertex_p'
    name1(MEAS_VTX_R + 4*NXY*NXY + 1 : MEAS_VTX_R + 5*NXY*NXY) = 'vertex_pup'

    ! (k-space) Vertex
    do i=1, NK
        name1(MEAS_VTX_K + 0*NK + i) = 'vertex_s_k'
        name1(MEAS_VTX_K + 1*NK + i) = 'vertex_se_k'
        name1(MEAS_VTX_K + 2*NK + i) = 'vertex_d_k'
        name1(MEAS_VTX_K + 3*NK + i) = 'vertex_p_k'
        name1(MEAS_VTX_K + 4*NK + i) = 'vertex_pup_k'
    enddo

    name1(Meas21+1:Meas21+NT*Nxy*Nxy)='curcor'
    name1(Meas22+1:Meas22+NWM*Nxy)='curqw'
    name1(Meas4+1)='nsign'
    name1(Meas4+2)='acc'
    name1(Meas4+3)='negs'
    
    name2(1:NXYS*NXYS)='gtvt'
    name2(1*NXYS*NXYS+1:2*NXYS*NXYS)='g2'
    name2(2*NXYS*NXYS+1:3*NXYS*NXYS)='g1'
    name2(MEAC1+0*NK+1:MEAC1+1*NK)='g2kuu'
    name2(MEAC1+1*NK+1:MEAC1+2*NK)='g2kdd'
    name2(MEAC1+2*NK+1:MEAC1+3*NK)='g2kud'
    name2(MEAC1+3*NK+1:MEAC1+4*NK)='g2kdu'
    name2(MEAC1+4*NK+1:MEAC1+5*NK)='g1kuu'
    name2(MEAC1+5*NK+1:MEAC1+6*NK)='g1kdd'
    name2(MEAC1+6*NK+1:MEAC1+7*NK)='g1kud'
    name2(MEAC1+7*NK+1:MEAC1+8*NK)='g1kdu'
    name2(MEAC2+0*(NT+1)*Nxys*Nxys+1:MEAC2+1*(NT+1)*Nxys*Nxys)='cplc0'
    name2(MEAC2+1*(NT+1)*Nxys*Nxys+1:MEAC2+2*(NT+1)*Nxys*Nxys)='clcp0'
    name3(1)='NT'
    name3(2)='dt'
    name3(3)='beta'
    name3(4)='U'
    name3(5)='mu'
    name3(6)='warms'
    name3(7)='runs'
    name3(8)='sweeps'
    name3(9)='lamda'
    name3(10)='NX'
    name3(11)='NY'
    name3(12)='NLA'
    name3(13)='h'
    name3(14)='tam'
    name3(15)='NWM'
    do i=1,MEAS1
        write(8,'(1 a,2 f25.15)') name1(i),mah(i)
    enddo
    do i=Meas22+1,Meas2
        write(8,'(1 a,2 f25.15)') name1(i),mah(i)
    enddo
    !do i=1,MEAS2
    !    write(8,'(1 a,2 f20.15)') name1(i),ram(i,1),ram(i,2)
    !enddo
    !do i=MEAS3,MEAS4
    !    write(8,'(1 a,2 f20.15)') name1(i),ram(i,1),ram(i,2)
    !enddo
    !do i=1,MEAC2
    !    write(8,'(1 a,4 f20.15)') name2(i),real(ramc(i,1)),aimag(ramc(i,1)),real(ramc(i,2)),aimag(ramc(i,2))
    !enddo
    !do i=1,MEAC3
    !    write(8,'(1 a,4 f20.15)') name2(i),real(ramc(i,1)),aimag(ramc(i,1)),real(ramc(i,2)),aimag(ramc(i,2))
    !enddo
    !write(8,'(1 a,2 f20.15)') name1(MEAS4+1),rsign(1),rsign(2)
    do i=MEAS4+1,MEAS5
        write(8,'(1 a,f25.15)') name1(i),mah(i)
    enddo
    
    write(8,'(1 a,i5)') name3(1),NT
    write(8,'(1 a,f20.15)') name3(2),dt
    write(8,'(1 a,f20.15)') name3(3),beta
    write(8,'(1 a,f20.15)') name3(4),ue
    write(8,'(1 a,f20.15)') name3(5),mu
    write(8,'(1 a,i5)') name3(6),warms
    write(8,'(1 a,i5)') name3(7),runs
    write(8,'(1 a,i5)') name3(8),sweeps
    write(8,'(1 a,2 f20.15)') name3(9),lamda
    write(8,'(1 a,i5)') name3(10),NX
    write(8,'(1 a,i5)') name3(11),NY
    write(8,'(1 a,i5)') name3(12),NLA
    write(8,'(1 a,f20.15)') name3(13),h
    write(8,'(1 a,f20.15)') name3(14),tam
    write(8,'(1 a,i5)') name3(15),NWM

    write(8, '("=== K-Space Pairing Data ===")')
    do i = 1, 10*NKK
        write(8, '(f25.15)') mah(MEAS14 + i)
    enddo

    write(8, '("=== K-Space Vertex Data ===")')
    do i = 1, 5*NK
        write(8, '(f25.15)') mah(MEAS_VTX_K + i)
    enddo

    write(8, '("=== R-Space Vertex Data ===")')
    do i = 1, 5*NXY*NXY
        write(8, '(f25.15)') mah(MEAS_VTX_R + i)
    enddo
    close(8)
    endsubroutine binOutput
    
    