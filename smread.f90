! These subroutines read/bin/interpolate soil moisture data for CCAM.
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine reads soil moisture data (i.e., assume resolution
! is >=10km so that we simply stream the data).  Replace with ecoveg
! style aggregate-on-read if there is a performance or memory issue.

Subroutine getsmdata(ncid,dataout,lsmask,grid,tlld,ccdim,wlev,lonlat,xyz,axyz,bxyz)

Use ccinterp

Implicit None

Include "netcdf.inc"

integer, intent(in) :: wlev,ncid
Integer, dimension(1:2), intent(in) :: ccdim
Integer, dimension(1:ccdim(1),1:ccdim(2)) :: countn,countm,counto
Integer, dimension(1:4) :: ncsize
Integer, dimension(1:4,1:2) :: arrsize
Integer, dimension(1:2) :: pxy
Integer ilat,ilon,lci,lcj,i,j,nface,k
Integer ncstatus,varid
Character*80, dimension(1:2) :: varname
Character*9 formout
character*2 chr
real, dimension(2), intent(in) :: lonlat
real, dimension(ccdim(1),ccdim(2),3), intent(in) :: xyz,axyz,bxyz
Real, dimension(1:ccdim(1),1:ccdim(2)), intent(in) :: lsmask
Real, dimension(1:ccdim(1),1:ccdim(2),1:63+4*wlev), intent(out) :: dataout
Real, dimension(1:ccdim(1),1:ccdim(2)), intent(in) :: grid
Real, dimension(1:ccdim(1),1:ccdim(2),1:2), intent(in) :: tlld
Real, dimension(1:ccdim(1),1:ccdim(2),1:2) :: rlld
Real, dimension(1:2,1:3) :: smlonlat
Real, dimension(:,:,:), allocatable :: coverout
Real aglon,aglat,alci,alcj,serlon,serlat
Real minlon,maxlon
real zonx,zony,zonz,polenz,poleny,polenx,costh,sinth,den
real, dimension(ccdim(1)*ccdim(2)) :: x,y,z
real, dimension(1:wlev) :: uzon,vmer
Logical, dimension(:,:), allocatable :: sermask

dataout=0. ! store data
countn=0   ! count for land points
countm=0   ! count for seaice points
counto=0.  ! count for water points

! set-up rotations for currents
polenx=-cos(lonlat(2)*3.1415926536/180.)
poleny=0.
polenz=sin(lonlat(1)*3.1415926536/180.)

! Get sm dimensions
Call getncdims(ncid,ncsize)
Call getnclonlat(ncid,smlonlat(1:2,1:2))
smlonlat(1,3)=Real(ncsize(1))
smlonlat(2,3)=Real(ncsize(2))

minlon=Min(smlonlat(1,1),smlonlat(1,2))
maxlon=minlon+360.
rlld=tlld
Do while (Any(rlld(:,:,1).LT.minlon))
  Where (rlld(:,:,1).LT.minlon)
    rlld(:,:,1)=rlld(:,:,1)+360.
  End where
End do
Do while (Any(rlld(:,:,1).GT.maxlon))
  Where (rlld(:,:,1).GT.maxlon)
    rlld(:,:,1)=rlld(:,:,1)-360.
  End where
End do

! Get size of slab
arrsize=1
arrsize(1:2,2)=ncsize(1:2)
Allocate(coverout(1:arrsize(1,2),1:arrsize(2,2),1:63+4*wlev))
coverout=0.

! Read data
If (nf_inq_varid(ncid,'wetfrac1',varid).EQ.nf_noerr) Then
  Write(6,*) "Reading wetfrac1/2/3/4/5/6"
  varname(1)='wetfrac1'
  varname(2)='none'
  Call getmeta(ncid,varname,coverout(:,:,1),arrsize)
  varname(1)='wetfrac2'
  varname(2)='none'
  Call getmeta(ncid,varname,coverout(:,:,2),arrsize)
  varname(1)='wetfrac3'
  varname(2)='none'
  Call getmeta(ncid,varname,coverout(:,:,3),arrsize)
  varname(1)='wetfrac4'
  varname(2)='none'
  Call getmeta(ncid,varname,coverout(:,:,4),arrsize)
  varname(1)='wetfrac5'
  varname(2)='none'
  Call getmeta(ncid,varname,coverout(:,:,5),arrsize)
  varname(1)='wetfrac6'
  varname(2)='none'
  Call getmeta(ncid,varname,coverout(:,:,6),arrsize)
Else If ((nf_inq_varid(ncid,'wetfrac3',varid).EQ.nf_noerr).AND.(nf_inq_varid(ncid,'wetfrac5',varid).EQ.nf_noerr)) Then
  Write(6,*) "Reading wetfrac3/5"
  varname(1)='wetfrac3'
  varname(2)='none'
  Call getmeta(ncid,varname,coverout(:,:,1),arrsize)
  coverout(:,:,2)=coverout(:,:,1)
  coverout(:,:,3)=coverout(:,:,1)
  varname(1)='wetfrac5'
  varname(2)='none'
  Call getmeta(ncid,varname,coverout(:,:,4),arrsize)
  coverout(:,:,5)=coverout(:,:,4)
  coverout(:,:,6)=coverout(:,:,4)
Else If (nf_inq_varid(ncid,'wetfrac',varid).EQ.nf_noerr) Then
  Write(6,*) "Reading wetfrac"
  varname(1)='wetfrac'
  varname(2)='none'
  Call getmeta(ncid,varname,coverout(:,:,1),arrsize)
  coverout(:,:,2)=coverout(:,:,1)
  coverout(:,:,3)=coverout(:,:,1)
  coverout(:,:,4)=coverout(:,:,1)
  coverout(:,:,5)=coverout(:,:,1)
  coverout(:,:,6)=coverout(:,:,1)
Else
  Write(6,*) "ERROR: Cannot find valid soil mositure data in input file"
  Stop
End If

If (nf_inq_varid(ncid,'tgg1',varid).EQ.nf_noerr) Then
  Write(6,*) "Reading tgg1/2/3/4/5/6"
  varname(1)='tgg1'
  varname(2)='K'
  Call getmeta(ncid,varname,coverout(:,:,7),arrsize)
  varname(1)='tgg2'
  varname(2)='K'
  Call getmeta(ncid,varname,coverout(:,:,8),arrsize)
  varname(1)='tgg3'
  varname(2)='K'
  Call getmeta(ncid,varname,coverout(:,:,9),arrsize)
  varname(1)='tgg4'
  varname(2)='K'
  Call getmeta(ncid,varname,coverout(:,:,10),arrsize)
  varname(1)='tgg5'
  varname(2)='K'
  Call getmeta(ncid,varname,coverout(:,:,11),arrsize)
  varname(1)='tgg6'
  varname(2)='K'
  Call getmeta(ncid,varname,coverout(:,:,12),arrsize)
Else If ((nf_inq_varid(ncid,'tgg3',varid).EQ.nf_noerr).AND.(nf_inq_varid(ncid,'tgg5',varid).EQ.nf_noerr)) Then
  Write(6,*) "Reading tgg3/5"
  varname(1)='tgg3'
  varname(2)='K'
  Call getmeta(ncid,varname,coverout(:,:,7),arrsize)
  coverout(:,:,8)=coverout(:,:,7)
  coverout(:,:,9)=coverout(:,:,7)
  varname(1)='tgg5'
  varname(2)='K'
  Call getmeta(ncid,varname,coverout(:,:,10),arrsize)
  coverout(:,:,11)=coverout(:,:,10)
  coverout(:,:,12)=coverout(:,:,10)
Else
  Write(6,*) "ERROR: Cannot find valid soil temperature data in input file"
  Stop
End If

If (nf_inq_varid(ncid,'rooftgg4',varid).EQ.nf_noerr) Then
  Write(6,*) "Reading rooftgg1/2/3/4"
  varname(1)='rooftgg1'
  varname(2)='K'
  Call getmeta(ncid,varname,coverout(:,:,13),arrsize)
  varname(1)='rooftgg2'
  varname(2)='K'
  Call getmeta(ncid,varname,coverout(:,:,14),arrsize)
  varname(1)='rooftgg3'
  varname(2)='K'
  Call getmeta(ncid,varname,coverout(:,:,15),arrsize)
  varname(1)='rooftgg4'
  varname(2)='K'
  Call getmeta(ncid,varname,coverout(:,:,16),arrsize)
  
  Write(6,*) "Reading waletgg1/2/3/4"
  varname(1)='waletgg1'
  varname(2)='K'
  Call getmeta(ncid,varname,coverout(:,:,17),arrsize)
  varname(1)='waletgg2'
  varname(2)='K'
  Call getmeta(ncid,varname,coverout(:,:,18),arrsize)
  varname(1)='waletgg3'
  varname(2)='K'
  Call getmeta(ncid,varname,coverout(:,:,19),arrsize)
  varname(1)='waletgg4'
  varname(2)='K'
  Call getmeta(ncid,varname,coverout(:,:,20),arrsize)
  
  Write(6,*) "Reading walwtgg1/2/3/4"
  varname(1)='walwtgg1'
  varname(2)='K'
  Call getmeta(ncid,varname,coverout(:,:,21),arrsize)
  varname(1)='walwtgg2'
  varname(2)='K'
  Call getmeta(ncid,varname,coverout(:,:,22),arrsize)
  varname(1)='walwtgg3'
  varname(2)='K'
  Call getmeta(ncid,varname,coverout(:,:,23),arrsize)
  varname(1)='walwtgg4'
  varname(2)='K'
  Call getmeta(ncid,varname,coverout(:,:,24),arrsize)
  
  Write(6,*) "Reading roadtgg1/2/3/4"
  varname(1)='roadtgg1'
  varname(2)='K'
  Call getmeta(ncid,varname,coverout(:,:,25),arrsize)
  varname(1)='roadtgg2'
  varname(2)='K'
  Call getmeta(ncid,varname,coverout(:,:,26),arrsize)
  varname(1)='roadtgg3'
  varname(2)='K'
  Call getmeta(ncid,varname,coverout(:,:,27),arrsize)
  varname(1)='roadtgg4'
  varname(2)='K'
  Call getmeta(ncid,varname,coverout(:,:,28),arrsize)
End If

If (nf_inq_varid(ncid,'urbnsmc',varid).EQ.nf_noerr) Then
  Write(6,*) "Reading urbansm"
  varname(1)='urbnsmc'
  varname(2)='m3/m3'
  Call getmeta(ncid,varname,coverout(:,:,29),arrsize)
  varname(1)='urbnsmr'
  varname(2)='m3/m3'
  Call getmeta(ncid,varname,coverout(:,:,30),arrsize)
End If

If (nf_inq_varid(ncid,'snd',varid).EQ.nf_noerr) Then
  Write(6,*) "Reading snd"
  varname(1)='snd'
  varname(2)='mm'
  Call getmeta(ncid,varname,coverout(:,:,31),arrsize)
End If

If (nf_inq_varid(ncid,'smass1',varid).EQ.nf_noerr) Then
  Write(6,*) "Reading smass1/2/3"
  varname(1)='smass1'
  varname(2)='K'
  Call getmeta(ncid,varname,coverout(:,:,32),arrsize)
  varname(1)='smass2'
  varname(2)='K'
  Call getmeta(ncid,varname,coverout(:,:,33),arrsize)
  varname(1)='smass3'
  varname(2)='K'
  Call getmeta(ncid,varname,coverout(:,:,34),arrsize)

  Write(6,*) "Reading ssdn1/2/3"
  varname(1)='ssdn1'
  varname(2)='K'
  Call getmeta(ncid,varname,coverout(:,:,35),arrsize)
  varname(1)='ssdn2'
  varname(2)='K'
  Call getmeta(ncid,varname,coverout(:,:,36),arrsize)
  varname(1)='ssdn3'
  varname(2)='K'
  Call getmeta(ncid,varname,coverout(:,:,37),arrsize)

  Write(6,*) "Reading wbice1/2/3/4/5/6"
  varname(1)='wbice1'
  varname(2)='m3/m3'
  Call getmeta(ncid,varname,coverout(:,:,41),arrsize)
  varname(1)='wbice2'
  varname(2)='m3/m3'
  Call getmeta(ncid,varname,coverout(:,:,42),arrsize)
  varname(1)='wbice3'
  varname(2)='m3/m3'
  Call getmeta(ncid,varname,coverout(:,:,43),arrsize)
  varname(1)='wbice4'
  varname(2)='m3/m3'
  Call getmeta(ncid,varname,coverout(:,:,44),arrsize)
  varname(1)='wbice5'
  varname(2)='m3/m3'
  Call getmeta(ncid,varname,coverout(:,:,45),arrsize)
  varname(1)='wbice6'
  varname(2)='m3/m3'
  Call getmeta(ncid,varname,coverout(:,:,46),arrsize)

  Write(6,*) "Reading snage"
  varname(1)='snage'
  varname(2)='none'
  Call getmeta(ncid,varname,coverout(:,:,47),arrsize)

  Write(6,*) "Reading sflag"
  varname(1)='sflag'
  varname(2)='none'
  Call getmeta(ncid,varname,coverout(:,:,48),arrsize)
End If

coverout(:,:,34:36)=272.2
If (nf_inq_varid(ncid,'tggsn1',varid).EQ.nf_noerr) Then
  Write(6,*) "Reading tggsn1/2/3"
  varname(1)='tggsn1'
  varname(2)='K'
  Call getmeta(ncid,varname,coverout(:,:,38),arrsize)
  varname(1)='tggsn2'
  varname(2)='K'
  Call getmeta(ncid,varname,coverout(:,:,39),arrsize)
  varname(1)='tggsn3'
  varname(2)='K'
  Call getmeta(ncid,varname,coverout(:,:,40),arrsize)
End If

If (nf_inq_varid(ncid,'swater',varid).EQ.nf_noerr) Then
  Write(6,*) "Reading swater"
  varname(1)='swater'
  varname(2)='mm'
  Call getmeta(ncid,varname,coverout(:,:,49),arrsize)
End If

If (nf_inq_varid(ncid,'ssalin',varid).EQ.nf_noerr) Then
  Write(6,*) "Reading ssalin"
  varname(1)='ssalin'
  varname(2)='PSU'
  Call getmeta(ncid,varname,coverout(:,:,50),arrsize)
End If

If (nf_inq_varid(ncid,'fracice',varid).EQ.nf_noerr) Then
  Write(6,*) "Reading fracice"
  varname(1)='fracice'
  varname(2)='none'
  Call getmeta(ncid,varname,coverout(:,:,51),arrsize)
  Write(6,*) "Reading siced"
  varname(1)='siced'
  varname(2)='m'
  Call getmeta(ncid,varname,coverout(:,:,52),arrsize)
End If

if (nf_inq_varid(ncid,'ocndepth',varid).EQ.nf_noerr) then
  write(6,*) "Reading ocndepth"
  varname(1)='ocndepth'
  varname(2)='m'
  call getmeta(ncid,varname,coverout(:,:,53),arrsize)
  varname(1)='ocheight'
  varname(2)='m'
  call getmeta(ncid,varname,coverout(:,:,54),arrsize)
  write(6,*) "Reading tgg,sal,uoc,voc"
  coverout(:,:,55:60)=coverout(:,:,7:12)
  varname(2)='K'
  do i=7,wlev
    write(chr,'(I2.2)') i
    varname(1)='tgg'//chr
    call getmeta(ncid,varname,coverout(:,:,54+i),arrsize)
  end do
  varname(2)='PSU'
  do i=1,wlev
    write(chr,'(I2.2)') i
    varname(1)='sal'//chr
    call getmeta(ncid,varname,coverout(:,:,54+wlev+i),arrsize)
  end do
  varname(2)='m/s'
  do i=1,wlev
    varname(1)='uoc'//chr
    call getmeta(ncid,varname,coverout(:,:,54+2*wlev+i),arrsize)
    varname(1)='voc'//chr
    call getmeta(ncid,varname,coverout(:,:,54+3*wlev+i),arrsize)          
  end do
  write(6,*) "Reading tggsn4,sto"
  coverout(:,:,55+4*wlev)=coverout(:,:,38) ! tggsn1
  coverout(:,:,56+4*wlev)=coverout(:,:,39) ! tggsn2
  coverout(:,:,57+4*wlev)=coverout(:,:,40) ! tggsn3
  varname(1)='tggsn4'
  varname(2)='K'
  call getmeta(ncid,varname,coverout(:,:,58+4*wlev),arrsize)
  coverout(:,:,59+4*wlev)=coverout(:,:,31) ! snd    
  varname(1)='sto'
  varname(2)='J/m2'
  call getmeta(ncid,varname,coverout(:,:,60+4*wlev),arrsize)
end if

if (nf_inq_varid(ncid,'uic',varid).eq.nf_noerr) then
  write(6,*) "Reading uic,vic"
  varname(1)='uic'
  varname(2)='m/s'
  call getmeta(ncid,varname,coverout(:,:,61+4*wlev),arrsize)
  varname(1)='vic'
  varname(2)='m/s'
  call getmeta(ncid,varname,coverout(:,:,62+4*wlev),arrsize)
end if
if (nf_inq_varid(ncid,'icesal',varid).eq.nf_noerr) then
  write(6,*) "Reading icesal"
  varname(1)='icesal'
  varname(2)='PSU'
  call getmeta(ncid,varname,coverout(:,:,63+4*wlev),arrsize)
end if

ncstatus=nf_close(ncid)

! clean up input data
Do ilon=1,ncsize(1)
  Do ilat=1,ncsize(2)
    If (any(coverout(ilon,ilat,1:6).lt.0.)) coverout(ilon,ilat,1)=0.
    If (any((coverout(ilon,ilat,7:12).gt.400.).or.(coverout(ilon,ilat,7:12).lt.200.))) coverout(ilon,ilat,1)=0.
    If (any((coverout(ilon,ilat,51:52).lt.0.))) coverout(ilon,ilat,51)=0.
    if (any((coverout(ilon,ilat,55:54+wlev).lt.200.).or.(coverout(ilon,ilat,55:54+wlev).gt.400.))) coverout(ilon,ilat,53)=0.
    if (any((coverout(ilon,ilat,55+wlev:54+2*wlev).lt.0.).or.(coverout(ilon,ilat,55+wlev:54+2*wlev).gt.100.))) &
        coverout(ilon,ilat,53)=0.
    if (coverout(ilon,ilat,1).eq.0.) then
      coverout(ilon,ilat,1:50)=0.
    end if
    if (coverout(ilon,ilat,51).eq.0.) then
      coverout(ilon,ilat,51:52)=0.
    end if
    if (coverout(ilon,ilat,53).eq.0.) then
      coverout(ilon,ilat,53:63+4*wlev)=0.
    end if
  End Do
End Do

! Bin data (not ocean points since each ocean point has a different vertical level structure)
Write(6,*) 'Bin'
Do ilat=1,ncsize(2)
  aglat=smlonlat(2,1)+real(ilat-1)*(smlonlat(2,2)-smlonlat(2,1))/(smlonlat(2,3)-1.)
  Do ilon=1,ncsize(1)

    If ((coverout(ilon,ilat,1).ne.0.).or.(coverout(ilon,ilat,51).ne.0.)) then

      aglon=smlonlat(1,1)+real(ilon-1)*(smlonlat(1,2)-smlonlat(1,1))/(smlonlat(1,3)-1.)

      Call lltoijmod(aglon,aglat,alci,alcj,nface)
      lci = nint(alci)
      lcj = nint(alcj)
      lcj = lcj+nface*ccdim(1)
  
      dataout(lci,lcj,1:52)=dataout(lci,lcj,1:52)+coverout(ilon,ilat,1:52)
      if (coverout(ilon,ilat,1).ne.0.) countn(lci,lcj)=countn(lci,lcj)+1
      if (coverout(ilon,ilat,51).ne.0.) countm(lci,lcj)=countm(lci,lcj)+1
    End if
    
  End Do
End Do

! assign null values over land/sea
Do lcj=1,ccdim(2)
  Do lci=1,ccdim(1)
    if (1-nint(lsmask(lci,lcj)).EQ.0) then
      dataout(lci,lcj,1:50)=0.
      countn(lci,lcj)=1
    end if
    if (1-nint(lsmask(lci,lcj)).EQ.1) then
      dataout(lci,lcj,51:52)=0.
      countm(lci,lcj)=1
      dataout(lci,lcj,53:63+4*wlev)=0.
      counto(lci,lcj)=1
    end If
  End Do
End Do

Write (formout,'("(",i3,"l1)" )') ccdim(1)
Do j=1,ccdim(2)
  Write(6,formout) (countn(:,j).NE.0).and.(countm(:,j).NE.0).and.(counto(:,j).ne.0)
End Do


! Interpolate/fill data
If (any(countn.LT.1).or.any(countm.lt.1).or.any(counto.lt.1)) then
  Write(6,*) 'Fill'  

  Do lcj=1,ccdim(2)
    Do lci=1,ccdim(1)
      aglon=rlld(lci,lcj,1)
      aglat=rlld(lci,lcj,2)
	      
      serlon=1.+(aglon-smlonlat(1,1))*(smlonlat(1,3)-1.)/(smlonlat(1,2)-smlonlat(1,1))
      serlat=1.+(aglat-smlonlat(2,1))*(smlonlat(2,3)-1.)/(smlonlat(2,2)-smlonlat(2,1))
      i=nint(serlon)
      j=nint(serlat)
      
      if (i.gt.ncsize(1)) i=i-ncsize(1)
      
      if (countn(lci,lcj).eq.0.and.coverout(i,j,1)/=0.) then
        dataout(lci,lcj,1:50)=coverout(i,j,1:50)
        countn(lci,lcj)=1
      end if
      
      if (countm(lci,lcj).eq.0) then
        dataout(lci,lcj,51:52)=coverout(i,j,51:52)
        countm(lci,lcj)=1
      end if
        
      if (counto(lci,lcj).eq.0.and.coverout(i,j,53)/=0.) then
        dataout(lci,lcj,53:63+4*wlev)=coverout(i,j,53:63+4*wlev)
        counto(lci,lcj)=1
      end if

    End Do
  End Do
    
  Write (formout,'("(",i3,"l1)" )') ccdim(1)
  Do j=1,ccdim(2)
    Write(6,formout) (countn(:,j).NE.0).and.(countm(:,j).NE.0).and.(counto(:,j).ne.0)
  End Do

End if

Deallocate(coverout)

! nearest nbr
If (Any(countn.LT.1)) then
  Write(6,*) "Use near nbr for land"
  allocate(sermask(ccdim(1),ccdim(2)))
  sermask(:,:)=(countn(:,:).GT.0).AND.(dataout(:,:,1).ne.0.)
  If (.NOT.any(sermask)) Then
    Write(6,*) "ERROR: Cannot find valid soil data"
    Stop
  End If
  Do lcj=1,ccdim(2)
    Do lci=1,ccdim(1)
      If (countn(lci,lcj).EQ.0) then
        call findnear(pxy,lci,lcj,sermask,rlld,ccdim)
	    dataout(lci,lcj,1:50)=dataout(pxy(1),pxy(2),1:50)
	    countn(lci,lcj)=countn(pxy(1),pxy(2))
      End if
    End do
  End do  
  deallocate(sermask)

  Write (formout,'("(",i3,"l1)" )') ccdim(1)
  Do j=1,ccdim(2)
    Write(6,formout) (countn(:,j).NE.0).and.(countm(:,j).NE.0).and.(counto(:,j).ne.0)
  End Do
  
End if

If (Any(counto.LT.1)) then
  Write(6,*) "Use near nbr for water"
  allocate(sermask(ccdim(1),ccdim(2)))
  sermask(:,:)=(counto(:,:).GT.0).AND.(dataout(:,:,53).ne.0.)
  If (.NOT.any(sermask)) Then
    Write(6,*) "WARN: Cannot find valid ocean data"
    dataout(:,:,53)=0.
    counto=1
  Else
    Do lcj=1,ccdim(2)
      Do lci=1,ccdim(1)
        If (counto(lci,lcj).EQ.0) then
          call findnear(pxy,lci,lcj,sermask,rlld,ccdim)
	      dataout(lci,lcj,53:63+4*wlev)=dataout(pxy(1),pxy(2),53:63+4*wlev)
	      counto(lci,lcj)=counto(pxy(1),pxy(2))
        End if
      End do
    End do  
  end if
  deallocate(sermask)

  Write (formout,'("(",i3,"l1)" )') ccdim(1)
  Do j=1,ccdim(2)
    Write(6,formout) (countn(:,j).NE.0).and.(countm(:,j).NE.0).and.(counto(:,j).ne.0)
  End Do
  
End if

If (Any(countn.LT.1).or.Any(countm.LT.1).or.any(counto.lt.1)) then
  Write(6,*) 'ERROR: unassigned points detected'
  Stop
End If

Do k=1,50
  dataout(:,:,k)=dataout(:,:,k)/Real(countn)
End do
Do k=51,52
  dataout(:,:,k)=dataout(:,:,k)/Real(countm)
End do
Do k=53,63+4*wlev
  dataout(:,:,k)=dataout(:,:,k)/Real(counto)
End do

do k=1,6
  where(dataout(:,:,53).ne.0.) ! water
    dataout(:,:,6+k)=dataout(:,:,54+k)
  end where
end do
do k=1,3
  where(dataout(:,:,53).ne.0.) ! water
    dataout(:,:,37+k)=dataout(:,:,54+4*wlev+k)
  end where
end do
where(dataout(:,:,53).ne.0.) ! water
  dataout(:,:,31)=dataout(:,:,59+4*wlev)
end where

do lcj=1,ccdim(2)
  do lci=1,ccdim(1)
    if (dataout(lci,lcj,53).ne.0.) then
      ! rotate current to CC coordinates
      zonx=                     -polenz*xyz(lci,lcj,2)
      zony=polenz*xyz(lci,lcj,1)-polenx*xyz(lci,lcj,3)
      zonz=polenx*xyz(lci,lcj,2)
      den=sqrt( max(zonx**2 + zony**2 + zonz**2,1.e-7) )  ! allow for poles
      costh= (zonx*axyz(lci,lcj,1)+zony*axyz(lci,lcj,2)+zonz*axyz(lci,lcj,3))/den
      sinth=-(zonx*bxyz(lci,lcj,1)+zony*bxyz(lci,lcj,2)+zonz*bxyz(lci,lcj,3))/den
      uzon(:)=dataout(lci,lcj,55+2*wlev:54+3*wlev)
      vmer(:)=dataout(lci,lcj,55+3*wlev:54+4*wlev)
      dataout(lci,lcj,55+2*wlev:54+3*wlev)= costh*uzon(:)+sinth*vmer(:)
      dataout(lci,lcj,55+3*wlev:54+4*wlev)=-sinth*uzon(:)+costh*vmer(:)
      uzon(1)=dataout(lci,lcj,61+4*wlev)
      vmer(1)=dataout(lci,lcj,62+4*wlev)
      dataout(lci,lcj,61+4*wlev)= costh*uzon(1)+sinth*vmer(1)
      dataout(lci,lcj,62+4*wlev)=-sinth*uzon(1)+costh*vmer(1)
    end if
  end do
end do

Return
End


