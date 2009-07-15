! These subroutines read/bin/interpolate soil moisture data for CCAM.
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine reads soil moisture data (i.e., assume resolution
! is >=10km so that we simply stream the data).  Replace with ecoveg
! style aggregate-on-read if there is a performance or memory issue.

Subroutine getsmdata(climfile,dataout,lsmask,grid,tlld,ccdim,wlev)

Use ccinterp

Implicit None

Include "netcdf.inc"

integer, intent(in) :: wlev
Integer, dimension(1:2), intent(in) :: ccdim
Integer, dimension(1:ccdim(1),1:ccdim(2)) :: countn,countm,counto
Integer, dimension(1:4) :: ncsize
Integer, dimension(1:4,1:2) :: arrsize
Integer, dimension(1:2) :: pxy
Integer ilat,ilon,lci,lcj,i,j,nface,k
Integer ncstatus,ncid,varid
Character(len=*), intent(in) :: climfile
Character*80, dimension(1:2) :: varname
Character*9 formout
character*2 chr
Real, dimension(1:ccdim(1),1:ccdim(2)), intent(in) :: lsmask
Real, dimension(1:ccdim(1),1:ccdim(2),1:45+4*wlev), intent(out) :: dataout
Real, dimension(1:ccdim(1),1:ccdim(2)), intent(in) :: grid
Real, dimension(1:ccdim(1),1:ccdim(2),1:2), intent(in) :: tlld
Real, dimension(1:ccdim(1),1:ccdim(2),1:2) :: rlld
Real, dimension(1:2,1:3) :: smlonlat
Real, dimension(:,:,:), allocatable :: coverout
Real aglon,aglat,alci,alcj,serlon,serlat
Real minlon,maxlon
Logical, dimension(:,:), allocatable :: sermask

dataout=0.
countn=0
countm=0
counto=0.

! Open nc moisture file

ncstatus=nf_open(climfile,nf_nowrite,ncid)
If (ncstatus.NE.nf_noerr) Then
  Write(6,*) "ERROR: Error opening NetCDF file ",trim(climfile)," (",ncstatus,")"
  Stop
End If

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
Allocate(coverout(1:arrsize(1,2),1:arrsize(2,2),1:45+4*wlev))
coverout=0.

! Read data
If (nf_inq_varid(ncid,'wetfrac1',varid).EQ.nf_noerr) Then
  Write(6,*) "Reading wetfrac1/2/3/4/5/6"
  varname=(/ 'wetfrac1', 'none' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,1),arrsize)
  varname=(/ 'wetfrac2', 'none' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,2),arrsize)
  varname=(/ 'wetfrac3', 'none' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,3),arrsize)
  varname=(/ 'wetfrac4', 'none' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,4),arrsize)
  varname=(/ 'wetfrac5', 'none' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,5),arrsize)
  varname=(/ 'wetfrac6', 'none' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,6),arrsize)
Else If ((nf_inq_varid(ncid,'wetfrac3',varid).EQ.nf_noerr).AND.(nf_inq_varid(ncid,'wetfrac5',varid).EQ.nf_noerr)) Then
  Write(6,*) "Reading wetfrac3/5"
  varname=(/ 'wetfrac3', 'none' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,1),arrsize)
  coverout(1:arrsize(1,2),:,2)=coverout(1:arrsize(1,2),:,1)
  coverout(1:arrsize(1,2),:,3)=coverout(1:arrsize(1,2),:,1)
  varname=(/ 'wetfrac5', 'none' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,4),arrsize)
  coverout(1:arrsize(1,2),:,5)=coverout(1:arrsize(1,2),:,4)
  coverout(1:arrsize(1,2),:,6)=coverout(1:arrsize(1,2),:,4)
Else If (nf_inq_varid(ncid,'wetfrac',varid).EQ.nf_noerr) Then
  Write(6,*) "Reading wetfrac"
  varname=(/ 'wetfrac', 'none' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,1),arrsize)
  coverout(1:arrsize(1,2),:,2)=coverout(1:arrsize(1,2),:,1)
  coverout(1:arrsize(1,2),:,3)=coverout(1:arrsize(1,2),:,1)
  coverout(1:arrsize(1,2),:,4)=coverout(1:arrsize(1,2),:,1)
  coverout(1:arrsize(1,2),:,5)=coverout(1:arrsize(1,2),:,1)
  coverout(1:arrsize(1,2),:,6)=coverout(1:arrsize(1,2),:,1)
Else
  Write(6,*) "ERROR: Cannot find valid soil mositure data in ",trim(climfile)
  Stop
End If

If (nf_inq_varid(ncid,'tgg1',varid).EQ.nf_noerr) Then
  Write(6,*) "Reading tgg1/2/3/4/5/6"
  varname=(/ 'tgg1', 'K' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,7),arrsize)
  varname=(/ 'tgg2', 'K' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,8),arrsize)
  varname=(/ 'tgg3', 'K' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,9),arrsize)
  varname=(/ 'tgg4', 'K' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,10),arrsize)
  varname=(/ 'tgg5', 'K' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,11),arrsize)
  varname=(/ 'tgg6', 'K' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,12),arrsize)
Else If ((nf_inq_varid(ncid,'tgg3',varid).EQ.nf_noerr).AND.(nf_inq_varid(ncid,'tgg5',varid).EQ.nf_noerr)) Then
  Write(6,*) "Reading tgg3/5"
  varname=(/ 'tgg3', 'K' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,7),arrsize)
  coverout(1:arrsize(1,2),:,8)=coverout(1:arrsize(1,2),:,7)
  coverout(1:arrsize(1,2),:,9)=coverout(1:arrsize(1,2),:,7)
  varname=(/ 'tgg5', 'K' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,10),arrsize)
  coverout(1:arrsize(1,2),:,11)=coverout(1:arrsize(1,2),:,10)
  coverout(1:arrsize(1,2),:,12)=coverout(1:arrsize(1,2),:,10)
Else
  Write(6,*) "ERROR: Cannot find valid soil temperature data in ",trim(climfile)
  Stop
End If

If (nf_inq_varid(ncid,'rooftgg1',varid).EQ.nf_noerr) Then
  Write(6,*) "Reading rooftgg1/2/3"
  varname=(/ 'rooftgg1', 'K' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,13),arrsize)
  varname=(/ 'rooftgg2', 'K' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,14),arrsize)
  varname=(/ 'rooftgg3', 'K' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,15),arrsize)

  Write(6,*) "Reading waletgg1/2/3"
  varname=(/ 'waletgg1', 'K' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,16),arrsize)
  varname=(/ 'waletgg2', 'K' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,17),arrsize)
  varname=(/ 'waletgg3', 'K' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,18),arrsize)

  Write(6,*) "Reading walwtgg1/2/3"
  varname=(/ 'walwtgg1', 'K' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,19),arrsize)
  varname=(/ 'walwtgg2', 'K' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,20),arrsize)
  varname=(/ 'walwtgg3', 'K' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,21),arrsize)

  Write(6,*) "Reading roadtgg1/2/3"
  varname=(/ 'roadtgg1', 'K' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,22),arrsize)
  varname=(/ 'roadtgg2', 'K' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,23),arrsize)
  varname=(/ 'roadtgg3', 'K' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,24),arrsize)
End If

If (nf_inq_varid(ncid,'snd',varid).EQ.nf_noerr) Then
  Write(6,*) "Reading snd"
  varname=(/ 'snd', 'mm' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,25),arrsize)
End If

If (nf_inq_varid(ncid,'smass1',varid).EQ.nf_noerr) Then
  Write(6,*) "Reading smass1/2/3"
  varname=(/ 'smass1', 'K' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,26),arrsize)
  varname=(/ 'smass2', 'K' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,27),arrsize)
  varname=(/ 'smass3', 'K' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,28),arrsize)

  Write(6,*) "Reading ssdn1/2/3"
  varname=(/ 'ssdn1', 'K' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,29),arrsize)
  varname=(/ 'ssdn2', 'K' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,30),arrsize)
  varname=(/ 'ssdn3', 'K' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,31),arrsize)

  Write(6,*) "Reading tggsn1/2/3"
  varname=(/ 'tggsn1', 'K' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,32),arrsize)
  varname=(/ 'tggsn2', 'K' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,33),arrsize)
  varname=(/ 'tggsn3', 'K' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,34),arrsize) 

  Write(6,*) "Reading wbice1/2/3/4/5/6"
  varname=(/ 'wbice1', 'm3/m3' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,35),arrsize)
  varname=(/ 'wbice2', 'm3/m3' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,36),arrsize)
  varname=(/ 'wbice3', 'm3/m3' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,37),arrsize)
  varname=(/ 'wbice4', 'm3/m3' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,38),arrsize)
  varname=(/ 'wbice5', 'm3/m3' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,39),arrsize)
  varname=(/ 'wbice6', 'm3/m3' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,40),arrsize)

  Write(6,*) "Reading snage"
  varname=(/ 'snage', 'none' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,41),arrsize)

  Write(6,*) "Reading sflag"
  varname=(/ 'sflag', 'none' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,42),arrsize)
End If

If (nf_inq_varid(ncid,'fracice',varid).EQ.nf_noerr) Then
  Write(6,*) "Reading fracice"
  varname=(/ 'fracice', 'none' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,43),arrsize)
  Write(6,*) "Reading siced"
  varname=(/ 'siced', 'm' /)
  Call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,44),arrsize)
End If

if (nf_inq_varid(ncid,'ocndepth',varid).EQ.nf_noerr) then
  write(6,*) "Reading ocndepth"
  varname=(/ 'ocndepth', 'm' /)
  call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,45),arrsize)
  write(6,*) "Reading tgg,sal,uoc,voc"
  coverout(1:arrsize(1,2),:,46:51)=coverout(1:arrsize(1,2),:,7:12)
  do i=7,wlev
    write(chr,'(I2.2)') i
    varname(1)='tgg'//chr
    varname(2)='K'
    call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,45+i),arrsize)
  end do
  do i=1,wlev
    write(chr,'(I2.2)') i
    varname(1)='sal'//chr
    varname(2)='PSU'
    call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,45+wlev+i),arrsize)
    varname(1)='uoc'//chr
    varname(2)='m/s'
    call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,45+2*wlev+i),arrsize)
    varname(1)='voc'//chr
    varname(2)='m/s'
    call getmeta(ncid,varname,coverout(1:arrsize(1,2),:,45+3*wlev+i),arrsize)          
  end do
end if

ncstatus=nf_close(ncid)

Do ilon=1,ncsize(1)
  Do ilat=1,ncsize(2)
    If (any((coverout(ilon,ilat,1:6).GT.2).OR.(coverout(ilon,ilat,1:6).lt.0.))) coverout(ilon,ilat,1)=0.
    If (any((coverout(ilon,ilat,7:12).GT.350.).OR.(coverout(ilon,ilat,7:12).LT.250.))) coverout(ilon,ilat,1)=0.
    If (any((coverout(ilon,ilat,43:44).lt.0.))) coverout(ilon,ilat,43)=0.
    if (any((coverout(ilon,ilat,46:45+wlev).lt.250.).or.(coverout(ilon,ilat,46:45+wlev).gt.350.))) coverout(ilon,ilat,45)=0.
    if (any((coverout(ilon,ilat,46+wlev:45+2*wlev).lt.0.).or.(coverout(ilon,ilat,46+wlev:45+2*wlev).gt.100.))) coverout(ilon,ilat,45)=0.
    if (any((coverout(ilon,ilat,46+2*wlev:45+3*wlev).lt.-10.).or.(coverout(ilon,ilat,46+2*wlev:45+3*wlev).gt.10.))) coverout(ilon,ilat,45)=0.
    if (any((coverout(ilon,ilat,46+3*wlev:45+4*wlev).lt.-10.).or.(coverout(ilon,ilat,46+3*wlev:45+4*wlev).gt.10.))) coverout(ilon,ilat,45)=0.
    if (coverout(ilon,ilat,1).eq.0.) coverout(ilon,ilat,1:42)=0.
    if (coverout(ilon,ilat,43).eq.0.) coverout(ilon,ilat,43:44)=0.
    if (coverout(ilon,ilat,45).eq.0.) coverout(ilon,ilat,45:45+4*wlev)=0.
  End Do
End Do

! Bin data (not ocean points)
Write(6,*) 'Bin'

Do ilat=1,ncsize(2)
  aglat=smlonlat(2,1)+real(ilat-1)*(smlonlat(2,2)-smlonlat(2,1))/(smlonlat(2,3)-1.)
  Do ilon=1,ncsize(1)

    If ((coverout(ilon,ilat,1).GT.0.).or.(coverout(ilon,ilat,43).GT.0.)) then

      aglon=smlonlat(1,1)+real(ilon-1)*(smlonlat(1,2)-smlonlat(1,1))/(smlonlat(1,3)-1.)

      Call lltoijmod(aglon,aglat,alci,alcj,nface)
      lci = nint(alci)
      lcj = nint(alcj)
      lcj = lcj+nface*ccdim(1)
  
      dataout(lci,lcj,1:44)=dataout(lci,lcj,1:44)+coverout(ilon,ilat,1:44)
      if (coverout(ilon,ilat,1).GT.0.) countn(lci,lcj)=countn(lci,lcj)+1
      if (coverout(ilon,ilat,43).GT.0.) countm(lci,lcj)=countm(lci,lcj)+1
    End if
    
  End Do
End Do


Do lcj=1,ccdim(2)
  Do lci=1,ccdim(1)
    if (1-nint(lsmask(lci,lcj)).EQ.0) then
      dataout(lci,lcj,1:42)=0.
      countn(lci,lcj)=1
    end if
    if (1-nint(lsmask(lci,lcj)).EQ.1) then
      dataout(lci,lcj,43:44)=0.
      countm(lci,lcj)=1
      dataout(lci,lcj,45:45+4*wlev)=0.
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

      if (countn(lci,lcj).eq.0.and.coverout(i,j,1).gt.0.) then
        dataout(lci,lcj,1:42)=coverout(i,j,1:42)
        countn(lci,lcj)=1
      end if
      
      if (countm(lci,lcj).eq.0) then
        dataout(lci,lcj,43:44)=coverout(i,j,43:44)
        countm(lci,lcj)=1
      end if
        
      if (counto(lci,lcj).eq.0.and.coverout(i,j,45).gt.0.) then
        dataout(lci,lcj,45:45+4*wlev)=coverout(i,j,45:45+4*wlev)
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
  sermask(:,:)=(countn(:,:).GT.0).AND.(dataout(:,:,1).GT.0.)
  If (.NOT.any(sermask)) Then
    Write(6,*) "ERROR: Cannot find valid soil data"
    Stop
  End If
  Do lcj=1,ccdim(2)
    Do lci=1,ccdim(1)
      If (countn(lci,lcj).EQ.0) then
        call findnear(pxy,lci,lcj,sermask,rlld,ccdim)
	    dataout(lci,lcj,1:42)=dataout(pxy(1),pxy(2),1:42)
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
  sermask(:,:)=(counto(:,:).GT.0).AND.(dataout(:,:,45).GT.0.)
  If (.NOT.any(sermask)) Then
    Write(6,*) "WARN: Cannot find valid ocean data"
    dataout(:,:,45)=0.
    counto=1
  Else
    Do lcj=1,ccdim(2)
      Do lci=1,ccdim(1)
        If (counto(lci,lcj).EQ.0) then
          call findnear(pxy,lci,lcj,sermask,rlld,ccdim)
	      dataout(lci,lcj,45:45+4*wlev)=dataout(pxy(1),pxy(2),45:45+4*wlev)
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

Do k=1,42
  dataout(:,:,k)=dataout(:,:,k)/Real(countn)
End do
Do k=43,44
  dataout(:,:,k)=dataout(:,:,k)/Real(countm)
End do
Do k=45,45+4*wlev
  dataout(:,:,k)=dataout(:,:,k)/Real(counto)
End do

do k=1,6
  where(dataout(:,:,45).gt.0.) ! water
    dataout(:,:,6+k)=dataout(:,:,45+k)
  end where
end do

Return
End


