Program smclim

! This code interpolates a soil moisture climatology onto the Conformal
! cubic grid (i.e., for stand-alone NWP forecasting)

Implicit None

Integer :: nopts
Character*80, dimension(:,:), allocatable :: options

Write(6,*) 'SMCLIM - interpolate/bin climatology soil moisture/temperature on'
Write(6,*) '         conformal-cubic grid (MAR-13)'

! Read switches
nopts=3
Allocate (options(nopts,2))
options(:,1) = (/ '-t', '-c', '-o' /)
options(:,2) = ''

Call readswitch(options,nopts)
Call defaults(options,nopts)

Call createsm(options,nopts)

Deallocate(options)

Stop
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine displays the help message
!

Subroutine help()

Implicit None

Write(6,*)
Write(6,*) "Usage:"
Write(6,*) "  smclim [-t topofile] -c climfile -o outputfile"
Write(6,*)
Write(6,*) "Options:"
Write(6,*) "  topofile      = input topography filename."
Write(6,*) "                  (DEFAULT = topout)"
Write(6,*) "  climfile      = input climatology filename (i.e.,"
Write(6,*) "                  appropriate for the current month)."
Write(6,*) "  outputfile    = output analysis filename"
Write(6,*) "                  (must be pre-existing)"
Write(6,*)
Stop

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine determins the default values for the switches
!

Subroutine defaults(options,nopts)

Implicit None

Integer nopts
Character(len=*), dimension(nopts,2), intent(inout) :: options
Integer climfile,outfile,topofile
Integer locate

climfile=locate('-c',options(:,1),nopts)
outfile=locate('-o',options(:,1),nopts)
topofile=locate('-t',options(:,1),nopts)

If (options(topofile,2).EQ.'') then
  options(topofile,2)='topout'
End if

If (options(climfile,2).EQ.'') then
  Write(6,*) "ERROR: No soil moisture file specified"
  Stop
End if

If (options(outfile,2).EQ.'') then
  options(outfile,2)='sm'
End if

Return
End


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine processes the soil moisture data
!

Subroutine createsm(options,nopts)

Use ccinterp

Implicit None

include 'netcdf.inc'

Integer, intent(in) :: nopts
Integer, dimension(1:2) :: ccdim
integer :: wlev = 20
integer ncid,ncstatus,olevid
Character(len=*), dimension(nopts,2), intent(in) :: options
Character*80 climfile,outfile,topofile,bintmp,returnoption
Character*9 formout
Character*47 header
Real, dimension(:,:,:), allocatable :: rlld,xyz,axyz,bxyz
Real, dimension(:,:), allocatable :: lsdata,gridout
Real, dimension(:,:,:), allocatable :: rawdata
Real, dimension(1:2) :: lonlat
Real schmidt,dsx,ds

climfile=returnoption('-c',options,nopts)
outfile=returnoption('-o',options,nopts)
topofile=returnoption('-t',options,nopts)

! Get CC grid coordinates
Call readtopography(1,topofile,ccdim,lonlat,schmidt,dsx,header)
Write(6,*) "Dimension : ",ccdim
Write(6,*) "lon0,lat0 : ",lonlat
Write(6,*) "Schmidt   : ",schmidt

! Open nc moisture file
ncstatus=nf_open(climfile,nf_nowrite,ncid)
If (ncstatus.NE.nf_noerr) Then
  Write(6,*) "ERROR: Error opening NetCDF file ",trim(climfile)," (",ncstatus,")"
  Stop
End If
ncstatus = nf_inq_dimid(ncid,'olev',olevid)
if (ncstatus.eq.nf_noerr) Then
  ncstatus = nf_inq_dimlen(ncid,olevid,wlev)
  If (ncstatus.NE.nf_noerr) Then
    Write(6,*) "ERROR: Cannot determine olev len (",ncstatus,")"
    Stop
  End If
else
  write(6,*) "WARN: Cannot locate olev in input file"
  write(6,*) "Using wlev ",wlev
end if

Allocate(gridout(1:ccdim(1),1:ccdim(2)),rlld(1:ccdim(1),1:ccdim(2),1:2))
Allocate(lsdata(1:ccdim(1),1:ccdim(2)),rawdata(1:ccdim(1),1:ccdim(2),1:63+4*wlev))
allocate(xyz(ccdim(1),ccdim(2),3),axyz(ccdim(1),ccdim(2),3))
allocate(bxyz(ccdim(1),ccdim(2),3))

! Determine lat/lon to CC mapping
Call getcc(rlld,gridout,xyz,axyz,bxyz,ccdim,lonlat,schmidt,ds)

! Get land/sea mask
Write(6,*) "Read land/sea mask"
Call gettopols(1,topofile,lsdata,ccdim)

! Get climate soil moisture
Call getsmdata(ncid,rawdata,lsdata,gridout,rlld,ccdim,wlev,lonlat,xyz,axyz,bxyz)

! Write output soil moisture data (currently disabled)
Call storesm(outfile,rawdata,ccdim,wlev)

Deallocate(gridout,rlld,lsdata,rawdata,xyz,axyz,bxyz)

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine appends the soil moisture data to the analysis file
!

Subroutine storesm(outfile,outdata,ccdim,wlev)

Implicit None

Include "netcdf.inc"

integer, intent(in) :: wlev
integer olevid
Integer, dimension(1:2), intent(in) :: ccdim
Integer, dimension(0:4) :: ncidarr
Integer, dimension(1:4) :: dimnum
Integer ncstatus,varid,dimtype,i
Character(len=*), intent(in) :: outfile
Character*80, dimension(1:3) :: elemdesc
Character*80 outname
character*2 chr
Real, dimension(1:ccdim(1),1:ccdim(2),1:63+4*wlev), intent(in) :: outdata
Real, dimension(1:ccdim(1),1:ccdim(2)) :: outdatab

Write(6,*) "Appending data to ",trim(outfile)

ncidarr=0

ncstatus=nf_open(outfile,nf_write,ncidarr(0))
If (ncstatus.NE.nf_noerr) Then
  Write(6,*) "ERROR: Error opening NetCDF file ",trim(outfile)," (",ncstatus,")"
  Stop
End If

Call ncfinddimid(ncidarr(0),"lon",outname,ncidarr(1))
Call ncfinddimid(ncidarr(0),"lat",outname,ncidarr(2))
Call ncfinddimid(ncidarr(0),"lev",outname,ncidarr(3))
Call ncfinddimid(ncidarr(0),"time",outname,ncidarr(4))

dimnum(1:2)=ccdim(1:2)
dimnum(3:4)=1

ncstatus=nf_redef(ncidarr(0))

ncstatus=nf_def_dim(ncidarr(0),"olev",wlev,olevid)
If (ncstatus /= nf_noerr) Then
  if (ncstatus.eq.nf_enameinuse) then
    write(6,*) "WARN: Dimension olev already exists"
  else
    Write(6,*) "ERROR: Error defining dim in NetCDF file (",ncstatus,"): olev"
    Stop
  end if
End If

elemdesc(1)='wetfrac1'
elemdesc(2)='Wetness fraction layer 1'
elemdesc(3)='none'
Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.0002,0.)
elemdesc(1)='wetfrac2'
elemdesc(2)='Wetness fraction layer 2'
elemdesc(3)='none'
Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.0002,0.)
elemdesc(1)='wetfrac3'
elemdesc(2)='Wetness fraction layer 3'
elemdesc(3)='none'
Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.0002,0.)
elemdesc(1)='wetfrac4'
elemdesc(2)='Wetness fraction layer 4'
elemdesc(3)='none'
Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.0002,0.)
elemdesc(1)='wetfrac5'
elemdesc(2)='Wetness fraction layer 5'
elemdesc(3)='none'
Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.0002,0.)
elemdesc(1)='wetfrac6'
elemdesc(2)='Wetness fraction layer 6'
elemdesc(3)='none'
Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.0002,0.)

elemdesc(1)='tgg1'
elemdesc(2)='Soil temperature lev 1'
elemdesc(3)='K'
Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
elemdesc(1)='tgg2'
elemdesc(2)='Soil temperature lev 2'
elemdesc(3)='K'
Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
elemdesc(1)='tgg3'
elemdesc(2)='Soil temperature lev 3'
elemdesc(3)='K'
Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
elemdesc(1)='tgg4'
elemdesc(2)='Soil temperature lev 4'
elemdesc(3)='K'
Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
elemdesc(1)='tgg5'
elemdesc(2)='Soil temperature lev 5'
elemdesc(3)='K'
Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
elemdesc(1)='tgg6'
elemdesc(2)='Soil temperature lev 6'
elemdesc(3)='K'
Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)

if (any(outdata(:,:,13).ne.0.)) then
  elemdesc(1)='rooftgg1'
  elemdesc(2)='roof temperature lev 1'
  elemdesc(3)='K'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
  elemdesc(1)='rooftgg2'
  elemdesc(2)='roof temperature lev 2'
  elemdesc(3)='K'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
  elemdesc(1)='rooftgg3'
  elemdesc(2)='roof temperature lev 3'
  elemdesc(3)='K'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
  elemdesc(1)='rooftgg4'
  elemdesc(2)='roof temperature lev 4'
  elemdesc(3)='K'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
  elemdesc(1)='waletgg1'
  elemdesc(2)='wall temperature lev 1'
  elemdesc(3)='K'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
  elemdesc(1)='waletgg2'
  elemdesc(2)='wall temperature lev 2'
  elemdesc(3)='K'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
  elemdesc(1)='waletgg3'
  elemdesc(2)='wall temperature lev 3'
  elemdesc(3)='K'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
  elemdesc(1)='waletgg4'
  elemdesc(2)='wall temperature lev 4'
  elemdesc(3)='K'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
  elemdesc(1)='walwtgg1'
  elemdesc(2)='wall temperature lev 1'
  elemdesc(3)='K'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
  elemdesc(1)='walwtgg2'
  elemdesc(2)='wall temperature lev 2'
  elemdesc(3)='K'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
  elemdesc(1)='walwtgg3'
  elemdesc(2)='wall temperature lev 3'
  elemdesc(3)='K'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
  elemdesc(1)='walwtgg4'
  elemdesc(2)='wall temperature lev 4'
  elemdesc(3)='K'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
  elemdesc(1)='roadtgg1'
  elemdesc(2)='road temperature lev 1'
  elemdesc(3)='K'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
  elemdesc(1)='roadtgg2'
  elemdesc(2)='road temperature lev 2'
  elemdesc(3)='K'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
  elemdesc(1)='roadtgg3'
  elemdesc(2)='road temperature lev 3'
  elemdesc(3)='K'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
  elemdesc(1)='roadtgg4'
  elemdesc(2)='road temperature lev 4'
  elemdesc(3)='K'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
end if

if (any(outdata(:,:,29).ne.0.)) then
  elemdesc(1)='urbnsmc'
  elemdesc(2)='urban canyon soil mositure'
  elemdesc(3)='m3/m3'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.0001,0.5)
  elemdesc(1)='urbnsmr'
  elemdesc(2)='urban roof soil mositure'
  elemdesc(3)='m3/m3'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.0001,0.5)
end if

if (any(outdata(:,:,31).ne.0.)) then
  elemdesc(1)='snd'
  elemdesc(2)='Snow depth (liquid water)'
  elemdesc(3)='mm'
  Call ncaddvargen(ncidarr,elemdesc,nf_float,3,varid,1.,0.)
end if

if (any(outdata(:,:,32).ne.0.)) then
  elemdesc(1)='smass1'
  elemdesc(2)='Snow mass lev 1'
  elemdesc(3)='K'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.01,250.)
  elemdesc(1)='smass2'
  elemdesc(2)='Snow mass lev 2'
  elemdesc(3)='K'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.01,250.)
  elemdesc(1)='smass3'
  elemdesc(2)='Snow mass lev 3'
  elemdesc(3)='K'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.01,250.)
  elemdesc(1)='ssdn1'
  elemdesc(2)='Snow density lev 1'
  elemdesc(3)='K'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.01,250.)
  elemdesc(1)='ssdn2'
  elemdesc(2)='Snow density lev 2'
  elemdesc(3)='K'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.01,250.)
  elemdesc(1)='ssdn3'
  elemdesc(2)='Snow density lev 3'
  elemdesc(3)='K'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.01,250.)
  elemdesc(1)='wbice1'
  elemdesc(2)='Soil ice lev 1'
  elemdesc(3)='m3/m3'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.0001,0.)  
  elemdesc(1)='wbice2'
  elemdesc(2)='Soil ice lev 2'
  elemdesc(3)='m3/m3'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.0001,0.)  
  elemdesc(1)='wbice3'
  elemdesc(2)='Soil ice lev 3'
  elemdesc(3)='m3/m3'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.0001,0.)  
  elemdesc(1)='wbice4'
  elemdesc(2)='Soil ice lev 4'
  elemdesc(3)='m3/m3'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.0001,0.)  
  elemdesc(1)='wbice5'
  elemdesc(2)='Soil ice lev 5'
  elemdesc(3)='m3/m3'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.0001,0.)  
  elemdesc(1)='wbice6'
  elemdesc(2)='Soil ice lev 6'
  elemdesc(3)='m3/m3'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.0001,0.)  
  elemdesc(1)='snage'
  elemdesc(2)='Snow age'
  elemdesc(3)='none'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.002,0.)
  elemdesc(1)='sflag'
  elemdesc(2)='Snow flag'
  elemdesc(3)='none'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.0004,0.)
end if

if (any(outdata(:,:,38).ne.0.)) then
  elemdesc(1)='tggsn1'
  elemdesc(2)='Snow temperature lev 1'
  elemdesc(3)='K'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
  elemdesc(1)='tggsn2'
  elemdesc(2)='Snow temperature lev 2'
  elemdesc(3)='K'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.) 
  elemdesc(1)='tggsn3'
  elemdesc(2)='Snow temperature lev 3'
  elemdesc(3)='K'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
end if

if (any(outdata(:,:,49).ne.0.)) then
  elemdesc(1)='swater'
  elemdesc(2)='Surface water depth'
  elemdesc(3)='mm'
  Call ncaddvargen(ncidarr,elemdesc,nf_float,3,varid,1.,0.)
end if

if (any(outdata(:,:,50).ne.0.)) then
  elemdesc(1)='ssalin'
  elemdesc(2)='Surface water salinity'
  elemdesc(3)='PSU'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.002,65.)
end if

if (any(outdata(:,:,51).ne.0.)) then
  elemdesc(1)='fracice'
  elemdesc(2)='Sea ice fraction'
  elemdesc(3)='none'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.00005,0.5)
  elemdesc(1)='siced'
  elemdesc(2)='Sea ice depth'
  elemdesc(3)='m'
  Call ncaddvargen(ncidarr,elemdesc,nf_float,3,varid,0.0008,25.)
end if

if (any(outdata(:,:,53).ne.0.)) then
  elemdesc(1)='ocndepth'
  elemdesc(2)='Ocean depth'
  elemdesc(3)='m'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,1.,0.)
  elemdesc(1)='ocheight'
  elemdesc(2)='Ocean surface height'
  elemdesc(3)='m'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.002,0.)
  do i=7,wlev
    write(chr,'(I2.2)') i
    elemdesc(1)='tgg'//chr
    elemdesc(2)='Soil/Ocean temperature lev '//chr
    elemdesc(3)='K'
    Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
  end do
  do i=1,wlev
    write(chr,'(I2.2)') i  
    elemdesc(1)='sal'//chr
    elemdesc(2)='Salinity lev '//chr
    elemdesc(3)='PSU'
    Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,50.)
    elemdesc(1)='uoc'//chr
    elemdesc(2)='x-component current lev '//chr
    elemdesc(3)='m/s'
    Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,0.)
    elemdesc(1)='voc'//chr
    elemdesc(2)='y-component current lev '//chr
    elemdesc(3)='m/s'
    Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,0.)        
  end do
  elemdesc(1)='tggsn4'
  elemdesc(2)='Ice temperature lev 4'
  elemdesc(3)='K'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
  elemdesc(1)='sto'
  elemdesc(2)='Ice storage'
  elemdesc(3)='J'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.1,0.)
  elemdesc(1)='uic'
  elemdesc(2)='x-component ice'
  elemdesc(3)='m/s'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,0.)
  elemdesc(1)='vic'
  elemdesc(2)='y-component ice'
  elemdesc(3)='m/s'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,0.)
  elemdesc(1)='icesal'
  elemdesc(2)='Ice salinity'
  elemdesc(3)='PSU'
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,50.)
end if

ncstatus=nf_enddef(ncidarr(0))

Write(6,*) "Replacing wetfrac1/2/3/4/5/6"
outdatab=max(-2.,min(outdata(:,:,1),2.))
Call ncfindvarid(ncidarr(0),"wetfrac1",outname,varid) ! top
Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
outdatab=max(-2.,min(outdata(:,:,2),2.))
Call ncfindvarid(ncidarr(0),"wetfrac2",outname,varid) ! top
Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
outdatab=max(-2.,min(outdata(:,:,3),2.))
Call ncfindvarid(ncidarr(0),"wetfrac3",outname,varid) ! top
Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
outdatab=max(-2.,min(outdata(:,:,4),2.))
Call ncfindvarid(ncidarr(0),"wetfrac4",outname,varid) ! bottom
Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
outdatab=max(-2.,min(outdata(:,:,5),2.))
Call ncfindvarid(ncidarr(0),"wetfrac5",outname,varid) ! bottom
Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
outdatab=max(-2.,min(outdata(:,:,6),2.))
Call ncfindvarid(ncidarr(0),"wetfrac6",outname,varid) ! bottom
Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)

Write(6,*) "Replacing tgg1/2/3/4/5/6"
outdatab=max(200.,min(outdata(:,:,7),400.))
Call ncfindvarid(ncidarr(0),"tgg1",outname,varid) ! top
Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
outdatab=max(200.,min(outdata(:,:,8),400.))
Call ncfindvarid(ncidarr(0),"tgg2",outname,varid) ! top
Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
outdatab=max(200.,min(outdata(:,:,9),400.))
Call ncfindvarid(ncidarr(0),"tgg3",outname,varid) ! top
Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
outdatab=max(200.,min(outdata(:,:,10),400.))
Call ncfindvarid(ncidarr(0),"tgg4",outname,varid) ! bottom
Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
outdatab=max(200.,min(outdata(:,:,11),400.))
Call ncfindvarid(ncidarr(0),"tgg5",outname,varid) ! bottom
Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
outdatab=max(200.,min(outdata(:,:,12),400.))
Call ncfindvarid(ncidarr(0),"tgg6",outname,varid) ! bottom
Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)

if (any(outdata(:,:,13).ne.0.)) then
  Write(6,*) "Replacing rooftgg1/2/3/4"
  outdatab=max(200.,min(outdata(:,:,13),400.))
  Call ncfindvarid(ncidarr(0),"rooftgg1",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(200.,min(outdata(:,:,14),400.))
  Call ncfindvarid(ncidarr(0),"rooftgg2",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(200.,min(outdata(:,:,15),400.))
  Call ncfindvarid(ncidarr(0),"rooftgg3",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(200.,min(outdata(:,:,16),400.))
  Call ncfindvarid(ncidarr(0),"rooftgg4",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)

  Write(6,*) "Replacing waletgg1/2/3/4"
  outdatab=max(200.,min(outdata(:,:,17),400.))
  Call ncfindvarid(ncidarr(0),"waletgg1",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(200.,min(outdata(:,:,18),400.))
  Call ncfindvarid(ncidarr(0),"waletgg2",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(200.,min(outdata(:,:,19),400.))
  Call ncfindvarid(ncidarr(0),"waletgg3",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(200.,min(outdata(:,:,20),400.))
  Call ncfindvarid(ncidarr(0),"waletgg4",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)

  Write(6,*) "Replacing walwtgg1/2/3/4"
  outdatab=max(200.,min(outdata(:,:,21),400.))
  Call ncfindvarid(ncidarr(0),"walwtgg1",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(200.,min(outdata(:,:,22),400.))
  Call ncfindvarid(ncidarr(0),"walwtgg2",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(200.,min(outdata(:,:,23),400.))
  Call ncfindvarid(ncidarr(0),"walwtgg3",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(200.,min(outdata(:,:,24),400.))
  Call ncfindvarid(ncidarr(0),"walwtgg4",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  
  Write(6,*) "Replacing roadtgg1/2/3/4"
  outdatab=max(200.,min(outdata(:,:,25),400.))
  Call ncfindvarid(ncidarr(0),"roadtgg1",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(200.,min(outdata(:,:,26),400.))
  Call ncfindvarid(ncidarr(0),"roadtgg2",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(200.,min(outdata(:,:,27),400.))
  Call ncfindvarid(ncidarr(0),"roadtgg3",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(200.,min(outdata(:,:,28),400.))
  Call ncfindvarid(ncidarr(0),"roadtgg4",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
end if

if (any(outdata(:,:,29).ne.0.)) then
  Write(6,*) "Replacing urbansm"
  outdatab=max(0.,min(outdata(:,:,29),1.))
  Call ncfindvarid(ncidarr(0),"urbnsmc",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(0.,min(outdata(:,:,30),1.))
  Call ncfindvarid(ncidarr(0),"urbnsmr",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
end if

if (any(outdata(:,:,31).ne.0.)) then
  Write(6,*) "Replacing snd"
  outdatab=max(0.,min(outdata(:,:,31),5.))
  Call ncfindvarid(ncidarr(0),"snd",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
end if

if (any(outdata(:,:,32).ne.0.)) then
  Write(6,*) "Replacing smass1/2/3"
  outdatab=max(0.,min(outdata(:,:,32),400.))
  Call ncfindvarid(ncidarr(0),"smass1",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(0.,min(outdata(:,:,33),400.))
  Call ncfindvarid(ncidarr(0),"smass2",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(0.,min(outdata(:,:,34),400.))
  Call ncfindvarid(ncidarr(0),"smass3",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  
  Write(6,*) "Replacing ssdn1/2/3"
  outdatab=max(0.,min(outdata(:,:,35),400.))
  Call ncfindvarid(ncidarr(0),"ssdn1",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(0.,min(outdata(:,:,36),400.))
  Call ncfindvarid(ncidarr(0),"ssdn2",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(0.,min(outdata(:,:,37),400.))
  Call ncfindvarid(ncidarr(0),"ssdn3",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)

  Write(6,*) "Replacing wbice1/2/3/4/5/6"
  outdatab=max(0.,min(outdata(:,:,41),1.))
  Call ncfindvarid(ncidarr(0),"wbice1",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(0.,min(outdata(:,:,42),1.))
  Call ncfindvarid(ncidarr(0),"wbice2",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(0.,min(outdata(:,:,43),1.))
  Call ncfindvarid(ncidarr(0),"wbice3",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(0.,min(outdata(:,:,44),1.))
  Call ncfindvarid(ncidarr(0),"wbice4",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(0.,min(outdata(:,:,45),1.))
  Call ncfindvarid(ncidarr(0),"wbice5",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(0.,min(outdata(:,:,46),1.))
  Call ncfindvarid(ncidarr(0),"wbice6",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)

  Write(6,*) "Replacing snage"
  outdatab=max(0.,min(outdata(:,:,47),20.))
  Call ncfindvarid(ncidarr(0),"snage",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)

  Write(6,*) "Replacing sflag"
  outdatab=max(0.,min(outdata(:,:,48),4.))
  Call ncfindvarid(ncidarr(0),"sflag",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)      
end if

if (any(outdata(:,:,38).ne.0.)) then
  Write(6,*) "Replacing tggsn1/2/3"
  outdatab=max(100.,min(outdata(:,:,38),400.))
  Call ncfindvarid(ncidarr(0),"tggsn1",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(100.,min(outdata(:,:,39),400.))
  Call ncfindvarid(ncidarr(0),"tggsn2",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(100.,min(outdata(:,:,40),400.))
  Call ncfindvarid(ncidarr(0),"tggsn3",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
end if

if (any(outdata(:,:,49).ne.0.)) then
  Write(6,*) "Replacing swater"
  outdatab=max(0.,min(outdata(:,:,49),6.5E3))
  Call ncfindvarid(ncidarr(0),"swater",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
end if

if (any(outdata(:,:,50).ne.0.)) then
  Write(6,*) "Replacing ssalin"
  outdatab=max(0.,min(outdata(:,:,50),130.))
  Call ncfindvarid(ncidarr(0),"ssalin",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
end if

if (any(outdata(:,:,51).ne.0.)) then
  Write(6,*) "Replacing fracice,siced"
  outdatab=max(0.,min(outdata(:,:,51),1.))
  Call ncfindvarid(ncidarr(0),"fracice",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(0.,min(outdata(:,:,52),50.))
  Call ncfindvarid(ncidarr(0),"siced",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
end if

if (any(outdata(:,:,53).ne.0.)) then
  Write(6,*) "Replacing ocndepth,ocheight,tgg,sal,uoc,voc"
  outdatab=max(0.,min(outdata(:,:,53),10000.))
  Call ncfindvarid(ncidarr(0),"ocndepth",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(-60.,min(outdata(:,:,54),60.))
  Call ncfindvarid(ncidarr(0),"ocheight",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  do i=7,wlev
    write(chr,'(I2.2)') i
    outdatab=max(200.,min(outdata(:,:,54+i),400.))
    Call ncfindvarid(ncidarr(0),"tgg"//chr,outname,varid)
    Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)  
  end do
  do i=1,wlev
    write(chr,'(I2.2)') i
    outdatab=max(0.,min(outdata(:,:,54+wlev+i),100.))
    Call ncfindvarid(ncidarr(0),"sal"//chr,outname,varid)
    Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
    outdatab=max(-100.,min(outdata(:,:,54+2*wlev+i),100.))
    Call ncfindvarid(ncidarr(0),"uoc"//chr,outname,varid)
    Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
    outdatab=max(-100.,min(outdata(:,:,54+3*wlev+i),100.))
    Call ncfindvarid(ncidarr(0),"voc"//chr,outname,varid)
    Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  end do
  Write(6,*) "Replacing tggsn4,sto"
  outdatab=max(100.,min(outdata(:,:,58+4*wlev),400.))
  Call ncfindvarid(ncidarr(0),"tggsn4",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(0.,min(outdata(:,:,60+4*wlev),1000.))
  Call ncfindvarid(ncidarr(0),"sto",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)  
  Write(6,*) "Replacing uic,vic"
  outdatab=max(-100.,min(outdata(:,:,61+4*wlev),100.))
  Call ncfindvarid(ncidarr(0),"uic",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(-100.,min(outdata(:,:,62+4*wlev),100.))
  Call ncfindvarid(ncidarr(0),"vic",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(0.,min(outdata(:,:,63+4*wlev),100.))
  Call ncfindvarid(ncidarr(0),"icesal",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
end if

ncstatus=nf_close(ncidarr(0))

Return
End
