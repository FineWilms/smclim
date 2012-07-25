Program smclim

! This code interpolates a soil moisture climatology onto the Conformal
! cubic grid (i.e., for stand-alone NWP forecasting)

Implicit None

Integer :: nopts
Character*80, dimension(:,:), allocatable :: options

Write(6,*) 'SMCLIM - interpolate/bin climatology soil moisture/temperature on'
Write(6,*) '         conformal-cubic grid (JUL-12)'

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
Allocate(lsdata(1:ccdim(1),1:ccdim(2)),rawdata(1:ccdim(1),1:ccdim(2),1:58+4*wlev))
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
Real, dimension(1:ccdim(1),1:ccdim(2),1:58+4*wlev), intent(in) :: outdata
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

elemdesc=(/ 'wetfrac1', 'Wetness fraction layer 1', 'none' /)
Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.0002,0.)
elemdesc=(/ 'wetfrac2', 'Wetness fraction layer 2', 'none' /)
Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.0002,0.)
elemdesc=(/ 'wetfrac3', 'Wetness fraction layer 3', 'none' /)
Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.0002,0.)
elemdesc=(/ 'wetfrac4', 'Wetness fraction layer 4', 'none' /)
Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.0002,0.)
elemdesc=(/ 'wetfrac5', 'Wetness fraction layer 5', 'none' /)
Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.0002,0.)
elemdesc=(/ 'wetfrac6', 'Wetness fraction layer 6', 'none' /)
Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.0002,0.)

elemdesc=(/ 'tgg1', 'Soil temperature lev 1', 'K' /)
Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
elemdesc=(/ 'tgg2', 'Soil temperature lev 2', 'K' /)
Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
elemdesc=(/ 'tgg3', 'Soil temperature lev 3', 'K' /)
Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
elemdesc=(/ 'tgg4', 'Soil temperature lev 4', 'K' /)
Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
elemdesc=(/ 'tgg5', 'Soil temperature lev 5', 'K' /)
Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
elemdesc=(/ 'tgg6', 'Soil temperature lev 6', 'K' /)
Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)

if (any(outdata(:,:,13).ne.0.)) then
  elemdesc=(/ 'rooftgg1', 'roof temperature lev 1', 'K' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
  elemdesc=(/ 'rooftgg2', 'roof temperature lev 2', 'K' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
  elemdesc=(/ 'rooftgg3', 'roof temperature lev 3', 'K' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
  elemdesc=(/ 'waletgg1', 'wall temperature lev 1', 'K' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
  elemdesc=(/ 'waletgg2', 'wall temperature lev 2', 'K' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
  elemdesc=(/ 'waletgg3', 'wall temperature lev 3', 'K' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
  elemdesc=(/ 'walwtgg1', 'wall temperature lev 1', 'K' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
  elemdesc=(/ 'walwtgg2', 'wall temperature lev 2', 'K' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
  elemdesc=(/ 'walwtgg3', 'wall temperature lev 3', 'K' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
  elemdesc=(/ 'roadtgg1', 'road temperature lev 1', 'K' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
  elemdesc=(/ 'roadtgg2', 'road temperature lev 2', 'K' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
  elemdesc=(/ 'roadtgg3', 'road temperature lev 3', 'K' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
end if

if (any(outdata(:,:,25).ne.0.)) then
  elemdesc=(/ 'urbnsmc', 'urban canyon soil mositure', 'm3/m3' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.0001,0.5)
  elemdesc=(/ 'urbnsmr', 'urban roof soil mositure', 'm3/m3' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.0001,0.5)
end if

if (any(outdata(:,:,27).ne.0.)) then
  elemdesc=(/ 'snd', 'Snow depth (liquid water)', 'mm' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_float,3,varid,1.,0.)
end if

if (any(outdata(:,:,28).ne.0.)) then
  elemdesc=(/ 'smass1', 'Snow mass lev 1', 'K' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.01,250.)
  elemdesc=(/ 'smass2', 'Snow mass lev 2', 'K' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.01,250.)
  elemdesc=(/ 'smass3', 'Snow mass lev 3', 'K' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.01,250.)
  elemdesc=(/ 'ssdn1', 'Snow density lev 1', 'K' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.01,250.)
  elemdesc=(/ 'ssdn2', 'Snow density lev 2', 'K' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.01,250.)
  elemdesc=(/ 'ssdn3', 'Snow density lev 3', 'K' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.01,250.)
  elemdesc=(/ 'wbice1', 'Soil ice lev 1', 'm3/m3' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.0001,0.)  
  elemdesc=(/ 'wbice2', 'Soil ice lev 2', 'm3/m3' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.0001,0.)  
  elemdesc=(/ 'wbice3', 'Soil ice lev 3', 'm3/m3' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.0001,0.)  
  elemdesc=(/ 'wbice4', 'Soil ice lev 4', 'm3/m3' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.0001,0.)  
  elemdesc=(/ 'wbice5', 'Soil ice lev 5', 'm3/m3' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.0001,0.)  
  elemdesc=(/ 'wbice6', 'Soil ice lev 6', 'm3/m3' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.0001,0.)  
  elemdesc=(/ 'snage', 'Snow age', 'none' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.002,0.)
  elemdesc=(/ 'sflag', 'Snow flag', 'none' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.0004,0.)
end if

if (any(outdata(:,:,34).ne.0.)) then
  elemdesc=(/ 'tggsn1', 'Snow temperature lev 1', 'K' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
  elemdesc=(/ 'tggsn2', 'Snow temperature lev 2', 'K' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.) 
  elemdesc=(/ 'tggsn3', 'Snow temperature lev 3', 'K' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
end if

if (any(outdata(:,:,45).ne.0.)) then
  elemdesc=(/ 'swater', 'Surface water', 'mm' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.002,65.)
end if

if (any(outdata(:,:,46).ne.0.)) then
  elemdesc=(/ 'fracice', 'Sea ice fraction', 'none' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.00005,0.5)
  elemdesc=(/ 'siced', 'Sea ice depth', 'm' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_float,3,varid,0.0008,25.)
end if

if (any(outdata(:,:,48).ne.0.)) then
  elemdesc=(/ 'ocndepth', 'Ocean depth', 'm' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,1.,0.)
  elemdesc=(/ 'ocheight', 'Ocean surface height', 'm' /)
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
  elemdesc=(/ 'tggsn4', 'Ice temperature lev 4', 'K' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,250.)
  elemdesc=(/ 'sto', 'Ice storage', 'J' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.1,0.)
  elemdesc=(/ 'uic', 'x-component ice', 'm/s' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,0.)
  elemdesc=(/ 'vic', 'y-component ice', 'm/s' /)
  Call ncaddvargen(ncidarr,elemdesc,nf_short,3,varid,0.005,0.)
  elemdesc=(/ 'icesal', 'Ice salinity', 'PSU' /)
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
  Write(6,*) "Replacing rooftgg1/2/3"
  outdatab=max(200.,min(outdata(:,:,13),400.))
  Call ncfindvarid(ncidarr(0),"rooftgg1",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(200.,min(outdata(:,:,14),400.))
  Call ncfindvarid(ncidarr(0),"rooftgg2",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(200.,min(outdata(:,:,15),400.))
  Call ncfindvarid(ncidarr(0),"rooftgg3",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)

  Write(6,*) "Replacing waletgg1/2/3"
  outdatab=max(200.,min(outdata(:,:,16),400.))
  Call ncfindvarid(ncidarr(0),"waletgg1",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(200.,min(outdata(:,:,17),400.))
  Call ncfindvarid(ncidarr(0),"waletgg2",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(200.,min(outdata(:,:,18),400.))
  Call ncfindvarid(ncidarr(0),"waletgg3",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)

  Write(6,*) "Replacing walwtgg1/2/3"
  outdatab=max(200.,min(outdata(:,:,19),400.))
  Call ncfindvarid(ncidarr(0),"walwtgg1",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(200.,min(outdata(:,:,20),400.))
  Call ncfindvarid(ncidarr(0),"walwtgg2",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(200.,min(outdata(:,:,21),400.))
  Call ncfindvarid(ncidarr(0),"walwtgg3",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)

  Write(6,*) "Replacing roadtgg1/2/3"
  outdatab=max(200.,min(outdata(:,:,22),400.))
  Call ncfindvarid(ncidarr(0),"roadtgg1",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(200.,min(outdata(:,:,23),400.))
  Call ncfindvarid(ncidarr(0),"roadtgg2",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(200.,min(outdata(:,:,24),400.))
  Call ncfindvarid(ncidarr(0),"roadtgg3",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
end if

if (any(outdata(:,:,25).ne.0.)) then
  Write(6,*) "Replacing urbansm"
  outdatab=max(0.,min(outdata(:,:,25),1.))
  Call ncfindvarid(ncidarr(0),"urbnsmc",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(0.,min(outdata(:,:,26),1.))
  Call ncfindvarid(ncidarr(0),"urbnsmr",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
end if

if (any(outdata(:,:,27).ne.0.)) then
  Write(6,*) "Replacing snd"
  outdatab=max(0.,min(outdata(:,:,27),5.))
  Call ncfindvarid(ncidarr(0),"snd",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
end if

if (any(outdata(:,:,28).ne.0.)) then
  Write(6,*) "Replacing smass1/2/3"
  outdatab=max(0.,min(outdata(:,:,28),400.))
  Call ncfindvarid(ncidarr(0),"smass1",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(0.,min(outdata(:,:,29),400.))
  Call ncfindvarid(ncidarr(0),"smass2",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(0.,min(outdata(:,:,30),400.))
  Call ncfindvarid(ncidarr(0),"smass3",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  
  Write(6,*) "Replacing ssdn1/2/3"
  outdatab=max(0.,min(outdata(:,:,31),400.))
  Call ncfindvarid(ncidarr(0),"ssdn1",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(0.,min(outdata(:,:,32),400.))
  Call ncfindvarid(ncidarr(0),"ssdn2",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(0.,min(outdata(:,:,33),400.))
  Call ncfindvarid(ncidarr(0),"ssdn3",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)

  Write(6,*) "Replacing wbice1/2/3/4/5/6"
  outdatab=max(0.,min(outdata(:,:,37),1.))
  Call ncfindvarid(ncidarr(0),"wbice1",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(0.,min(outdata(:,:,38),1.))
  Call ncfindvarid(ncidarr(0),"wbice2",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(0.,min(outdata(:,:,39),1.))
  Call ncfindvarid(ncidarr(0),"wbice3",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(0.,min(outdata(:,:,40),1.))
  Call ncfindvarid(ncidarr(0),"wbice4",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(0.,min(outdata(:,:,41),1.))
  Call ncfindvarid(ncidarr(0),"wbice5",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(0.,min(outdata(:,:,42),1.))
  Call ncfindvarid(ncidarr(0),"wbice6",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)

  Write(6,*) "Replacing snage"
  outdatab=max(0.,min(outdata(:,:,43),20.))
  Call ncfindvarid(ncidarr(0),"snage",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)

  Write(6,*) "Replacing sflag"
  outdatab=max(0.,min(outdata(:,:,44),4.))
  Call ncfindvarid(ncidarr(0),"sflag",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)      
end if

if (any(outdata(:,:,34).ne.0.)) then
  Write(6,*) "Replacing tggsn1/2/3"
  outdatab=max(100.,min(outdata(:,:,34),400.))
  Call ncfindvarid(ncidarr(0),"tggsn1",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(100.,min(outdata(:,:,35),400.))
  Call ncfindvarid(ncidarr(0),"tggsn2",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(100.,min(outdata(:,:,36),400.))
  Call ncfindvarid(ncidarr(0),"tggsn3",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
end if

if (any(outdata(:,:,45).ne.0.)) then
  Write(6,*) "Replacing swater"
  outdatab=max(0.,min(outdata(:,:,45),130.))
  Call ncfindvarid(ncidarr(0),"swater",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
end if

if (any(outdata(:,:,46).ne.0.)) then
  Write(6,*) "Replacing fracice,siced"
  outdatab=max(0.,min(outdata(:,:,46),1.))
  Call ncfindvarid(ncidarr(0),"fracice",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(0.,min(outdata(:,:,47),50.))
  Call ncfindvarid(ncidarr(0),"siced",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
end if

if (any(outdata(:,:,48).ne.0.)) then
  Write(6,*) "Replacing ocndepth,ocheight,tgg,sal,uoc,voc"
  outdatab=max(0.,min(outdata(:,:,48),10000.))
  Call ncfindvarid(ncidarr(0),"ocndepth",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(-60.,min(outdata(:,:,49),60.))
  Call ncfindvarid(ncidarr(0),"ocheight",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  do i=7,wlev
    write(chr,'(I2.2)') i
    outdatab=max(200.,min(outdata(:,:,49+i),400.))
    Call ncfindvarid(ncidarr(0),"tgg"//chr,outname,varid)
    Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)  
  end do
  do i=1,wlev
    write(chr,'(I2.2)') i
    outdatab=max(0.,min(outdata(:,:,49+wlev+i),100.))
    Call ncfindvarid(ncidarr(0),"sal"//chr,outname,varid)
    Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
    outdatab=max(-100.,min(outdata(:,:,49+2*wlev+i),100.))
    Call ncfindvarid(ncidarr(0),"uoc"//chr,outname,varid)
    Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
    outdatab=max(-100.,min(outdata(:,:,49+3*wlev+i),100.))
    Call ncfindvarid(ncidarr(0),"voc"//chr,outname,varid)
    Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  end do
  Write(6,*) "Replacing tggsn4,sto"
  outdatab=max(100.,min(outdata(:,:,53+4*wlev),400.))
  Call ncfindvarid(ncidarr(0),"tggsn4",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(0.,min(outdata(:,:,55+4*wlev),1000.))
  Call ncfindvarid(ncidarr(0),"sto",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)  
  Write(6,*) "Replacing uic,vic"
  outdatab=max(-100.,min(outdata(:,:,56+4*wlev),100.))
  Call ncfindvarid(ncidarr(0),"uic",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(-100.,min(outdata(:,:,57+4*wlev),100.))
  Call ncfindvarid(ncidarr(0),"vic",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
  outdatab=max(0.,min(outdata(:,:,58+4*wlev),100.))
  Call ncfindvarid(ncidarr(0),"icesal",outname,varid)
  Call ncwritedatgen(ncidarr,outdatab,dimnum,varid)
end if

ncstatus=nf_close(ncidarr(0))

Return
End
