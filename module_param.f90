module variables
! Boundary conditions : ncl = 2 --> Dirichlet
! Boundary conditions : ncl = 1 --> Free-slip
! Boundary conditions : ncl = 0 --> Periodic
! l: power of 2,3,4,5 and 6
! if ncl = 1 or 2, --> n  = 2l+ 1 
!                  --> nm = n - 1 
!                  --> m  = n + 1
! If ncl = 0,      --> n  = 2*l
!                  --> nm = n  
!                  --> m  = n + 2
integer,parameter :: nx=257,ny=128,nz=1
integer,parameter :: nxm=nx-1,nym=ny,nzm=nz 
integer,parameter :: mx=nx+1,my=ny+2,mz=nz
integer,parameter :: nwork=512*max(mx,my,mz)
integer,parameter :: ntrigsX=4*nx,ntrigsY=4*ny,ntrigsZ=4*nz,nfax=19
integer,parameter :: nphi=1
!end module variables

integer,parameter :: nb_x=nxm/20, nb_y=nym/12
integer,parameter :: nb_front=2*nxm/20+2*nym/12
integer,parameter :: nb_src=2*(2*nxm/20+2*nym/12)/3+16

!module géométrie complexe
!new_IBM : 
integer,parameter               :: nraf=10,nobjmax=20
integer,parameter               :: nxraf=(nx-1)*nraf+1,nyraf=(ny-1)*nraf+1,nzraf=(nz-1)*nraf+1
integer                         :: ilag,npif,izap
integer,dimension(ny,nz)        :: nobjx
integer,dimension(nx,nz)        :: nobjy
integer,dimension(nx,ny)        :: nobjz
integer,dimension(nym,nz )      :: nobjxdy
integer,dimension(nxm,nz )      :: nobjydx
integer,dimension(nxm,ny )      :: nobjzdx
integer,dimension(nx ,nym)      :: nobjzdy
integer,dimension(nym,nzm)      :: nobjxdyz
integer,dimension(nxm,nzm)      :: nobjydxz
integer,dimension(nxm,nym)      :: nobjzdxy
integer,dimension(0:nobjmax,ny ,nz ) :: nxipif,nxfpif
integer,dimension(0:nobjmax,nx ,nz ) :: nyipif,nyfpif
integer,dimension(0:nobjmax,nx ,ny ) :: nzipif,nzfpif
real(8),dimension(  nobjmax,ny ,nz ) :: xi,xf
real(8),dimension(  nobjmax,nx ,nz ) :: yi,yf
real(8),dimension(  nobjmax,nx ,ny ) :: zi,zf
real(8),dimension(  nobjmax,nym,nz ) :: xidecy ,xfdecy
real(8),dimension(  nobjmax,nxm,nz ) :: yidecx ,yfdecx
real(8),dimension(  nobjmax,nxm,ny ) :: zidecx ,zfdecx
real(8),dimension(  nobjmax,nx ,nym) :: zidecy ,zfdecy
real(8),dimension(  nobjmax,nxm,nym) :: xidecyz,xfdecyz
real(8),dimension(  nobjmax,nxm,nzm) :: yidecxz,yfdecxz
real(8),dimension(  nobjmax,nxm,nym) :: zidecxy,zfdecxy

!module filtre
real(8),dimension(nx) :: fifx,ficx,fibx,fiffx,fibbx,fiz1x,fiz2x
real(8),dimension(nx,2) ::filax,filaxp
real(8),dimension(nx) :: fifxp,ficxp,fibxp,fiffxp,fibbxp
real(8),dimension(ny) :: fify,ficy,fiby,fiffy,fibby,fiz1y,fiz2y
real(8),dimension(ny,2) ::filay,filayp
real(8),dimension(ny) :: fifyp,ficyp,fibyp,fiffyp,fibbyp
real(8),dimension(nz) :: fifz,ficz,fibz,fiffz,fibbz,fiz1z,fiz2z
real(8),dimension(nz,2) ::filaz,filazp
real(8),dimension(nz) :: fifzp,ficzp,fibzp,fiffzp,fibbzp

!module derive
real(8),dimension(nx) :: ffx,fcx,fbx,sfx,scx,sbx,fsx,fwx,ssx,swx
real(8),dimension(nx) :: ffxp,fsxp,fwxp,sfxp,ssxp,swxp
real(8),dimension(ny) :: ffy,fcy,fby,sfy,scy,sby,fsy,fwy,ssy,swy
real(8),dimension(ny) :: ffyp,fsyp,fwyp,sfyp,ssyp,swyp
real(8),dimension(nz) :: ffz,fcz,fbz,sfz,scz,sbz,fsz,fwz,ssz,swz
real(8),dimension(nz) :: ffzp,fszp,fwzp,sfzp,sszp,swzp
real(8),dimension(ny,nz) :: sx,vx 
real(8),dimension(nx,nz) :: sy,vy
real(8),dimension(nx,ny) :: sz,vz

!module pression
real(8),dimension(ny,nz) :: dpdyx1,dpdyxn,dpdzx1,dpdzxn,dpdyx12,dpdyxn2,dpdzx12,dpdzxn2
real(8),dimension(nx,nz) :: dpdxy1,dpdxyn,dpdzy1,dpdzyn,dpdxy12,dpdxyn2,dpdzy12,dpdzyn2
real(8),dimension(nx,ny) :: dpdxz1,dpdxzn,dpdyz1,dpdyzn,dpdxz12,dpdxzn2,dpdyz12,dpdyzn2

!module inflow
real(8),dimension(ny,nz) :: bxx1,bxy1,bxz1,bxxn,bxyn,bxzn,bxo,byo,bzo
real(8),dimension(nx,nz) :: byx1,byy1,byz1,byxn,byyn,byzn
real(8),dimension(nx,ny) :: bzx1,bzy1,bzz1,bzxn,bzyn,bzzn

!module derpres
real(8),dimension(nxm) :: cfx6,ccx6,cbx6,cfxp6,ciwxp6,csxp6,cwxp6,csx6,cwx6,cifx6,cicx6,cisx6   
real(8),dimension(nxm) :: cibx6,cifxp6,cisxp6,ciwx6
real(8),dimension(nx) :: cfi6,cci6,cbi6,cfip6,csip6,cwip6,csi6,cwi6,cifi6,cici6,cibi6,cifip6  
real(8),dimension(nx) :: cisip6,ciwip6,cisi6,ciwi6 
real(8),dimension(nym) :: cfy6,ccy6,cby6,cfyp6,csyp6,cwyp6,csy6 
real(8),dimension(nym) :: cwy6,cify6,cicy6,ciby6,cifyp6,cisyp6,ciwyp6,cisy6,ciwy6 
real(8),dimension(ny) :: cfi6y,cci6y,cbi6y,cfip6y,csip6y,cwip6y,csi6y,cwi6y,cifi6y,cici6y  
real(8),dimension(ny) :: cibi6y,cifip6y,cisip6y,ciwip6y,cisi6y,ciwi6y  
real(8),dimension(nzm) :: cfz6,ccz6,cbz6,cfzp6,cszp6,cwzp6,csz6 
real(8),dimension(nzm) :: cwz6,cifz6,cicz6,cibz6,cifzp6,ciszp6,ciwzp6,cisz6,ciwz6 
real(8),dimension(nz) :: cfi6z,cci6z,cbi6z,cfip6z,csip6z,cwip6z,csi6z,cwi6z,cifi6z,cici6z  
real(8),dimension(nz) :: cibi6z,cifip6z,cisip6z,ciwip6z,cisi6z,ciwi6z 

!module waves
real(8),dimension(mx) :: xkx,xk2,exs
real(8),dimension(my) :: yky,yk2,eys 
real(8),dimension(mz) :: zkz,zk2,ezs

!module mesh
real(8),dimension(ny) :: ppy,pp2y,pp4y
real(8),dimension(ny) :: ppyi,pp2yi,pp4yi
real(8),dimension(ny) :: yp,ypi

! Variables for reading initial conditions of reference case
real(8),dimension(ny) :: uxin,uxout,uyin,uyout
real(8),dimension(nx) :: uytop,uybot,uxtop,uxbot

real(8),dimension(nx,ny)   :: ux_ana, uy_ana, pp_ana
real(8),dimension(nxm,nym) :: p_ana

! Phase averaged variables
real(8),dimension(nx,ny,nz)   ::  txx02d, txx02u, txxbm, txxbm2d, txxbm2u, txxtm, txxtm2d, txxtm2u
real(8),dimension(nx,ny,nz)   ::  tzz02d, tzz02u, tzzbm, tzzbm2d, tzzbm2u, tzztm, tzztm2d, tzztm2u
real(8),dimension(nx,ny,nz)   ::  ux02d, ux02u, uxbm, uxbm2d, uxbm2u, uxtm, uxtm2d, uxtm2u
real(8),dimension(nx,ny,nz)   ::  uy02d, uy02u, uybm, uybm2d, uybm2u, uytm, uytm2d, uytm2u

real(8),dimension(nx,ny,nz)   ::  uxx02d, uxx02u, uxxbm, uxxbm2d, uxxbm2u, uxxtm, uxxtm2d, uxxtm2u
real(8),dimension(nx,ny,nz)   ::  uxy02d, uxy02u, uxybm, uxybm2d, uxybm2u, uxytm, uxytm2d, uxytm2u
real(8),dimension(nx,ny,nz)   ::  uyy02d, uyy02u, uyybm, uyybm2d, uyybm2u, uyytm, uyytm2d, uyytm2u
real(8),dimension(nx,ny,nz)   ::  txx_ave, tzz_ave, ux_ave, uy_ave, uxx_ave, uyy_ave, uxy_ave, epsi_ave 

end module variables

module param
  integer :: nclx,ncly,nclz
  integer :: ifft,ifiltre, ivirtuel,istret,iforc_entree,iturb, ipen
  integer :: iecoule, iskew, ientree, nschema, idebut, ifin, iles
  integer :: isave,ilit,idebmod, imodulo, idemarre, icommence, irecord
  integer :: iscalaire
  integer :: iord
  integer :: nxboite, istat,iread,iavance_temps 
  real(8) :: xlx,yly,zlz,dx,dy,dz,dx2,dy2,dz2
  real(8) :: dt,xnu,bruit,pi,twopi,u1,u2,beta,sc,alpha
  real(8) :: t
  integer :: itr,itime
  character :: filecharge*80, filesauve*80, filebruit*80, nchamp*80,filepath*80, fileturb*80 
  real(8),dimension(5) :: adt,bdt,gdt,cdt,ddt
end module param

module magnets
real(8), dimension(100) :: xMag,yMag,zMag,angle,Br,a,b,thickness
real(8) :: jX,jY,jZ
real(8) :: timeIni,duration,nextTime
integer :: nBlocks,nMagnets,mem1
real(8) :: t1,t2
real(8) :: rho,uref,lref
end module magnets

module IBM
  real(8) :: cex,cey,cez,ra,alfa,c,re,rads,ubcx,ubcy
  integer :: i5,jj,j1,j2,i55
  integer,parameter :: lskin=500
  real(8) :: AoA,cx,cy,utrans,vtrans,omega,dAoAdt
  real(8),dimension(lskin,2) :: body,urot
  real(8),dimension(lskin,6) :: vertex
  real(8),dimension(lskin,2) :: cdcl
  real(8) :: amax, freq, length, thick, aspec, aspec2, distance
  real(8) :: amax2, freq2
  integer :: flap, angexc, geoah, flapex
end module IBM

module derivX
  real(8) :: alcaix6,acix6,bcix6
  real(8) :: ailcaix6, aicix6, bicix6,cicix6 
  real(8) :: alfa1x,af1x,bf1x,cf1x,df1x,alfa2x,af2x,alfanx,afnx,bfnx
  real(8) :: cfnx,dfnx,alfamx,afmx,alfaix,afix,bfix,alsa1x,as1x,bs1x
  real(8) :: cs1x,ds1x,alsa2x,as2x,alsanx,asnx,bsnx,csnx,dsnx,alsamx
  real(8) :: asmx,alsaix,asix,bsix,csix,as3x,bs3x,astx,bstx
end module derivX

module derivY 
  real(8) :: alcaiy6,aciy6,bciy6
  real(8) :: ailcaiy6, aiciy6, biciy6,ciciy6 
  real(8) :: alfa1y,af1y,bf1y,cf1y,df1y,alfa2y,af2y,alfany,afny,bfny
  real(8) :: cfny,dfny,alfamy,afmy,alfajy,afjy,bfjy,alsa1y,as1y,bs1y
  real(8) :: cs1y,ds1y,alsa2y,as2y,alsany,asny,bsny,csny,dsny,alsamy
  real(8) :: asmy,alsajy,asjy,bsjy,csjy,as3y,bs3y,asty,bsty 
end module derivY

module derivZ
  real(8) :: alcaiz6,aciz6,bciz6 
  real(8) :: ailcaiz6, aiciz6, biciz6,ciciz6 
  real(8) :: alfa1z,af1z,bf1z,cf1z,df1z,alfa2z,af2z,alfanz,afnz,bfnz
  real(8) :: cfnz,dfnz,alfamz,afmz,alfakz,afkz,bfkz,alsa1z,as1z,bs1z
  real(8) :: cs1z,ds1z,alsa2z,as2z,alsanz,asnz,bsnz,csnz,dsnz,alsamz
  real(8) :: asmz,alsakz,askz,bskz,cskz,as3z,bs3z,astz,bstz 
end module derivZ

module parfiX
  real(8) :: fia1x, fib1x, fic1x, fid1x, fie1x, fia2x, fib2x, fic2x, fid2x
  real(8) :: fie2x, fia3x, fib3x, fic3x, fid3x, fie3x, fianx, fibnx, ficnx, fidnx
  real(8) :: fienx, fiamx, fibmx, ficmx, fidmx, fiemx, fiapx, fibpx, ficpx, fidpx
  real(8) :: fiepx, fiaix, fibix, ficix, fidix, fialx, fibex, fih1x, fih2x, fih3x,fih4x 
end module parfiX
!
module parfiY
  real(8) :: fia1y, fib1y, fic1y, fid1y, fie1y, fia2y, fib2y, fic2y, fid2y
  real(8) :: fie2y, fia3y, fib3y, fic3y, fid3y, fie3y, fiany, fibny, ficny, fidny
  real(8) :: fieny, fiamy, fibmy, ficmy, fidmy, fiemy, fiapy, fibpy, ficpy, fidpy
  real(8) :: fiepy, fiaiy, fibiy, ficiy, fidiy, fialy, fibey, fih1y, fih2y, fih3y,fih4y 
end module parfiY

module parfiZ
  real(8) :: fia1z, fib1z, fic1z, fid1z, fie1z, fia2z, fib2z, fic2z, fid2z
  real(8) :: fie2z, fia3z, fib3z, fic3z, fid3z, fie3z, fianz, fibnz, ficnz, fidnz
  real(8) :: fienz, fiamz, fibmz, ficmz, fidmz, fiemz, fiapz, fibpz, ficpz, fidpz
  real(8) :: fiepz, fiaiz, fibiz, ficiz, fidiz, fialz, fibez, fih1z, fih2z, fih3z,fih4z 
end module parfiZ

module aeroforce
  real(4) :: xforce,yforce,zforce,xmom,ymom,zmom
  real(4) :: f2x,f2y,f2z,xfront,f2xx
  real(4) :: xld,xud,yld,yud,zld,zud
  integer :: ild,iud,jld,jud,kld,kud,iaero
  integer :: counter
end module aeroforce


