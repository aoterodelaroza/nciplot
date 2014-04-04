! Copyright (c) 2013 Alberto Otero de la Roza <aoterodelaroza@ucmerced.edu>,
! Julia Conteras-Garcia <julia.contreras.garcia@gmail.com>, 
! Erin R. Johnson <ejohnson29@ucmerced.edu>, and Weitao Yang
! <weitao.yang@duke.edu>
!
! nciplot is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
program nciplot
  use param
  use tools_io
  use tools_math
  use reader
  use props
  implicit none

  integer, parameter :: mfiles = 100 !< only applies to rhom().

  character*(mline) :: argv(2), oname
  integer :: argc, nfiles, ifile, idx, istat, ntotal
  integer :: i, j, k, nn0, nnf, lp, n1
  character*(mline) :: filein, line, oline, word, wx, wc
  logical :: ok, ispromol
  real*8 :: rdum
  ! the molecular info
  type(molecule), allocatable :: m(:)
  ! logical units
  integer :: lugc, ludc, luvmd, ludat, luelf, luxc, luchk
  logical :: lchk
  ! cubes
  real*8, allocatable, dimension(:,:,:) :: crho, cgrad, celf, cxc
  ! ligand, intermolecular keywords
  logical :: ligand, inter, intra
  real*8 :: rthres
  integer :: udat0
  ! radius and cube keywords
  logical :: autor
  real*8 :: x(3), xinit(3), xmax(3), xinc(3)
  integer :: nstep(3)
  ! noutput keyword
  integer :: noutput
  ! cutoffs
  real*8 :: rhocut, dimcut
  ! cutplot
  real*8 :: rhoplot, isordg
  ! discarding rho parameter
  real*8 :: rhoparam, rhoparam2
  ! properties of rho
  real*8 :: rho, grad(3), dimgrad, grad2, hess(3,3), elf, exc
  integer, parameter :: mfrag = 100
  real*8 :: rhom(mfrag)
  ! eispack
  real*8 :: wk1(3), wk2(3), heigs(3), hvecs(3,3)
  ! elf
  logical :: doelf
  ! xc
  integer :: ixc(2)
  ! fragments
  integer :: nfrag
  logical :: autofrag
  ! chk file
  logical :: alcrho, alcgrad, alcelf, alcxc
  real*8 :: xinit0(3), xinc0(3)
  integer :: nstep0(3)

  ! initialize
  call param_init()

  ! I/O units, process arguments
  call getargs(argc,argv)
  if (argc >= 1) then
     open (uin,file=argv(1),status='old')
     if (argc == 2) then
        open(uout,file=argv(2),status='unknown')
     endif
  endif

  ! header
  call header()
  call tictac(' # Start')

  ! read files, define fragments
  read (uin,*) nfiles
  if (nfiles > mfiles) call error('nciplot','too many files, increase mfiles',faterr)
  allocate(m(nfiles),stat=istat)
  if (istat /= 0) call error('nciplot','could not allocate memory for molecules',faterr)
  ntotal = 0
  do ifile = 1, nfiles
     read (uin,'(a)') filein
     filein = trim(adjustl(filein))
     inquire(file=filein,exist=ok)
     if(.not.ok) &
        call error('nciplot','requested file does not exist: '//trim(filein),faterr)
     m(ifile) = readfile(filein)
     if (ifile == 1) oname = filein(1:index(filein,'.',.true.)-1)
     ntotal = ntotal + m(ifile)%n
     ! define fragments
     do i = 1, m(ifile)%n
        m(ifile)%ifrag(i) = ifile
     end do
  enddo
  nfrag = nfiles
  autofrag = .true.
  if (nfrag > mfrag) then
     call error('nciplot','too many fragments. Increase mfrag',faterr)
  end if

  if (any(m(:)%ifile == ifile_wfn) .and..not.all(m(:)%ifile == ifile_wfn)) then
     call error('nciplot','mixing xyz and wfn not allowed',faterr)
  end if
  ispromol = .not.all(m(:)%ifile == ifile_wfn)

  ! read density grids (props)
  call init_rhogrid(m,nfiles)

  ! by default, use density grids for heavier or charged atoms.
  do i = 1, nfiles
     if (m(i)%ifile == ifile_xyz .and. (any(m(i)%z > atomic_zmax) .or. any(m(i)%q > 0))) then
        m(i)%ifile = ifile_grd
        do j = 1, m(i)%n
           if (m(i)%z(j) > atomic_zmax .and..not.grd(iztype(m(i)%z(j),m(i)%q(j)))%init) then 
              call error('nciplot','Some atomic density grids for heavy atoms are needed but not initialized',faterr)
           end if
        end do
     end if
  end do

  ! default values
  rhocut = 0.2d0
  xinc = 0.1d0
  if (any(m(:)%ifile == ifile_wfn)) then
     dimcut = 2.0d0
     isordg=0.5d0
     rhoplot=0.05d0
  else
     dimcut = 1.0d0
     isordg=0.3d0
     rhoplot=0.07d0
  end if
  rhoparam = 0.95d0
  rhoparam2 = 0.75d0
  noutput = 3
  udat0 = 1
  autor = .true.
  ligand = .false.
  inter = .false.
  rthres = 2.d0
  doelf = .false.
  ixc = 0
  ! read the rest of (optional) keywords
  xinit = m(1)%x(:,1)
  xmax = m(1)%x(:,1)
  do i = 1, nfiles
     do j = 1, m(i)%n
        xinit = min(xinit,m(i)%x(:,j))
        xmax = max(xmax,m(i)%x(:,j))
     end do
  end do
  do while (.true.)
     read (uin,'(a)',end=11) line
     line = trim(adjustl(line))
     oline = line
     call upper(line)
     if (line(1:1) == "#") cycle ! skip comments
     if (len(trim(line)) < 1) cycle ! skip blank lines

     idx = index(line,' ')
     word = line(1:idx-1)
     line = line(idx:)
     oline = oline(idx:)
     select case(trim(word))
     case ("RTHRES")
        read (line,*) rthres
        rthres = rthres / bohrtoa
     case ("LIGAND")
        ligand = .true.
        inter = .true.
        read (line,*) udat0, rthres
        rthres = rthres / bohrtoa
     case ("INTERMOLECULAR")
        inter = .true.
     case ("RADIUS")
        autor=.false.
        read(line,*) x, rdum
        xinit = (x - rdum) / bohrtoa
        xmax = (x + rdum) / bohrtoa
     case ("ONAME")
        read(oline,'(a)') oname
        oname = trim(adjustl(oname))
        oname = oname(1:index(oname,' '))
     case ("OUTPUT")
        read(line,*) noutput
     case ("CUBE")
        autor = .false.
        read(line,*) xinit, xmax
        xinit = xinit / bohrtoa
        xmax = xmax / bohrtoa
     case ("ATCUBE")
        autor = .false.
        xinit = 1d40
        xmax = -1d40
        do while (.true.)
           read (uin,'(a)') line
           line = trim(adjustl(line))
           call upper(line)
           if (line(1:1) == "#") cycle ! skip comments
           if (len(trim(line)) < 1) cycle ! skip blank lines
           lp = 1
           idx = index(line,' ')
           word = line(1:idx-1)
           if (trim(word) /= "END" .and. trim(word) /= "ENDATCUBE") then
              ok = isinteger(ifile,line,lp)
              if (.not.ok) call error('nciplot','bad atcube syntax',faterr)
              if (ifile < 1.or. ifile > nfiles) call error('nciplot','atcube: wrong file number',faterr)
              ok = isinteger(n1,line,lp)
              do while (ok)
                 if (n1<1.or.n1>m(ifile)%n) call error('nciplot','atcube: wrong atom number',faterr)
                 xinit = min(xinit,m(ifile)%x(:,n1))
                 xmax = max(xmax,m(ifile)%x(:,n1))
                 ok = isinteger(n1,line,lp)
              end do
           else
              exit
           end if
        end do
        xinit = xinit - rthres
        xmax = xmax + rthres
     case ("FRAGMENT")
        if (autofrag) then
           nfrag = 0
           do ifile = 1, nfiles
              do i = 1, m(ifile)%n
                 m(ifile)%ifrag(i) = 0
              end do
           end do
        end if
        autofrag = .false.
        inter = .true.

        nfrag = nfrag + 1
        do while (.true.)
           read (uin,'(a)') line
           line = trim(adjustl(line))
           call upper(line)
           if (line(1:1) == "#") cycle ! skip comments
           if (len(trim(line)) < 1) cycle ! skip blank lines
           lp = 1
           idx = index(line,' ')
           word = line(1:idx-1)
           if (trim(word) /= "END" .and. trim(word) /= "ENDFRAGMENT") then
              ok = isinteger(ifile,line,lp)
              if (.not.ok) call error('nciplot','bad fragment syntax',faterr)
              if (ifile < 1.or. ifile > nfiles) call error('nciplot','fragment: wrong file number',faterr)
              ok = isinteger(n1,line,lp)
              do while (ok)
                 if (n1<1.or.n1>m(ifile)%n) call error('nciplot','fragment: wrong atom number',faterr)
                 m(ifile)%ifrag(n1) = nfrag
                 ok = isinteger(n1,line,lp)
              end do
           else
              exit
           end if
        end do
     case ("INCREMENTS")
        read(line,*) xinc
        xinc = xinc / bohrtoa
     case ("CUTOFFS")
        read(line,*) rhocut, dimcut
     case ("CUTPLOT")
        read(line,*) rhoplot
     case ("ISORDG")
        read(line,*) isordg
     case ("RHOCUT2")
        read(line,*) rhoparam, rhoparam2
     case ("ELF")
        doelf = .true.
     case ("EXC")
        ! exchange
        read (line,*) wx, wc
        if (trim(wx) == 'S' .or. trim(wx) == 'SLATER' .or. &
            trim(wx) == 'LDA') then
           ixc(1) = 1
        else if (trim(wx) == 'PBE') then
           ixc(1) = 2
        else if (trim(wx) == 'B88') then
           ixc(1) = 3
        else if (trim(wx) == '-') then
           ixc(1) = -1
        elseif (trim(wx) == 'TEST') then
           ixc(1) = 99
        else
           call error('nciplot','Unknown exchange functional',faterr)
        end if
        ! correlation
        if (trim(wc) == 'PW' .or. trim(wc) == 'LDA' ) then
           ixc(2) = 1
        else if (trim(wc) == 'PBE') then
           ixc(2) = 2
        else if (trim(wc) == 'B88') then
           ixc(2) = 3
        else if (trim(wc) == '-') then
           ixc(2) = -1
        elseif (trim(wc) == 'TEST') then
           ixc(2) = 99
        else
           call error('nciplot','Unknown correlation functional',faterr)
        end if
     case ("DGRID")
        do i = 1, nfiles
           if (m(i)%ifile == ifile_xyz) m(i)%ifile = ifile_grd
        end do
     !case ("CHECKDEN")
     !   call checkden()
     case default
        call error('nciplot','Don''t know what to do with '//trim(word)//' keyword',faterr)
     end select
  enddo
11 continue

  ! set grid limits and npts
  if (autor) then
     if (ligand) then
        xinit = m(udat0)%x(:,1)
        xmax = m(udat0)%x(:,1)
        do j = 1, m(udat0)%n
           xinit = min(xinit,m(udat0)%x(:,j))
           xmax = max(xmax,m(udat0)%x(:,j))
        end do
     end if
     xinit = xinit - rthres
     xmax = xmax + rthres
  end if
  nstep = ceiling((xmax - xinit) / xinc)

  ! punch info
  write(uout,124)
  if (inter) write(uout,126) 
  if (ligand) write(uout,123) trim(m(udat0)%name)
  write(uout,120)
  write(uout,110) 'RHO  THRESHOLD   (au):', rhocut
  write(uout,110) 'RDG  THRESHOLD   (au):', dimcut
  if (inter) write(uout,110) 'DISCARDING RHO PARAM :',rhoparam
  if (ligand) write(uout,110) 'RADIAL THRESHOLD  (A):',rthres
  write (uout,*)
  write(uout,121) xinit, xmax, xinc, nstep

  ! open output files
  if (doelf) then
     luelf = 9
     open(luelf,file=trim(oname)//"-elf.cube")
  end if
  if (all(ixc /= 0)) then
     luxc = 8
     open(luxc,file=trim(oname)//"-xc.cube")
  end if
  if (noutput >= 2) then
     lugc = 11
     ludc = 12
     luvmd = 13
     open(lugc,file=trim(oname)//"-grad.cube")
     open(ludc,file=trim(oname)//"-dens.cube")
     open(luvmd,file=trim(oname)//".vmd")
  else
     lugc = -1
     ludc = -1
     luvmd = -1
  end if
  if (noutput == 1 .or. noutput == 3) then
     ludat = 14
     open(ludat,file=trim(oname)//".dat")
  else
     ludat = -1
  end if
  write(uout,122) trim(oname)//"-grad.cube",&
     trim(oname)//"-dens.cube",&
     trim(oname)//"-elf.cube",&
     trim(oname)//"-xc.cube",&
     trim(oname)//".vmd",&
     trim(oname)//".dat",&
     trim(oname)//".ncichk"
     
  ! write cube headers
  if (lugc > 0) call write_cube_header(lugc,'grad_cube','3d plot, reduced density gradient')
  if (ludc > 0) call write_cube_header(ludc,'dens_cube','3d plot, density')
  if (doelf) call write_cube_header(luelf,'elf_cube','3d plot, electron localisation function')
  if (all(ixc /= 0)) call write_cube_header(luxc,'xc_cube','3d plot, xc energy density')

  ! allocate memory for density and gradient
  allocate(crho(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1),stat=istat)
  if (istat /= 0) call error('nciplot','could not allocate memory for density cube',faterr)
  allocate(cgrad(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1),stat=istat)
  if (istat /= 0) call error('nciplot','could not allocate memory for grad',faterr)
  if (doelf) then
     allocate(celf(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1),stat=istat)
     if (istat /= 0) call error('nciplot','could not allocate memory for elf',faterr)
  end if
  if (all(ixc /= 0)) then
     allocate(cxc(0:nstep(1)-1,0:nstep(2)-1,0:nstep(3)-1),stat=istat)
     if (istat /= 0) call error('nciplot','could not allocate memory for xc',faterr)
  end if

  ! calculate density, rdg, elf, xc... read from chkpoint if available
  inquire(file=trim(oname)//".ncichk",exist=lchk)
10 continue
  if (lchk) then
     open(luchk,file=trim(oname)//".ncichk",form="unformatted")
     read (luchk) alcrho, alcgrad, alcelf, alcxc
     read (luchk) xinit0, xinc0, nstep0
     if ((alcrho.neqv.allocated(crho)) .or. (alcgrad.neqv.allocated(cgrad)) .or.&
         (alcelf.neqv.allocated(celf)) .or. (alcxc.neqv.allocated(cxc)) .or.&
         any(abs(xinit0 - xinit) > 1d-12) .or. any(abs(xinc0 - xinc) > 1d-12) .or.&
         any((nstep - nstep0) /= 0)) then
        lchk = .false.
        close(luchk)
        goto 10
     endif

     write(uout,'(" Reading the checkpoint file: ",A/)') trim(oname)//".ncichk"
     read (luchk) crho, cgrad
     if (allocated(celf)) read (luchk) celf
     if (allocated(cxc)) read (luchk) cxc
     close(luchk)
  else
     if (ispromol) then
        !$omp parallel do private (x,rho,grad,hess,heigs,hvecs,wk1,wk2,istat,grad2,&
        !$omp dimgrad,intra,rhom,elf,exc) schedule(dynamic)
        do k = 0, nstep(3)-1
           do j = 0, nstep(2)-1
              do i = 0, nstep(1)-1
                 x = xinit + (/i,j,k/) * xinc

                 ! calculate properties at x
                 call calcprops_pro(x,m,nfiles,rho,rhom(1:nfrag),nfrag,autofrag,&
                    grad,hess,doelf,elf,ixc,exc)
                 call rs(3,3,hess,heigs,0,hvecs,wk1,wk2,istat)
                 rho = max(rho,1d-30)
                 grad2 = dot_product(grad,grad)
                 dimgrad = sqrt(grad2) / (const*rho**(4.D0/3.D0))           

                 intra = inter .and. ((any(rhom(1:nfrag) >= sum(rhom(1:nfrag))*rhoparam)) .or. &
                    (sum(rhom(1:nfrag)) < rhoparam2 * rho))
                 if (intra) dimgrad = -dimgrad

                 !$omp critical (cubewrite)
                 crho(i,j,k) = sign(rho,heigs(2))*100.D0
                 cgrad(i,j,k) = dimgrad
                 if (doelf) celf(i,j,k) = elf
                 if (all(ixc /= 0)) cxc(i,j,k) = exc
                 !$omp end critical (cubewrite)
              end do
           end do
        end do
        !$omp end parallel do
     else
        call calcprops_wfn(xinit,xinc,nstep,m,nfiles,crho,cgrad,doelf,celf,ixc,cxc)
        if (inter) then
           !$omp parallel do private (x,rho,grad,hess,intra,rhom,elf,exc) schedule(dynamic)
           do k = 0, nstep(3)-1
              do j = 0, nstep(2)-1
                 do i = 0, nstep(1)-1
                    x = xinit + (/i,j,k/) * xinc
                    call calcprops_pro(x,m,nfiles,rho,rhom(1:nfrag),nfrag,autofrag,&
                       grad,hess,doelf,elf,ixc,exc)
                    intra = ((any(rhom(1:nfrag) >= sum(rhom(1:nfrag))*rhoparam)) .or. &
                       (sum(rhom(1:nfrag)) < rhoparam2 * rho))
                    !$omp critical (cubewrite)
                    if (intra) cgrad(i,j,k) = -abs(cgrad(i,j,k))
                    !$omp end critical (cubewrite)
                 enddo
              enddo
           enddo
           !$omp end parallel do
        endif
     endif

     ! save the ncichk file
     write(uout,'(" Writing the checkpoint file: ",A/)') trim(oname)//".ncichk"
     open(luchk,file=trim(oname)//".ncichk",form="unformatted")
     write (luchk) allocated(crho), allocated(cgrad), allocated(celf), allocated(cxc)
     write (luchk) xinit, xinc, nstep
     write(luchk) crho, cgrad
     if (allocated(celf)) write (luchk) celf
     if (allocated(cxc)) write (luchk) cxc
     close(luchk)
  endif

  ! apply cutoffs
  do k = 0, nstep(3)-1
     do j = 0, nstep(2)-1
        do i = 0, nstep(1)-1
           ! fragments for the wfn case
           intra = (cgrad(i,j,k) < 0)
           cgrad(i,j,k) = abs(cgrad(i,j,k))
           dimgrad = cgrad(i,j,k)
           rho = crho(i,j,k) / 100d0

           ! write the dat file
           if (ludat>0 .and. .not.intra .and. (abs(rho) < rhocut) .and. (dimgrad < dimcut) .and.&
              abs(rho)>1d-30) then
              write(ludat,'(1p,E18.10,E18.10)') rho, dimgrad
           end if ! rhocut/dimcut
           
           ! write the cube files
           if (all(ixc /= 0) .and. intra) cxc(i,j,k) = 100d0
           if (abs(rho) > rhoplot .or. intra) cgrad(i,j,k) = 100d0
        end do
     end do
  end do

  ! write cubes
  if (ludc > 0) call write_cube_body(ludc,nstep,crho)
  if (lugc > 0) call write_cube_body(lugc,nstep,cgrad)
  if (doelf) call write_cube_body(luelf,nstep,celf)
  if (all(ixc /= 0)) call write_cube_body(luxc,nstep,cxc)

  ! deallocate grids and close files
  if (allocated(crho)) deallocate(crho)
  if (allocated(cgrad)) deallocate(cgrad)
  if (allocated(celf)) deallocate(celf)
  if (allocated(cxc)) deallocate(cxc)
  if (ludat > 0) close(ludat)

  ! write vmd script
  if (ligand) then
     nn0 = sum(m(1:udat0-1)%n) + 1
     nnf = sum(m(1:udat0-1)%n) + m(udat0)%n
  else
     nn0 = 1
     nnf = ntotal
  end if
  if (luvmd > 0) then
     write (luvmd,114) trim(oname)//"-dens.cube"
     write (luvmd,115) trim(oname)//"-grad.cube"
     if (doelf) write (luvmd,115) trim(oname)//"-elf.cube"
     if (all(ixc /= 0)) write (luvmd,115) trim(oname)//"-xc.cube"
     write (luvmd,116) nn0-1,nnf-1,isordg,2,2,2,-rhoplot*100D0,rhoplot*100D0,2,2
     close(luvmd)
  end if

  ! end
  call tictac('End')

  ! close files
  if (uin /= stdin) close(uin)
  if (uout /= stdout) close(uout)

110 format (A,F5.2)
114 format ('#!/usr/local/bin/vmd',/,&
     '# VMD script written by save_state $Revision: 1.10 $',/,&
     '# VMD version: 1.8.6            ',/,&
     'set viewplist            ',/,&
     'set fixedlist            ',/,&
     '# Display settings            ',/,&
     'display projection   Orthographic            ',/,&
     'display nearclip set 0.000000            ',/,&
     '# load new molecule         ',/,&
     'mol new ',a,' type cube first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all')
115 format ('mol addfile ',a,' type cube first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all')
116 format ('#',/,&
     '# representation of the atoms',/,&
     'mol delrep 0 top',/,&
     'mol representation Lines 1.00000',/,&
     'mol color Name',/,&
     'mol selection {all}',/,&
     'mol material Opaque',/,&
     'mol addrep top',/,&
     'mol representation CPK 1.000000 0.300000 118.000000 131.000000',/,&
     'mol color Name',/,&
     'mol selection {index ', i5, ' to ', i5, ' }',/,&
     'mol material Opaque',/,&
     'mol addrep top',/,&
     '#',/,&
     '# add representation of the surface',/,&
     'mol representation Isosurface ',f7.5,' 1 0 0 1 1',/,&
     'mol color Volume 0',/,&
     'mol selection {all}',/,&
     'mol material Opaque',/,&
     'mol addrep top',/,&
     'mol selupdate ',i1,' top 0',/,&
     'mol colupdate ',i1,' top 0',/,&
     'mol scaleminmax top ',i1,' ', f7.4, f7.4,/,&
     'mol smoothrep top ', i1,' 0',/,&
     'mol drawframes top ',i1,' {now}',/,&
     'color scale method BGR',/,&
     'set colorcmds {{color Name {C} gray}}',/,&
     '#some more',/)
120 format(/'-----------------------------------------------------'/&
            '      Calculation details:'/&
            '-----------------------------------------------------')
121 format(/,'-----------------------------------------------------'/&
             '      Operating grid and increments:'/&
             '-----------------------------------------------------'/&
             ' x0,y0,z0  = ',f8.4,' ',f8.4,' ',f8.4/&
             ' x1,y1,z1  = ',f8.4,' ',f8.4,' ',f8.4/&
             ' ix,iy,iz  = ',f5.2,'   ',f5.2,'   ',f5.2/&
             ' nx,ny,nz  = ',i4,'    ', i4,'    ', i4/)
122 format('-----------------------------------------------------'/&
           '      Writing output in the following units:'/&
           '-----------------------------------------------------'/&
           ' Reduced Density Gradient,RDG  = ',a,/&
           ' Sign(lambda2)xDensity,LS      = ',a,/&
           ' ELF cube file                 = ',a,/&
           ' XC energy density cube file   = ',a,/&
           ' LS x RDG                      = ',a,/&
           ' VMD script                    = ',a,/&
           ' NCI checkpoint                = ',a,/&
           '-----------------------------------------------------',/)
123 format('      Using ',a40,' as LIGAND')
124 format('-----------------------------------------------------'/&
           '      INPUT INFORMATION:'/&
           '-----------------------------------------------------')
126 format(/'      MIND YOU'/&
            '      ONLY ANALYZING INTERMOLECULAR INTERACTIONS     '/)

contains
  subroutine write_cube_header(lu,l1,l2)

    integer, intent(in) :: lu
    character*(*), intent(in) :: l1, l2

    integer :: i, j

    write(lu,*) trim(l1)
    write(lu,*) trim(l2)
    write(lu,'(I5,3(F12.6))') ntotal, xinit
    write(lu,'(I5,3(F12.6))') nstep(1), xinc(1), 0d0, 0d0
    write(lu,'(I5,3(F12.6))') nstep(2), 0d0, xinc(2), 0d0
    write(lu,'(I5,3(F12.6))') nstep(3), 0d0, 0d0, xinc(3)
    do i = 1, nfiles
       do j = 1, m(i)%n
          write(lu,'(I4,F5.1,F11.6,F11.6,F11.6)') m(i)%z(j), 0d0, m(i)%x(:,j)
       end do
    end do

  end subroutine write_cube_header

  subroutine write_cube_body(lu,n,c)
    
    integer, intent(in) :: lu
    integer, intent(in) :: n(3)
    real*8, intent(in) :: c(0:n(1)-1,0:n(2)-1,0:n(3)-1)

    integer :: i, j

    do i = 0, n(1)-1
       do j = 0, n(2)-1
          write (lu,'(6(1x,e12.5))') (c(i,j,k),k=0,n(3)-1)
       enddo
    enddo
    close(lu)

  end subroutine write_cube_body

end program
