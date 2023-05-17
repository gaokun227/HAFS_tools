!========================================================================================
  subroutine hafsvi_preproc_ic(in_dir, in_date, nestdoms, radius, res, out_file)

!-----------------------------------------------------------------------------
! HAFS DA tool - hafsvi_preproc
! Yonghui Weng, 20211210
!
! This subroutine read hafs restart files and output hafsvi needed input.
! Variables needed:
!      WRITE(IUNIT) NX,NY,NZ
!      WRITE(IUNIT) lon1,lat1,lon2,lat2,cen_lon,cen_lat
!      WRITE(IUNIT) (((pf1(i,j,k),i=1,nx),j=1,ny),k=nz,1,-1)                   ! 3D, NZ
!      WRITE(IUNIT) (((tmp(i,j,k),i=1,nx),j=1,ny),k=nz,1,-1)
!      WRITE(IUNIT) (((spfh(i,j,k),i=1,nx),j=1,ny),k=nz,1,-1)
!      WRITE(IUNIT) (((ugrd(i,j,k),i=1,nx),j=1,ny),k=nz,1,-1)
!      WRITE(IUNIT) (((vgrd(i,j,k),i=1,nx),j=1,ny),k=nz,1,-1)
!      WRITE(IUNIT) (((dzdt(i,j,k),i=1,nx),j=1,ny),k=nz,1,-1)
!!      WRITE(IUNIT) hgtsfc                ! 2D
!      WRITE(IUNIT) (((z1(i,j,k),i=1,nx),j=1,ny),k=nz1,1,-1)
!      WRITE(IUNIT) glon,glat,glon,glat   ! 2D
!      WRITE(IUNIT) (((ph1(i,j,k),i=1,nx),j=1,ny),k=nz1,1,-1)                ! 3D, NZ+1
!      WRITE(IUNIT) pressfc1              ! 2D
!      WRITE(IUNIT) ak
!      WRITE(IUNIT) bk
!      WRITE(IUNIT) land                  ! =A101 = land sea mask, B101 = ZNT
!      WRITE(IUNIT) sfcr                  ! =B101 = Z0
!      WRITE(IUNIT) C101                  ! =C101 = (10m wind speed)/(level 1 wind speed)
!

!-----------------------------------------------------------------------------

  use constants
  use netcdf
  use module_mpi
  use var_type

  implicit none

  character (len=*), intent(in) :: in_dir, in_date, radius, res, out_file
  integer, intent(in)           :: nestdoms
!--- in_dir,  HAFS_restart_folder, which holds grid_spec.nc, fv_core.res.tile1.nc,
!             fv_srf_wnd.res.tile1.nc, fv_tracer.res.tile1.nc, phy_data.nc, sfc_data.nc
!--- in_date, HAFS_restart file date, like 20200825.120000
!--- radius,  to cut a square, default value is 40, which means a 40deg x 40deg square.
!--- out_file: output file, default is bin format, if the file name is *.nc, then output nc format.
!--- nestdoms: total nest domain number: 0-no nesting
!---                                     1-nest02.tile2 + 0
!---                                     2-nest03.tile3 + 1

  character (len=2500)   :: indir, infile
  character (len=2500)   :: infile_fvcore, infile_core, infile_uv, infile_tracer, infile_phy, &
                            infile_sfc, infile_grid, infile_grid2, infile_vertial, infile_oro
  type(grid2d_info)      :: dstgrid  ! rot-ll grid for output
  type(grid2d_info)      :: ingrid   ! hafs restart grid
  real     :: radiusf
  logical  :: file_exist

!----for hafs restart
  integer  :: ix, iy, iz, kz, ndom, nd
  character (len=50) :: nestfl, tilefl, tempfl
                        ! grid_spec.nc : grid_spec.nest02.tile2.nc
                        ! fv_core.res.tile1.nc : fv_core.res.nest02.tile2.nc
                        ! phy_data.nc  : phy_data.nest02.tile2.nc

!----for hafsvi
  integer  :: nx, ny, nz, filetype  ! filetype: 1=bin, 2=nc
  real     :: lon1,lat1,lon2,lat2,cen_lat,cen_lon,dlat,dlon
  real, allocatable, dimension(:,:) :: glon,glat

  integer  :: i, j, k, flid_in, flid_out, ncid, ndims, nrecord
  real     :: rot_lon, rot_lat, ptop
  integer, dimension(nf90_max_var_dims) :: dims
  real, allocatable, dimension(:,:,:,:) :: dat4, dat41, dat42, dat43, u, v
  real, allocatable, dimension(:,:,:)   :: dat3, dat31
  real, allocatable, dimension(:,:)     :: dat2, dat21, sfcp
  real, allocatable, dimension(:)       :: dat1

  !real, allocatable, dimension(:)       :: pfull, phalf
  real, allocatable, dimension(:,:)     :: cangu, sangu, cangv, sangv
  real    :: cputime1, cputime2, cputime3
  integer :: io_proc, nm, ks, ke, nv


!------------------------------------------------------------------------------
! 1 --- arg process
  io_proc=nprocs-1
  !io_proc=0
!
! 1.1 --- ndom
  ndom=nestdoms+1

! 1.2 --- input_dir
  if (len_trim(in_dir) < 2 .or. trim(in_dir) == 'w' .or. trim(in_dir) == 'null') then
     indir='.'
  else
     indir=trim(in_dir)
  endif

  if (trim(radius) == 'w' .or. trim(radius) == 'null') then
     radiusf = 40.  !deg
  else
     read(radius,*)i
     radiusf = real(i)
     if ( radiusf < 3. .or. radiusf > 70. ) then
        if ( my_proc_id == 0 ) write(*,'(a)')'!!! hafsvi cut radius number wrong: '//trim(radius)
        if ( my_proc_id == 0 ) write(*,'(a)')'!!! please call with --vortexradius=40 (75< 3)'
        stop 'hafsvi_preproc'
     endif
  endif

  if (trim(res) == 'w' .or. trim(res) == 'null') then
     dlat=0.02
  else
     read(res,*)dlat
  endif
  dlon=dlat

!------------------------------------------------------------------------------
! 2 --- set dstgrid: rot-ll grid
! 2.1 --- define rot-ll grid
  cen_lat = tc%lat
  cen_lon = tc%lon
  nx = int(radiusf/2.0/dlon+0.5)*2+1
  ny = int(radiusf/2.0/dlat+0.5)*2+1
  lon1 = - radiusf/2.0
  lat1 = - radiusf/2.0
  lon2 = radiusf/2.0
  lat2 = radiusf/2.0
  !!--- get rot-ll grid
  allocate(glon(nx,ny), glat(nx,ny))
  !$omp parallel do &
  !$omp& private(i,j,rot_lon,rot_lat)
  do j = 1, ny; do i = 1, nx
     rot_lon = lon1 + dlon*(i-1)
     rot_lat = lat1 + dlat*(j-1)
     call rtll(rot_lon, rot_lat, glon(i,j), glat(i,j), cen_lon, cen_lat)
  enddo; enddo
  if ( my_proc_id == 0 ) write(*,'(a)')'---rot-ll grid: nx, ny, cen_lon, cen_lat, dlon, dlat, lon1, lon2, lat1, lat2'
  if ( my_proc_id == 0 ) write(*,'(15x,2i5,8f10.5)')    nx, ny, cen_lon, cen_lat, dlon, dlat, lon1, lon2, lat1, lat2
  !write(*,'(a,4f10.5)')'---rot-ll grid rot_lon:', glon(1,1), glon(1,ny), glon(nx,ny), glon(nx,1)
  !write(*,'(a,4f10.5)')'---rot-ll grid rot_lat:', glat(1,1), glat(1,ny), glat(nx,ny), glat(nx,1)

! 2.2 --- set dstgrid
  dstgrid%grid_x = nx
  dstgrid%grid_y = ny
  dstgrid%ntime  = 1
  dstgrid%grid_xt = nx
  dstgrid%grid_yt = ny
  allocate(dstgrid%grid_lon (dstgrid%grid_x,dstgrid%grid_y))
  allocate(dstgrid%grid_lont(dstgrid%grid_x,dstgrid%grid_y))
  dstgrid%grid_lon  = glon
  dstgrid%grid_lont = glon
  allocate(dstgrid%grid_lat (dstgrid%grid_x,dstgrid%grid_y))
  allocate(dstgrid%grid_latt(dstgrid%grid_x,dstgrid%grid_y))
  dstgrid%grid_lat  = glat
  dstgrid%grid_latt = glat

!------------------------------------------------------------------------------
! 3 --- process output file type: now is only for bin
!  i=len_trim(out_file)
!  if ( out_file(i-2:i) == '.nc' ) then
!     write(*,'(a)')' --- output to '//trim(out_file)
!     filetype=2
!     call nccheck(nf90_open(trim(out_file), nf90_write, flid), 'wrong in open '//trim(out_file), .true.)
!  else
!     filetype=1
!     flid=71
!     open(unit=flid,file=trim(out_file),form='unformatted',status='unknown')
!  endif

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! --- domain loop: from inner domain to outer domain, so the number is from max to 1
  do_nestdom_loop: do nd = 1, ndom

     !-------------------------------------------------------------------------
     ! 3 --- initialization: clean ingrid, weight
     ! ingrid%grid_x=-99; ingrid%grid_y=-99; ingrid%grid_xt=-99; ingrid%grid_yt=-99

     !-------------------------------------------------------------------------
     ! 4 --- input grid info
     !       read from grid file grid_spec.nc:
     !       nestfl, tilefl: infile_core, infile_tracer, infile_grid, infile_vertial, infile_oro
     write(nestfl,'(a4,i2.2)')'nest',nd
     write(tilefl,'(a4,i0)')'tile',nd

     ! KGao - name changes
     !if ( nd == 0 ) then
     infile_grid=trim(indir)//'/grid_spec.nc'
     !infile_oro =trim(indir)//'/oro_data.nc'
     infile_vertical=trim(indir)//'/gfs_vertical.nc' ! ak, bk, pfull, phalf
     infile_core=trim(indir)//'/gfs_data.nc'
     infile_tracer=trim(indir)//'/gfs_data.nc'
     infile_sfc =trim(indir)//'/sfc_data.nc'
     infile_uv=trim(indir)//'/gfs_uv_agrid.nc'
     !endif

     ! KGao
     !inquire(file=infile_grid2, exist=file_exist)
     !if ( file_exist ) infile_grid = infile_grid2

     if ( debug_level > 10 .and. my_proc_id == 0 ) write(*,'(a)')' --- read grid info from '//trim(infile_grid)
     call rd_grid_spec_data(trim(infile_grid), ingrid)
     ix=ingrid%grid_xt
     iy=ingrid%grid_yt
     if ( debug_level > 10 .and. my_proc_id == 0 ) then
        write(*,'(a,i,1x,i,1x,f,1x,f,1x,f,1x,f)')' --- ingrid info: ', ix, iy, &
              ingrid%grid_lon(int(ix/2), int(iy/2)), ingrid%grid_lat(int(ix/2), int(iy/2)), &
              ingrid%grid_lont(int(ix/2), int(iy/2)), ingrid%grid_latt(int(ix/2), int(iy/2))
     endif

     !---to add the test if the tc was inside of the domain

     ! KGao
     ! call FV3-grid cos and sin
     !allocate( cangu(ix,iy+1),sangu(ix,iy+1),cangv(ix+1,iy),sangv(ix+1,iy) )
     !call cal_uv_coeff_fv3(ix, iy, ingrid%grid_lat, ingrid%grid_lon, cangu, sangu, cangv, sangv)

     !-------------------------------------------------------------------------
     ! 5 --- calculate output-grid in input-grid's positions (xin, yin), and each grid's weight to dst
     if ( debug_level > 10 ) then
        if ( my_proc_id == 0 ) write(*,'(a)')' --- call cal_src_dst_grid_weight'
        write(*,'(i,a,2(i,1x),2(f,1x))')my_proc_id,' --- dstgrid: ', nx, ny, &
             dstgrid%grid_lont(int(nx/2),int(ny/2)), dstgrid%grid_latt(int(nx/2),int(ny/2))
     endif
     call cal_src_dst_grid_weight(ingrid, dstgrid)

     !-------------------------------------------------------------------------
     ! 6 --- dst files
     if ( my_proc_id == io_proc ) then
        flid_in=71   !inner domain rot-ll file
        flid_out=72  !current domain rot-ll file
        if ( nd == 0 ) then
           open(unit=flid_out,file=trim(out_file),form='unformatted',status='unknown')
        endif
     endif

     !-------------------------------------------------------------------------
     ! 7 --- output
     do_out_var_loop: do nrecord = 1, 17
        !write(*,*)my_proc_id, '=== nrecord =',nrecord
        !-----------------------------
        !---7.1 record 1: nx, ny, nz
        !---nx, ny, nz, & lon1,lat1,lon2,lat2,cen_lon,cen_lat
        call cpu_time(cputime1)
        if ( my_proc_id == io_proc ) write(*,'(a,i3,f)')' --- record start cputime: ', nrecord, cputime1
        if ( nrecord == 1 ) then

           ! KGao: get nz dim
           !call get_var_dim(trim(infile_vertial), 'pfull', ndims, dims)
           !nz=dims(1)
           nz=128

           if ( my_proc_id == io_proc ) write(*,'(a,3i6)')'=== record1: ',nx, ny, nz
           if ( my_proc_id == io_proc ) write(flid_out) nx, ny, nz
           if ( nd > 1 .and. my_proc_id == io_proc ) read(flid_in)
        endif

        iz=nz   !same vertical levels
        if ( nrecord == 12 .or. nrecord == 15 .or. nrecord == 16 .or. nrecord ==17 ) iz=1
        if ( nrecord ==  9 .or. nrecord == 11 ) iz=nz+1

        !-----------------------------
        !---7.2 record 2: lon1,lat1,lon2,lat2,cen_lon,cen_lat
        if ( nrecord == 2 ) then
           write(*,'(i,a,6f8.3)')my_proc_id, '=== record2: ',lon1,lat1,lon2,lat2,cen_lon,cen_lat
           if ( my_proc_id == io_proc ) write(*,'(a,6f8.3)')'=== record2: ',lon1,lat1,lon2,lat2,cen_lon,cen_lat
           if ( my_proc_id == io_proc ) write(flid_out) lon1,lat1,lon2,lat2,cen_lon,cen_lat
           if ( nd > 1 .and. my_proc_id == io_proc ) read(flid_in)
           write(*,*)'==== finished record 2 at ', my_proc_id
        endif

        !-----------------------------
        !---7.3 record 3: (((pf1(i,j,k),i=1,nx),j=1,ny),k=nz,1,-1)
        !---     hafs-VI/read_hafs_out.f90 pf:
        !---          ph(k) = ak(k) + bk(k)*p_s --> half level pressure
        !---          pf(k) = (ph(k+1) - ph(k)) / log(ph(k+1)/ph(k)) --> full level pressure
        !---
        !---seem pf1 is pressure on full level, use
        !---    pf1(k) = phalf(1) + sum(delp(1:k))
        if ( nrecord == 3 .or. nrecord == 9 .or. nrecord == 11 ) then
           if ( nrecord == 3 ) then
              !allocate(dat4(iz+1,1,1,1))
              !call get_var_data(trim(infile_vertial), 'phalf', iz+1, 1, 1, 1, dat4)
              !ptop=dat4(1,1,1,1)*100.  !phalf:units = "mb" ;
              !deallocate(dat4)

              ! KGao - get model top
              ptop = 0*100 !  
           endif

           if ( my_proc_id == io_proc ) then
              if ( nrecord == 3 ) then
                 allocate(dat4(ix, iy, iz,1))
                 allocate(dat41(ix, iy, iz,1))
                 allocate(dat2(ix, iy))
                 !write(*,'(a,3i5)')'delp: ',ix, iy, iz
                 call get_var_data(trim(infile_core), 'delp', ix, iy, iz,1, dat4)
                 dat2(:,:)=ptop
                 do k = 1, iz
                    dat41(:,:,k,1)=dat2(:,:)+dat4(:,:,k,1)/2.0
                    dat2(:,:)=dat2(:,:)+dat4(:,:,k,1)
                 enddo
                 allocate(sfcp(ix, iy))
                 sfcp=dat41(:,:,iz,1)
                 deallocate(dat2, dat4)
              else if ( nrecord == 9 ) then
                 !allocate(dat4(ix, iy, 1,1))
                 !call get_var_data(trim(infile_core), 'phis', ix, iy, 1, 1, dat4)
                 !allocate(dat41(ix, iy, iz, 1))
                 !dat41(:,:,iz,1)=dat4(:,:,1,1)/g ! KGao note: surface height
                 !deallocate(dat4)
                 !allocate(dat4(ix, iy, iz-1, 1))
                 !call get_var_data(trim(infile_core), 'DZ', ix, iy, iz-1, 1, dat4)
                 !do k = iz-1, 1, -1
                 !   dat41(:,:,k,1)=dat41(:,:,k+1,1)-dat4(:,:,k,1) ! KGao note: layer height above surface
                 !enddo
                 !deallocate(dat4)

                 allocate(dat41(ix, iy, iz, 1))
                 call get_var_data(trim(infile_core), 'zh', ix, iy, iz, 1, dat41)


              else if ( nrecord == 11 ) then
                 allocate(dat4(ix, iy, iz-1, 1), dat41(ix, iy, iz, 1))
                 call get_var_data(trim(infile_core), 'delp', ix, iy, iz-1, 1, dat4)
                 dat41(:,:,1,1)=ptop
                 do k = 2, iz
                    dat41(:,:,k,1)=dat41(:,:,k-1,1)+dat4(:,:,k-1,1)
                 enddo
                 deallocate(dat4)
              endif  !if ( nrecord == 3 ) then

              !---broadcast dat41 to each computing-core
              if ( nprocs > 1 ) then
                 nm=max(1,int((iz+nprocs-1)/nprocs))  !devide iz to each processor
                 do k = 0, nprocs-1
                    ks=k*nm+1              !k-start
                    ke=k*nm+nm            !k-end
                    if ( ke > iz ) ke=iz
                    if ( ks >= 1 .and. ks <= iz .and. ke >= 1 .and. ke <= iz ) then
                       allocate(dat42(ix, iy, ke-ks+1,1))
                       dat42(:,:,:,1)=dat41(:,:,ks:ke,1)
                       if ( k /= io_proc ) then
                          call mpi_send(dat42(1,1,1,1), size(dat42), mpi_real, k, 3000+ks, comm, ierr)
                       else
                          allocate(dat43(ix, iy, ke-ks+1,1))
                          dat43=dat42
                       endif
                       deallocate(dat42)
                    endif
                 enddo
              else  !if ( nprocs > 1 ) then
                 allocate(dat43(ix, iy, iz,1))
                 dat43=dat41
              endif
              deallocate(dat41)
           else ! if ( my_proc_id == io_proc ) then
              !---receive dat43
              nm=max(1,int((iz+nprocs-1)/nprocs))
              ks=min(iz,my_proc_id*nm+1)
              ke=min(iz,my_proc_id*nm+nm)
              if ( ks >= 1 .and. ks <= iz .and. ke >= 1 .and. ke <= iz ) then
                 allocate(dat43(ix, iy, ke-ks+1,1))
                 !call mpi_recv(dat43(1,1,1,1), size(dat43), mpi_real, io_proc, mpi_any_tag, comm, status, ierr)
                 call mpi_recv(dat43(1,1,1,1), size(dat43), mpi_real, io_proc, 3000+ks, comm, status, ierr)
              endif
           endif
        endif

        !-----------------------------
        !---7.4 record 4: (((tmp(i,j,k),i=1,nx),j=1,ny),k=nz,1,-1)
        !---7.5 record 5: (((spfh(i,j,k),i=1,nx),j=1,ny),k=nz,1,-1)
        !---7.8 record 8: (((dzdt(i,j,k),i=1,nx),j=1,ny),k=nz,1,-1)
        if ( nrecord == 4 .or. nrecord == 5 .or. nrecord == 8 ) then
           if ( my_proc_id == io_proc ) then
              !---read in
              allocate(dat4(ix, iy, iz,1))
              ! KGao: 'T' and 'W' to 't' and 'w'
              if ( nrecord == 4 ) call get_var_data(trim(infile_core), 't', ix, iy, iz,1, dat4)
              if ( nrecord == 5 ) call get_var_data(trim(infile_tracer), 'sphum', ix, iy, iz,1, dat4)
              if ( nrecord == 8 ) call get_var_data(trim(infile_core), 'w', ix, iy, iz,1, dat4)
              !---send to other core
              if ( nprocs > 1 ) then
                 nm=max(1,int((iz+nprocs-1)/nprocs))
                 do k = 0, nprocs-1
                    ks=k*nm+1
                    ke=k*nm+nm            !k-end
                    if ( ke > iz ) ke=iz
                    if ( ks >= 1 .and. ks <= iz .and. ke >= 1 .and. ke <= iz ) then
                       allocate(dat42(ix, iy, ke-ks+1,1))
                       dat42(:,:,:,1)=dat4(:,:,ks:ke,1)
                       if ( k /= io_proc ) then
                          call mpi_send(dat42(1,1,1,1), size(dat42), mpi_real, k, 4000+ks, comm, ierr)
                       else
                          allocate(dat43(ix, iy, ke-ks+1,1))
                          dat43=dat42
                       endif
                       deallocate(dat42)
                    endif
                 enddo
              else
                 allocate(dat43(ix, iy, iz,1))
                 dat43=dat4
              endif
              deallocate(dat4)
           else  !if ( my_proc_id == io_proc ) then
              !---receive dat43
              nm=max(1,int((iz+nprocs-1)/nprocs))
              ks=min(iz,my_proc_id*nm+1)
              ke=min(iz,my_proc_id*nm+nm)
              if ( ks >= 1 .and. ks <= iz .and. ke >= 1 .and. ke <= iz ) then
                 allocate(dat43(ix, iy, ke-ks+1,1))
                 call mpi_recv(dat43(1,1,1,1), size(dat43), mpi_real, io_proc, 4000+ks, comm, status, ierr)
              endif
           endif  !if ( my_proc_id == io_proc ) then
        endif

        !-----------------------------
        !---7.6 record 6: (((ugrd(i,j,k),i=1,nx),j=1,ny),k=nz,1,-1)
        !---7.7 record 7: (((vgrd(i,j,k),i=1,nx),j=1,ny),k=nz,1,-1)
        if ( nrecord == 6 ) then
           !---get u,v from restart
           if ( my_proc_id == io_proc ) then
              do nv = 1, 2
                 if (nv==1) then
                    ! KGao
                    !allocate(dat4(ix, iy+1, iz,1))
                    allocate(dat4(ix, iy, iz,1))
                    call get_var_data(trim(infile_uv), 'u', ix, iy, iz, 1, dat4)
                 else if (nv==2) then
                    ! KGao
                    !allocate(dat4(ix+1, iy, iz,1))
                    allocate(dat4(ix, iy, iz,1))
                    call get_var_data(trim(infile_uv), 'v', ix, iy, iz, 1, dat4)
                 endif
                 !---send to other core
                 if ( nprocs > 1 ) then
                    nm=max(1,int((iz+nprocs-1)/nprocs))
                    do k = 0, nprocs-1
                       ks=k*nm+1
                       ke=k*nm+nm            !k-end
                       if ( ke > iz ) ke=iz
                       if ( ks >= 1 .and. ks <= iz .and. ke >= 1 .and. ke <= iz ) then
                          if (nv==1) then
                             ! KGao
                             !allocate(dat41(ix, iy+1, ke-ks+1,1))
                             allocate(dat41(ix, iy, ke-ks+1,1))
                          else if (nv==2) then
                             ! KGao
                             !allocate(dat41(ix+1, iy, ke-ks+1,1))
                             allocate(dat41(ix, iy, ke-ks+1,1))
                          endif
                          dat41(:,:,:,1)=dat4(:,:,ks:ke,1)
                          if ( k /= io_proc ) then
                             call mpi_send(dat41(1,1,1,1), size(dat41), mpi_real, k, 200*nv+ks, comm, ierr)
                          else
                             if (nv==1) then
                                ! KGao
                                !allocate(dat42(ix, iy+1, ke-ks+1,1))
                                allocate(dat42(ix, iy, ke-ks+1,1))
                                dat42=dat41
                             else if (nv==2) then
                                ! KGao
                                !allocate(dat43(ix+1, iy, ke-ks+1,1))
                                allocate(dat43(ix, iy, ke-ks+1,1))
                                dat43=dat41
                             endif
                          endif
                          deallocate(dat41)
                       endif
                    enddo
                 else  !if ( nprocs > 1 ) then
                    if (nv==1) then
                       ! KGao
                       !allocate(dat42(ix, iy+1, iz,1))
                       allocate(dat42(ix, iy, iz,1))
                       dat42=dat4
                    else if (nv==2) then
                       ! KGao
                       !allocate(dat43(ix+1, iy, iz,1))
                       allocate(dat43(ix, iy, iz,1))
                       dat43=dat4
                    endif
                 endif
                 deallocate(dat4)
              enddo  !do nv = 1, 2
           else  !if ( my_proc_id == io_proc ) then
              !---receive dat42 dat43
              nm=max(1,int((iz+nprocs-1)/nprocs))
              ks=min(iz,my_proc_id*nm+1)
              ke=min(iz,my_proc_id*nm+nm)
              if ( ks >= 1 .and. ks <= iz .and. ke >= 1 .and. ke <= iz ) then

                 ! KGao
                 !allocate(dat42(ix, iy+1, ke-ks+1,1), dat43(ix+1, iy, ke-ks+1,1))
                 allocate(dat42(ix, iy, ke-ks+1,1), dat43(ix, iy, ke-ks+1,1))

                 call mpi_recv(dat42(1,1,1,1), size(dat42), mpi_real, io_proc, 200*1+ks, comm, status, ierr)
                 call mpi_recv(dat43(1,1,1,1), size(dat43), mpi_real, io_proc, 200*2+ks, comm, status, ierr)
              endif
           endif  !if ( my_proc_id == io_proc ) then
           !write(*,*)'===w11 distributed u,v @ dat42,dat43'
        endif

        !!-----------------------------
        !!---7.9 record 9: (((z1(i,j,k),i=1,nx),j=1,ny),k=nz1,1,-1)
        !!---     hafs-VI/read_hafs_out.f90 z1:
        !!---          z1(I,J,K)=z1(I,J,K+1)+rdgas1*tmp(i,j,k)*(1.+0.608*spfh(i,j,k))*ALOG(ph1(i,j,k+1)/ph1(i,j,k))
        !!--- hgt: phis/g-sum(DZ)
        !if ( nrecord == 9 ) then
        !   allocate(dat4(ix, iy, 1,1))
        !   call get_var_data(trim(infile_core), 'phis', ix, iy, 1, 1, dat4)
        !   allocate(dat41(ix, iy, iz+1, 1))
        !   dat41(:,:,iz+1,1)=dat4(:,:,1,1)/g
        !   deallocate(dat4)

        !   allocate(dat4(ix, iy, iz, 1))
        !   call get_var_data(trim(infile_core), 'DZ', ix, iy, iz, 1, dat4)
        !   do k = iz, 1, -1
        !      dat41(:,:,k,1)=dat41(:,:,k+1,1)-dat4(:,:,k,1)
        !   enddo
        !   !write(*,'(a,200f)')'z1: ',dat41(int(ix/2),int(iy/2),:,1)
        !   deallocate(dat4)
        !endif

        !-----------------------------
        !---7.10 record 10: glon,glat,glon,glat   ! 2D
        !--- glat=grid_yt*180./pi, grid_yt=1:2160, what is this?
        if ( nrecord == 10 ) then
           if ( my_proc_id == io_proc ) then
              write(*,'(a,4f8.3)')'=== record10: ',glon(1,1), glat(1,1), glon(nx,ny), glat(nx,ny)

! New change: Jul 2022 JH Shin-------------------------------------------------------------

          !For WPAC storms located near the west of international date line, add 360 to
          !longitude value in the eastern side of IDL (western hemisphere) so that all
          !longitude values have positive values in the VI domain.
          !if ( cen_lon > 0. ) then
          ! do j = 1, ny; do i = 1, nx
          !  if(glon(i,j).lt.0.0) glon(i,j)=glon(i,j)+360.0
          ! enddo; enddo
          !endif
          if ( cen_lon > 0. ) where ( glon < 0.) glon=glon+360.

          !For CPAC storms located near the east of international date line, subrtact 360
          !from longitude value in the western side of IDL (eastern hemisphere) so
          !that all longitude values have negative values in the VI domain.
          !if ( cen_lon < -140. )then
          ! do j = 1, ny; do i = 1, nx
          !  if(glon(i,j).gt.0.0) glon(i,j)=glon(i,j)-360.0
          ! enddo; enddo
          !endif
          if ( cen_lon < -140. ) where ( glon > 0. ) glon=glon-360.

! New change: Jul 2022 JH Shin-------------------------------------------------------------

              write(flid_out) glon,glat,glon,glat
              if ( nd > 1 ) read(flid_in)
           endif
        endif

        !!-----------------------------
        !!---7.11 record 11: (((ph1(i,j,k),i=1,nx),j=1,ny),k=nz1,1,-1)
        !!---     hafs-VI/read_hafs_out.f90 ph:
        !!---       ph(k) = ak(k) + bk(k)*p_s --> pressure in pa
        !!---       64.270-->100570
        !!---seem ph1 is pressure on half level, use
        !!---    pf1(k) = phalf(1) + sum(delp(1:k))
        !if ( nrecord == 11 ) then
        !   allocate(dat4(iz+1,1,1,1))
        !   call get_var_data(trim(infile_vertial), 'phalf', iz+1, 1, 1, 1, dat4)
        !   ptop=dat4(1,1,1,1)*100.  !phalf:units = "mb" ;
        !   deallocate(dat4)

        !   allocate(dat4(ix, iy, iz,1))
        !   allocate(dat41(ix, iy, iz+1,1))
        !   call get_var_data(trim(infile_core), 'delp', ix, iy, iz,1, dat4)
        !   dat41(:,:,1,1)=ptop
        !   do k = 1, iz
        !      dat41(:,:,k+1,1)=dat41(:,:,k,1)+dat4(:,:,k,1)
        !   enddo
        !   deallocate(dat4)
        !endif

        !-----------------------------
        !---7.12 record 12: pressfc1              ! 2D
        !--- use lowest-level pressure?
        if ( nrecord == 12 ) then
           if ( my_proc_id == io_proc ) then
              allocate(dat43(ix, iy, 1, 1))
              dat43(:,:,1,1)=sfcp(:,:)
              deallocate(sfcp)
           endif
        endif

        !-----------------------------
        !---7.13 record 13: ak
        if ( nrecord == 13 ) then
           if ( my_proc_id == io_proc ) then
              allocate(dat4(iz+1,2,1,1))
              ! KGao: ak->vcoord
              call get_var_data(trim(infile_vertical), 'vcoord', iz+1, 2, 1, 1, dat4)
              write(*,'(a,200f12.1)')'=== record13: ', (dat4(k,1,1,1),k=1,iz+1)
              write(flid_out) (dat4(k,1,1,1),k=1,iz+1)
              if ( nd > 1 ) read(flid_in)
              deallocate(dat4)
           endif
        endif

        !-----------------------------
        !---7.14 record 14: bk
        if ( nrecord == 14 ) then
           if ( my_proc_id == io_proc ) then
              ! KGao: bk-> vcoord
              allocate(dat4(iz+1,2,1,1))
              call get_var_data(trim(infile_vertical), 'vcoord', iz+1, 2, 1, 1, dat4)
              write(*,'(a,200f10.3)')'=== record14: ', (dat4(k,2,1,1),k=1,iz+1)
              write(flid_out) (dat4(k,2,1,1),k=1,iz+1)
              if ( nd > 1 ) read(flid_in)
              deallocate(dat4)
           endif
        endif

        !-----------------------------
        !---7.15 record 15: land                  ! =A101 = land sea mask, B101 = ZNT
        !---     hafs-VI/read_hafs_out.f90 land:long_name = "sea-land-ice mask (0-sea, 1-land, 2-ice)" ;
        !--- sfc_data.nc: slmsk
        if ( nrecord == 15 ) then
           if ( my_proc_id == io_proc ) then
              allocate(dat43(ix, iy, 1,1))
              call get_var_data(trim(infile_sfc), 'slmsk', ix, iy, 1, 1, dat43)
           endif
        endif

        !-----------------------------
        !---7.16 record 16: sfcr                  ! =B101 = Z0
        !---surface roughness
        if ( nrecord == 16 ) then
           if ( my_proc_id == io_proc ) then
              allocate(dat43(ix, iy, 1,1))
              call get_var_data(trim(infile_sfc), 'zorl', ix, iy, 1, 1, dat43)
              ! convert from cm to m
              dat43=dat43/100.
           endif
        endif

        !-----------------------------
        !---7.17 record 17: C101                  ! =C101 = (10m wind speed)/(level 1 wind speed)
        !---                                      ! =C101 = f10m (in the sfc_data.nc)
        if ( nrecord == 17 ) then
           if ( my_proc_id == io_proc ) then
              allocate(dat43(ix, iy, 1, 1))
              call get_var_data(trim(infile_sfc), 'f10m', ix, iy, 1, 1, dat43)
           endif
        endif

        !-----------------------------
        !---7.18 output 3d
        if ( nrecord == 3 .or. nrecord == 4 .or. nrecord == 5 .or. &
             nrecord == 8 .or. nrecord == 9 .or. nrecord ==11 )then
           call mpi_barrier(comm,ierr)
           kz=nz
           if ( nrecord ==  9 .or. nrecord == 11 ) then
              kz=nz+1
           endif
           !--- map fv3 grid to rot-ll grid: ingrid-->dstgrid
           !call cpu_time(cputime2)
           !write(*,'(a,i3,f)')' --- read rot-ll grid for 1 record ', nrecord, cputime2

           if ( nprocs == 1 ) then  !--no mpi
              !----only 1-core
              allocate(dat41(nx,ny,kz,1), dat42(nx,ny,kz,1))
              if ( nd > 1 ) then
                 read(flid_in)dat42
                 do k = 1, kz
                     dat41(:,:,k,1)=dat42(:,:,kz-k+1,1)
                 enddo
                 dat42=-999999.
              else
                 dat41=-999999.
              endif
              call combine_grids_for_remap(ix,iy,kz,1,dat43,nx,ny,kz,1,dat41,gwt%gwt_t,dat42)

              !--- output
              !write(*,'(a,i2.2,a,200f)')'=== record',nrecord,': ', dat42(int(nx/2),int(ny/2),:,1)
              write(flid_out) (((dat42(i,j,k,1),i=1,nx),j=1,ny),k=kz,1,-1)
              deallocate(dat41, dat42, dat43)
           else
              !----mpi: 0 is for IO, >0 is for computing

              !---when nd>1, get previous data at 0, and then send to other cores
              if ( nd > 1 ) then
                 if ( my_proc_id == io_proc ) then
                    allocate(dat4(nx,ny,kz,1), dat42(nx,ny,kz,1))
                    read(flid_in)dat42
                    do k = 1, kz
                        dat4(:,:,k,1)=dat42(:,:,kz-k+1,1)
                    enddo
                    deallocate(dat42)
                    nm=max(1,int((kz+nprocs-1)/nprocs))
                    do k = 0, nprocs-1
                       ks=k*nm+1
                       ke=k*nm+nm
                       if ( ke > iz ) ke=iz
                       if ( ks >= 1 .and. ks <= iz .and. ke >= 1 .and. ke <= iz ) then
                          allocate(dat42(nx, ny, ke-ks+1,1))
                          dat42(:,:,:,1)=dat4(:,:,ks:ke,1)
                          if ( k /= io_proc ) then
                             call mpi_send(dat42(1,1,1,1),size(dat42),mpi_real, k, 4000+ks, comm, ierr)
                          else
                             allocate(dat41(nx, ny, ke-ks+1,1))
                             dat41=dat42
                          endif
                          deallocate(dat42)
                       endif
                    enddo
                    deallocate(dat4)
                 else  !if ( my_proc_id == io_proc ) then
                    !---receive dat42 for each core
                    nm=max(1,int((kz+nprocs-1)/nprocs))
                    ks=my_proc_id*nm+1
                    ke=my_proc_id*nm+nm
                    if ( ke > kz ) ke=kz
                    if ( ks >= 1 .and. ks <= kz .and. ke >= 1 .and. ke <= kz ) then
                       allocate(dat41(nx, ny, ke-ks+1,1))
                       call mpi_recv(dat41(1,1,1,1), size(dat41), mpi_real, io_proc, 4000+ks, comm, status, ierr)
                    endif
                 endif  !if ( my_proc_id == io_proc ) then
              else  !if ( nd > 1 ) then
                 nm=max(1,int((kz+nprocs-1)/nprocs))
                 ks=my_proc_id*nm+1
                 ke=my_proc_id*nm+nm
                 if ( ke > kz ) ke=kz
                 if ( ks >= 1 .and. ks <= kz .and. ke >= 1 .and. ke <= kz ) then
                    allocate(dat41(nx, ny, ke-ks+1,1))
                    dat41=-999999.
                 endif
              endif  !if ( nd > 1 ) then

              !---combine dat43+dat41 --> dat42
              !call mpi_barrier(comm,ierr)
              nm=max(1,int((kz+nprocs-1)/nprocs))
              ks=my_proc_id*nm+1
              ke=my_proc_id*nm+nm
              if ( ke > kz ) ke=kz
              if ( ks >= 1 .and. ks <= kz .and. ke >= 1 .and. ke <= kz ) then
                 allocate(dat42(nx, ny, ke-ks+1,1))
                 dat42=-999999.
                 call combine_grids_for_remap(ix,iy,ke-ks+1,1,dat43,nx,ny,ke-ks+1,1,dat41,gwt%gwt_t,dat42)
                 if ( my_proc_id /= io_proc ) then
                    call mpi_send(dat42(1,1,1,1),size(dat42),mpi_real, io_proc, 5000+ks, comm, ierr)
                    deallocate(dat42)
                 endif
                 deallocate(dat41, dat43)
              endif

              !---collect dat43 to io_proc, and output
              !call mpi_barrier(comm,ierr)
              if ( my_proc_id == io_proc ) then
                 allocate(dat43(nx,ny,kz,1))
                 nm=max(1,int((kz+nprocs-1)/nprocs))
                 do k = 0, nprocs-1
                    ks=k*nm+1
                    ke=k*nm+nm
                    if ( ke > kz ) ke=kz
                    if ( ks >= 1 .and. ks <= kz .and. ke >= 1 .and. ke <= kz ) then
                       if ( k /= io_proc ) then
                          allocate(dat41(nx, ny, ke-ks+1,1))
                          call mpi_recv(dat41(1,1,1,1), size(dat41), mpi_real, k, 5000+ks, comm, status, ierr)
                          dat43(:,:,ks:ke,1)=dat41(:,:,1:ke-ks+1,1)
                          deallocate(dat41)
                       else
                          dat43(:,:,ks:ke,1)=dat42(:,:,1:ke-ks+1,1)
                          deallocate(dat42)
                       endif
                    endif
                 enddo
                 !write(*,'(a,3i5,100f12.3)')'===w34 ', nx, ny, kz, (dat43(10,10,k,1),k=kz,1,-1)
                 write(flid_out) (((dat43(i,j,k,1),i=1,nx),j=1,ny),k=kz,1,-1)
                 deallocate(dat43)
              endif  !if ( my_proc_id == io_proc ) then
           endif  ! if ( nprocs == 1 ) then  !--no mpi
        else if ( nrecord == 6 ) then  !---u,v
           kz=nz
           !---convert u,v from fv3grid to earth
           nm=max(1,int((kz+nprocs-1)/nprocs))
           ks=my_proc_id*nm+1; ke=my_proc_id*nm+nm
           if ( ke > kz ) ke=kz
           if ( ks >= 1 .and. ks <= kz .and. ke >= 1 .and. ke <= kz ) then
              allocate(u(ix,iy,ke-ks+1,1), v(ix,iy,ke-ks+1,1))

              ! KGao
              !allocate(dat2(ix, iy+1), dat21(ix+1, iy))
              allocate(dat2(ix, iy), dat21(ix, iy))
              do k = 1, ke-ks+1

                 ! KGao 
                 !call fv3uv2earth(ix, iy, dat42(:,:,k,1), dat43(:,:,k,1), cangu, sangu, cangv, sangv, dat2, dat21)
                 !---destage: C-/D- grid to A-grid
                 !u(:,:,k,1) = (dat2 (:,1:iy)+dat2 (:,2:iy+1))/2.0
                 !v(:,:,k,1) = (dat21(1:ix,:)+dat21(2:ix+1,:))/2.0
                 u(:,:,k,1) = dat42 (:,:,k,1)
                 v(:,:,k,1) = dat43 (:,:,k,1)
              enddo
              deallocate(dat42, dat43, dat2, dat21, cangu, sangu, cangv, sangv)
           endif
           !write(*,*)'===w12 dat42,dat43 to earth wind u,v'

           !--- loop u,v
           do nv = 1, 2
              !--- get outer domain u,v
              if ( nd > 1 ) then
                 if ( my_proc_id == io_proc ) then
                    allocate(dat4(nx,ny,kz,1), dat42(nx,ny,kz,1))
                    read(flid_in)dat42
                    do k = 1, kz
                       dat4(:,:,k,1)=dat42(:,:,kz-k+1,1)
                    enddo
                    deallocate(dat42)
                    nm=max(1,int((kz+nprocs-1)/nprocs))
                    do k = 0, nprocs-1
                       ks=k*nm+1
                       ke=k*nm+nm
                       if ( ke > iz ) ke=iz
                       if ( ks >= 1 .and. ks <= iz .and. ke >= 1 .and. ke <= iz ) then
                          allocate(dat42(nx, ny, ke-ks+1,1))
                          dat42(:,:,:,1)=dat4(:,:,ks:ke,1)
                          if ( k /= io_proc ) then
                             call mpi_send(dat42(1,1,1,1),size(dat42),mpi_real, k, 300*nv+ks, comm, ierr)
                          else
                             allocate(dat41(nx, ny, ke-ks+1,1))
                             dat41=dat42
                          endif
                          deallocate(dat42)
                       endif
                    enddo
                    deallocate(dat4)
                 else  !if ( my_proc_id == io_proc ) then
                    !---receive dat42 for each core
                    nm=max(1,int((kz+nprocs-1)/nprocs))
                    ks=my_proc_id*nm+1
                    ke=my_proc_id*nm+nm
                    if ( ke > kz ) ke=kz
                    if ( ks >= 1 .and. ks <= kz .and. ke >= 1 .and. ke <= kz ) then
                       allocate(dat41(nx, ny, ke-ks+1,1))
                       call mpi_recv(dat41(1,1,1,1), size(dat41), mpi_real, io_proc, 300*nv+ks, comm, status, ierr)
                    endif
                 endif  !if ( my_proc_id == io_proc ) then
              else
                 nm=max(1,int((kz+nprocs-1)/nprocs))
                 ks=my_proc_id*nm+1
                 ke=my_proc_id*nm+nm
                 if ( ke > kz ) ke=kz
                 if ( ks >= 1 .and. ks <= kz .and. ke >= 1 .and. ke <= kz ) then
                    allocate(dat41(nx, ny, ke-ks+1,1))
                    dat41=-999999.
                 endif
              endif  !if ( nd > 1 ) then
              !write(*,*)'===w13 got outer domain dat41'

              !--- map u/v to rot-ll grid: ingrid-->dstgrid
              nm=max(1,int((kz+nprocs-1)/nprocs))
              ks=my_proc_id*nm+1
              ke=my_proc_id*nm+nm
              if ( ke > kz ) ke=kz
              if ( ks >= 1 .and. ks <= kz .and. ke >= 1 .and. ke <= kz ) then
                 allocate(dat42(nx, ny, ke-ks+1,1))
                 dat42=-999999.
                 if (nv==1) then
                    call combine_grids_for_remap(ix,iy,ke-ks+1,1,u,nx,ny,ke-ks+1,1,dat41,gwt%gwt_t,dat42)
                 else if (nv==2) then
                    call combine_grids_for_remap(ix,iy,ke-ks+1,1,v,nx,ny,ke-ks+1,1,dat41,gwt%gwt_t,dat42)
                 endif
                 if ( my_proc_id /= io_proc ) then
                    call mpi_send(dat42(1,1,1,1),size(dat42),mpi_real, io_proc, 400*nv+ks, comm, ierr)
                    deallocate(dat42)
                 endif
                 deallocate(dat41)
              endif
              !write(*,*)'===w14 got dat42'

              !---collect dat43 to io_proc, and output
              !call mpi_barrier(comm,ierr)
              if ( my_proc_id == io_proc ) then
                 allocate(dat43(nx,ny,kz,1))
                 nm=max(1,int((kz+nprocs-1)/nprocs))
                 do k = 0, nprocs-1
                    ks=k*nm+1
                    ke=k*nm+nm
                    if ( ke > kz ) ke=kz
                    if ( ks >= 1 .and. ks <= kz .and. ke >= 1 .and. ke <= kz ) then
                       if ( k /= io_proc ) then
                          allocate(dat41(nx, ny, ke-ks+1,1))
                          call mpi_recv(dat41(1,1,1,1), size(dat41), mpi_real, k, 400*nv+ks, comm, status, ierr)
                          dat43(:,:,ks:ke,1)=dat41(:,:,1:ke-ks+1,1)
                          deallocate(dat41)
                       else
                          dat43(:,:,ks:ke,1)=dat42(:,:,1:ke-ks+1,1)
                          deallocate(dat42)
                       endif
                    endif
                 enddo
                 !write(*,'(a,3i5,100f12.3)')'===w51 ', nx, ny, kz, (dat43(10,10,k,1),k=kz,1,-1)
                 write(flid_out) (((dat43(i,j,k,1),i=1,nx),j=1,ny),k=kz,1,-1)
                 deallocate(dat43)
              endif  !if ( my_proc_id == io_proc ) then
              if (nv==1) then
                 deallocate(u)
              else if (nv==2) then
                 deallocate(v)
              endif
           enddo  !do nv = 1, 2
        else if ( nrecord ==12 .or. nrecord ==15 .or. nrecord ==16 .or. nrecord ==17 ) then
           kz=1
           if ( my_proc_id == io_proc ) then
              allocate(dat41(nx,ny,kz,1), dat42(nx,ny,kz,1))
              if ( nd > 1 ) then
                 read(flid_in)dat42
                 do k = 1, kz
                     dat41(:,:,k,1)=dat42(:,:,kz-k+1,1)
                 enddo
                 dat42=-999999.
              else
                 dat41=-999999.
                 dat42=-999999.
              endif
              call combine_grids_for_remap(ix,iy,kz,1,dat43,nx,ny,kz,1,dat41,gwt%gwt_t,dat42)
              write(flid_out) (((dat42(i,j,k,1),i=1,nx),j=1,ny),k=kz,1,-1)
              deallocate(dat41, dat42, dat43)
           endif  !if ( my_proc_id == io_proc ) then
        endif

     enddo do_out_var_loop !: for nrecord = 1, 17

     !-------------------------------------------------------------------------
     ! 8 --- clean ingrid gwt
     deallocate( ingrid%grid_lon, ingrid%grid_lat, ingrid%grid_lont, ingrid%grid_latt)
     deallocate( gwt%gwt_t, gwt%gwt_u, gwt%gwt_v )

  enddo do_nestdom_loop !: do nd = 1, ndom
  write(*,*)' === finished hafsvi_preproc ==='

  return
  end subroutine hafsvi_preproc_ic
