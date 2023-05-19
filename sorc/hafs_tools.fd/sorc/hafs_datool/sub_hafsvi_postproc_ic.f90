!========================================================================================
  subroutine hafsvi_postproc_ic(in_file, in_date, out_dir, nestdoms)

!-----------------------------------------------------------------------------
! HAFS DA tool - hafsvi_postproc
! Yonghui Weng, 20220121
!
! This subroutine reads hafs_vi binary output file and merge it to hafs restart files.
! hafs_vi binary output:
!      WRITE(IUNIT) NX,NY,NZ,I360
!      WRITE(IUNIT) LON1,LAT1,LON2,LAT2,CENTRAL_LON,CENTRAL_LAT
!      WRITE(IUNIT) PMID1
!      WRITE(IUNIT) T1
!      WRITE(IUNIT) Q1
!      WRITE(IUNIT) U1
!      WRITE(IUNIT) V1
!      WRITE(IUNIT) DZDT
!      WRITE(IUNIT) Z1
!!     WRITE(IUNIT) GLON,GLAT
!      WRITE(IUNIT) HLON,HLAT,VLON,VLAT
!      WRITE(IUNIT) P1
!      WRITE(IUNIT) PD1
!      WRITE(IUNIT) ETA1
!      WRITE(IUNIT) ETA2
!
!      ALLOCATE ( T1(NX,NY,NZ),Q1(NX,NY,NZ) )
!      ALLOCATE ( U1(NX,NY,NZ),V1(NX,NY,NZ),DZDT(NX,NY,NZ) )
!      ALLOCATE ( Z1(NX,NY,NZ+1),P1(NX,NY,NZ+1) )
!      ALLOCATE ( GLON(NX,NY),GLAT(NX,NY) )
!      ALLOCATE ( PD1(NX,NY),ETA1(NZ+1),ETA2(NZ+1) )
!      ALLOCATE ( USCM(NX,NY),VSCM(NX,NY) )          ! Env. wind at new grids
!      ALLOCATE ( HLON(NX,NY),HLAT(NX,NY) )
!      ALLOCATE ( VLON(NX,NY),VLAT(NX,NY) )
!      ALLOCATE ( PMID1(NX,NY,NZ),ZMID1(NX,NY,NZ) )

!-----------------------------------------------------------------------------

  use constants
  use netcdf
  use module_mpi
  use var_type

  implicit none

  character (len=*), intent(in) :: in_file,  & ! The VI output binary file on 30x30degree
                                   in_date,  & ! HAFS_restart file date, like 20200825.120000
                                   out_dir     ! HAFS_restart_folder, which holds grid_spec.nc, fv_core.res.tile1.nc,
                                               ! fv_srf_wnd.res.tile1.nc, fv_tracer.res.tile1.nc, phy_data.nc, sfc_data.nc
  integer, intent(in)           :: nestdoms

  type(grid2d_info)             :: ingrid   ! hafs restart grid
  type(grid2d_info)             :: dstgrid  ! rot-ll grid for output

!----for hafs restart
  integer  :: ix, iy, iz, kz, ndom, nd
  character(len=2500)     :: ncfile
  character (len=2500)   :: ncfile_fvcore, ncfile_core, ncfile_uv, ncfile_tracer, ncfile_phy, ncfile_sfc, ncfile_grid, ncfile_grid2, ncfile_atmos, ncfile_oro
  character (len=50) :: nestfl, tilefl, tempfl
                        ! grid_spec.nc : grid_spec.nest02.tile2.nc
                        ! fv_core.res.tile1.nc : fv_core.res.nest02.tile2.nc
                        ! phy_data.nc  : phy_data.nest02.tile2.nc


!----for hafsvi
  integer  :: nx, ny, nz, i360, filetype  ! filetype: 1=bin, 2=nc
  real     :: lon1,lat1,lon2,lat2,cen_lat,cen_lon,dlat,dlon
  real, allocatable, dimension(:,:) :: hlon, hlat, vlon, vlat

  integer  :: i, j, k, n, flid, ncid, ndims, nrecord, iunit
  real, allocatable, dimension(:,:,:,:) :: dat4, dat41, dat42, dat43, dat44, phis1, phis2, sfcp1, sfcp2, u1, v1, u, v
  real, allocatable, dimension(:,:,:)   :: dat3, dat31
  real, allocatable, dimension(:,:)     :: dat2, dat21
  real, allocatable, dimension(:)       :: dat1
  real     :: ptop
  logical  :: file_exist

  real, allocatable, dimension(:,:)     :: cangu, sangu, cangv, sangv

  integer :: io_proc, nm, ks, ke, nv

!------------------------------------------------------------------------------
! 1 --- arg process
! 1.1 --- i/o processor
  io_proc=nprocs-1  !last processor as I/O
  !io_proc=0

! 1.2 --- ndom
  ndom=nestdoms+1

!------------------------------------------------------------------------------
! 2 --- input grid info

  !-----------------------------
  !---2.1 get input grid info from binary file
  iunit=36
  open(iunit, file=trim(in_file), form='unformatted')
  read(iunit) nx, ny, nz!, i360
  write(*,'(a,4i5)')'===w40 nx, ny, nz, i360 = ',nx, ny, nz, i360
  read(iunit) lon1,lat1,lon2,lat2,cen_lon,cen_lat
  write(*,'(a,6f10.3)')'lon1,lat1,lon2,lat2,cen_lon,cen_lat =', lon1,lat1,lon2,lat2,cen_lon,cen_lat

  !!---add to test vortex-replacement
  !tc%vortexrep=1
  !tc%lat=cen_lat
  !tc%lon=cen_lon
  !!---add to test vortex-replacement

  do i = 1, 7
     read(iunit)
  enddo
  allocate(hlon(nx,ny), hlat(nx,ny), vlon(nx,ny), vlat(nx,ny))
  read(iunit)hlon, hlat, vlon, vlat
  close(iunit)
! New change: Jul 2022 JH Shin
! -------------------------------------------------------------

! For WPAC storms, CONVERT positive values of western hemisphere within the VI domain
! into negative value, IF the portion of western hemisphere (e.g., the eastern side of IDL)
! is included in VI domain, because VI is completed
  if ( cen_lon > 0. ) then
     where ( hlon > 180. ) hlon=hlon-360.
     where ( vlon > 180. ) vlon=vlon-360.
  endif

! For CPAC storms located near the east of international date line
! CONVERT negative values of eastern hemisphere within the VI domain into positive value
! because VI is done
  if ( cen_lon < -140. )then
     where ( hlon <= -180. ) hlon=hlon+360.
     where ( vlon <= -180. ) vlon=vlon+360.
  endif

! New change: Jul 2022 JH Shin
! -------------------------------------------------------------


  if (my_proc_id==0) then
     write(*, '(a,8f10.3)')' hlon,hlat(1,1; nx,1; nx,ny; 1,ny) =', &
                        hlon(1,1), hlat(1,1), hlon(nx,1), hlat(nx,1), hlon(nx,ny), hlat(nx,ny), hlon(1,ny), hlat(1,ny)
     write(*, '(a,8f10.3)')' vlon,vlat(1,1; nx,1; nx,ny; 1,ny) =', &
                        vlon(1,1), vlat(1,1), vlon(nx,1), vlat(nx,1), vlon(nx,ny), vlat(nx,ny), vlon(1,ny), vlat(1,ny)
  endif

  !-----------------------------
  !---2.2 define input rot-ll grids
  ingrid%grid_x = nx
  ingrid%grid_y = ny
  ingrid%ntime  = 1
  ingrid%grid_xt = nx
  ingrid%grid_yt = ny
  allocate(ingrid%grid_lon (ingrid%grid_x,ingrid%grid_y))
  allocate(ingrid%grid_lont(ingrid%grid_x,ingrid%grid_y))
  ingrid%grid_lon  = hlon
  ingrid%grid_lont = vlon
  allocate(ingrid%grid_lat (ingrid%grid_x,ingrid%grid_y))
  allocate(ingrid%grid_latt(ingrid%grid_x,ingrid%grid_y))
  ingrid%grid_lat  = hlat
  ingrid%grid_latt = vlat

!------------------------------------------------------------------------------
! --- domain loop: from inner domain to outer domain, so the number is from max to 1
  do_nestdom_loop: do nd = ndom, ndom

     !-------------------------------------------------------------------------
     ! 3 --- input file
     !       nestfl, tilefl: ncfile_core, ncfile_tracer, ncfile_grid, ncfile_atmos, ncfile_oro
     ! KGao - name changes
     !if ( nd == 0 ) then
     ncfile_grid=trim(out_dir)//'/grid_spec.nc'
     ncfile_core=trim(out_dir)//'/gfs_data.nc'
     ncfile_uv=trim(out_dir)//'/gfs_uv_agrid.nc'
     ncfile_tracer=trim(out_dir)//'/gfs_data.nc'
     !endif

     !-------------------------------------------------------------------------
     ! 4 --- input grid info
     !---4.1 read from grid file grid_spec.nc:

     ! KGao
     !inquire(file=ncfile_grid2, exist=file_exist)
     !if ( file_exist ) ncfile_grid = ncfile_grid2

     if ( debug_level > 10 ) write(*,'(a)')' --- read grid info from '//trim(ncfile_grid)
     call rd_grid_spec_data(trim(ncfile_grid), dstgrid)
     ix=dstgrid%grid_xt
     iy=dstgrid%grid_yt
     if ( debug_level > 10 ) then
        write(*,'(a,i,1x,i,1x,f,1x,f,1x,f,1x,f)')' --- dstgrid info: ', ix, iy, &
              dstgrid%grid_lon(int(ix/2), int(iy/2)), dstgrid%grid_lat(int(ix/2), int(iy/2)), &
              dstgrid%grid_lont(int(ix/2), int(iy/2)), dstgrid%grid_latt(int(ix/2), int(iy/2))
     endif

     !-----------------------------
     !---4.2 call FV3-grid cos and sin

     ! KGao
     !allocate( cangu(ix,iy+1),sangu(ix,iy+1),cangv(ix+1,iy),sangv(ix+1,iy) )
     !call cal_uv_coeff_fv3(ix, iy, dstgrid%grid_lat, dstgrid%grid_lon, cangu, sangu, cangv, sangv)

     !-----------------------------
     !---4.3 calculate output-grid in input-grid's positions (xin, yin), and each grid's weight to dst

     call cal_src_dst_grid_weight(ingrid, dstgrid)

     !-------------------------------------------------------------------------
     ! 5 --- process record one-by-one
     do_record_loop: do nrecord = 1, 14

        !write(*,*)' nrecord', nrecord

        if ( my_proc_id == io_proc ) open(iunit, file=trim(in_file), form='unformatted')

        !-----------------------------
        !---5.1 read data and derive out the var for restart
        iz=-99

        if ( nrecord == 1 .or. nrecord == 2 .or. nrecord == 3 .or. nrecord == 10 .or. &
             nrecord == 12 .or. nrecord == 13 .or. nrecord == 14 ) then
           !---ignore these records
           !---record 1 : nx, ny, nz, i360
           !---record 2 : lon1,lat1,lon2,lat2,cen_lon,cen_lat
           !---record 3 : pmid1(nx,ny,nz): pressure on full level
           !---                ignore, we use p1 to derive delp.
           !---record 10: hlon, hlat, vlon, vlat
           !---record 12: pd1,  PD1(NX,NY): surface pressure
           !---record 13: eta1, ETA1(NZ+1)
           !---record 14: eta2, ETA2(NZ+1)
           if ( my_proc_id == io_proc ) then
              if ( nrecord == 12 ) then
                 allocate(dat2(nx,ny))
                 read(iunit)dat2
                 deallocate(dat2)
              elseif ( nrecord == 13 .or. nrecord == 14 ) then
                 allocate(dat1(nz+1))
                 read(iunit)dat1
                 deallocate(dat1)
              else
                 read(iunit)
              endif
           endif
        endif

        !  ALLOCATE ( T1(NX,NY,NZ),Q1(NX,NY,NZ) )
        !  ALLOCATE ( U1(NX,NY,NZ),V1(NX,NY,NZ),DZDT(NX,NY,NZ) )
        !  ALLOCATE ( Z1(NX,NY,NZ+1),P1(NX,NY,NZ+1) )
        if ( nrecord == 6 ) then   !u,v - 6,7
           !---record 6 : U1
           !---record 7 : V1
           iz=nz
           if ( my_proc_id == io_proc ) then
              nm=max(1,int((iz+nprocs-1)/nprocs))
              do nv = 1, 2
                 !---get data
                 allocate(dat3(nx,ny,iz), dat4(nx,ny,iz,1))
                 read(iunit) dat3
                 do k = 1, nz
                    dat4(:,:,nz-k+1,1)=dat3(:,:,k)
                 enddo
                 deallocate(dat3)

                 !---send
                 if ( nprocs > 1 ) then
                    do k = 0, nprocs-1
                       ks=k*nm+1
                       ke=k*nm+nm            !k-end
                       if ( ke > iz ) ke=iz
                       if ( ks >= 1 .and. ks <= iz .and. ke >= 1 .and. ke <= iz ) then
                          allocate(dat41(nx, ny, ke-ks+1,1))
                          dat41(:,:,:,1)=dat4(:,:,ks:ke,1)
                          if ( k /= io_proc ) then
                             call mpi_send(dat41(1,1,1,1), size(dat41), mpi_real, k, 200*nv+ks, comm, ierr)
                          else
                             if ( nv == 1 ) then
                                allocate(dat42(nx, ny, ke-ks+1,1))
                                dat42=dat41
                             else if ( nv == 2 ) then
                                allocate(dat43(nx, ny, ke-ks+1,1))
                                dat43=dat41
                             endif
                          endif
                          deallocate(dat41)
                       endif
                    enddo
                 else  !if ( nprocs > 1 ) then
                    if ( nv == 1 ) then
                       allocate(dat42(nx, ny, iz, 1))
                       dat42=dat4
                    else if ( nv == 2 ) then
                       allocate(dat43(nx, ny, iz, 1))
                       dat43=dat4
                    endif
                 endif  !if ( nprocs > 1 ) then
                 deallocate(dat4)
              enddo  !do nv = 1, 2
           else    !if ( my_proc_id == io_proc ) then
              !---receive dat42 dat43
              if ( nprocs > 1 ) then
                 nm=max(1,int((iz+nprocs-1)/nprocs))
                 ks=min(iz,my_proc_id*nm+1)
                 ke=min(iz,my_proc_id*nm+nm)
                 if ( ks >= 1 .and. ks <= iz .and. ke >= 1 .and. ke <= iz ) then
                    allocate(dat42(nx, ny, ke-ks+1,1),dat43(nx, ny, ke-ks+1,1))
                    call mpi_recv(dat42(1,1,1,1), size(dat42), mpi_real, io_proc, 200*1+ks, comm, status, ierr)
                    call mpi_recv(dat43(1,1,1,1), size(dat43), mpi_real, io_proc, 200*2+ks, comm, status, ierr)
                 endif
              endif
           endif  !if ( my_proc_id == io_proc ) then
        endif  !if ( nrecord == 6 ) then   !u,v - 6,7

        if ( nrecord == 4 .or. nrecord == 5 .or. nrecord == 8 .or. &
             nrecord == 9 .or. nrecord == 11 ) then
           !---record 4 : t1-->T
           !---record 5 : Q1
           !---record 8 : DZDT
           !---record 9 : z1 --> DZ
           !---record 11: p1-->delp, p1(nx,ny,nz+1): (((p1(i,j,k),i=1,nx),j=1,ny),k=nz+1,1,-1)
           !---           p1-->ps
           iz=nz
           if ( my_proc_id == io_proc ) then
              !---get data
              if ( nrecord == 9 .or. nrecord == 11 ) then
                 allocate(dat3(nx,ny,iz+1))
              else
                 allocate(dat3(nx,ny,iz))
              endif
              read(iunit) dat3

              allocate(dat41(nx,ny,iz,1))
              if ( nrecord == 9 .or. nrecord == 11 ) then  ! z1 to dz; p1 to delp
                 !---back pressure to delp on fv_core.res.tile1.nc
                 do k = 1, nz
                    dat41(:,:,nz-k+1,1)=dat3(:,:,k)-dat3(:,:,k+1)
                 enddo
                 !---phis
                 if ( nrecord == 9 ) then
                    allocate(phis1(nx,ny,1,1))
                    phis1(:,:,1,1)=dat3(:,:,1)*g
                 endif
              else
                 do k = 1, iz
                    dat41(:,:,iz-k+1,1)=dat3(:,:,k)
                 enddo
              endif
              deallocate(dat3)

              !---send data to other cores
              if ( nprocs > 1 ) then
                 nm=max(1,int((iz+nprocs-1)/nprocs))  !devide iz to each processor
                 do k = 0, nprocs-1
                    ks=k*nm+1              !k-start
                    ke=k*nm+nm            !k-end
                    if ( ke > iz ) ke=iz
                    if ( ks >= 1 .and. ks <= iz .and. ke >= 1 .and. ke <= iz ) then
                       allocate(dat42(nx, ny, ke-ks+1,1))
                       dat42(:,:,:,1)=dat41(:,:,ks:ke,1)
                       if ( k /= io_proc ) then
                          call mpi_send(dat42(1,1,1,1), size(dat42), mpi_real, k, 3000+ks, comm, ierr)
                       else
                          allocate(dat43(nx, ny, ke-ks+1,1))
                          dat43=dat42
                       endif
                       deallocate(dat42)
                    endif
                 enddo
              else  !if ( nprocs > 1 ) then
                 allocate(dat43(nx, ny, iz,1))
                 dat43=dat41
              endif
              deallocate(dat41)
           else ! if ( my_proc_id == io_proc ) then
              !---receive dat43
              nm=max(1,int((iz+nprocs-1)/nprocs))
              ks=min(iz,my_proc_id*nm+1)
              ke=min(iz,my_proc_id*nm+nm)
              if ( ks >= 1 .and. ks <= iz .and. ke >= 1 .and. ke <= iz ) then
                 allocate(dat43(nx, ny, ke-ks+1,1))
                 call mpi_recv(dat43(1,1,1,1), size(dat43), mpi_real, io_proc, 3000+ks, comm, status, ierr)
              endif
           endif
        endif   !if ( nrecord == 4 .or. nrecord == 5 .or. nrecord == 8 .or. &

        !-----------------------------
        !---5.2 merge hafs restart and update restart files
        !---    note: need to change nesting domain's filenames
        if ( nrecord == 6 ) then  !u and v
           !---get u,v
           iz=nz
           nm=max(1,int((iz+nprocs-1)/nprocs))
           do nv = 1, 2
              if ( my_proc_id == io_proc ) then

                 ! KGao
                 !allocate(dat4(ix+nv-1, iy+2-nv, iz, 1))  !u(ix, iy+1, iz, 1), v(ix+1, iy, iz, 1)
                 allocate(dat4(ix, iy, iz, 1)) 

                 if ( nv == 1 ) then
                    ! KGao
                    call get_var_data(trim(ncfile_uv), 'u', ix, iy, iz,1, dat4)
                 else if ( nv == 2 ) then
                    ! KGao
                    call get_var_data(trim(ncfile_uv), 'v', ix, iy, iz,1, dat4)
                 endif
                 !---send to other core
                 if ( nprocs > 1 ) then
                    do k = 0, nprocs-1
                       ks=k*nm+1
                       ke=k*nm+nm            !k-end
                       if ( ke > iz ) ke=iz
                       if ( ks >= 1 .and. ks <= iz .and. ke >= 1 .and. ke <= iz ) then
                          if ( k /= io_proc ) then
                             ! KGao
                             allocate(dat41(ix, iy, ke-ks+1, 1))
                             dat41(:,:,1:ke-ks+1,1)=dat4(:,:,ks:ke,1)
                             call mpi_send(dat41(1,1,1,1), size(dat41), mpi_real, k, 200*nv+ks, comm, ierr)
                             deallocate(dat41)
                          else
                             ! KGao
                             allocate(dat44(ix, iy, ke-ks+1, 1))
                             dat44(:,:,1:ke-ks+1,1)=dat4(:,:,ks:ke,1)
                          endif  !if ( k /= io_proc ) then
                       endif  !if ( ks >= 1 .and. ks <= iz .and. ke >= 1 .and. ke <= iz ) then
                    enddo  !do k = 0, nprocs-1
                 else
                    ! KGao
                    allocate(dat44(ix, iy, iz, 1))
                    dat44=dat4
                 endif  !if ( nprocs > 1 ) then
                 deallocate(dat4)
              else  !if ( my_proc_id == io_proc ) then
                 !---receive u,v
                 ks=min(iz,my_proc_id*nm+1)
                 ke=min(iz,my_proc_id*nm+nm)
                 if ( ks >= 1 .and. ks <= iz .and. ke >= 1 .and. ke <= iz ) then
                    ! KGao
                    allocate(dat44(ix, iy, ke-ks+1, 1))
                    ! KGao fix
                    !call mpi_recv(dat44, size(dat43), mpi_real, io_proc, 200*nv+ks, comm, status, ierr)
                    call mpi_recv(dat44, size(dat44), mpi_real, io_proc, 200*nv+ks, comm, status, ierr)
                 endif
              endif  !if ( my_proc_id == io_proc ) then

              ks=min(iz,my_proc_id*nm+1)
              ke=min(iz,my_proc_id*nm+nm)
              if ( nv == 1 ) then
                 ! KGao
                 allocate(u(ix, iy, ke-ks+1,1))
                 u(:,:,1:ke-ks+1,1)=dat44(:,:,1:ke-ks+1,1)
              else if ( nv == 2 ) then
                 ! KGao
                 allocate(v(ix, iy, ke-ks+1,1))
                 v(:,:,1:ke-ks+1,1)=dat44(:,:,1:ke-ks+1,1)
              endif
              deallocate(dat44)
           enddo  !do nv = 1, 2

           !---convert fv3grid to earth
           ks=min(iz,my_proc_id*nm+1)
           ke=min(iz,my_proc_id*nm+nm)

           ! KGao
           allocate(dat4 (ix, iy, ke-ks+1, 1), dat41(ix, iy, ke-ks+1, 1))
           !$omp parallel do &
           !$omp& private(k)
           do k = 1, ke-ks+1
              ! KGao
              !call fv3uv2earth(ix, iy, u(:,:,k,1), v(:,:,k,1), cangu, sangu, cangv, sangv, dat4(:,:,k,1), dat41(:,:,k,1))
              dat4(:,:,k,1)  = u(:,:,k,1)
              dat41(:,:,k,1) = v(:,:,k,1)
           enddo
           deallocate(u,v)

           !---merge

           ! KGao
           !allocate(u1(ix, iy+1, ke-ks+1, 1), v1(ix+1, iy, ke-ks+1, 1))
           !u1=0.; v1=0.
           !call combine_grids_for_remap(nx,ny,ke-ks+1,1,dat42,ix,iy+1,ke-ks+1,1,dat4,gwt%gwt_u,u1)
           !call combine_grids_for_remap(nx,ny,ke-ks+1,1,dat43,ix+1,iy,ke-ks+1,1,dat41,gwt%gwt_v,v1)
           allocate(u1(ix, iy, ke-ks+1, 1), v1(ix, iy, ke-ks+1, 1))
           u1=0.; v1=0.
           call combine_grids_for_remap(nx,ny,ke-ks+1,1,dat42,ix,iy,ke-ks+1,1,dat4,gwt%gwt_t,u1)
           call combine_grids_for_remap(nx,ny,ke-ks+1,1,dat43,ix,iy,ke-ks+1,1,dat41,gwt%gwt_t,v1)
           
           deallocate(dat42, dat43, dat4, dat41)

           ! KGao

           !---convert earth wind to fv3grid wind
           allocate(u(ix, iy+1, ke-ks+1, 1), v(ix+1, iy, ke-ks+1, 1)) ! not used; D-grid wind
           !u=-999999.; v=-99999999.;
           !$omp parallel do &
           !$omp& private(k)
           !do k = 1, ke-ks+1
           !   call earthuv2fv3(ix, iy, u1(:,:,k,1), v1(:,:,k,1), cangu, sangu, cangv, sangv, u(:,:,k,1), v(:,:,k,1))
           !enddo
           !deallocate(u1,v1,cangu, sangu, cangv, sangv)

           !---send and collect
           if ( nprocs == 1 ) then

              ! KGao
              !allocate(u1(ix, iy+1, iz, 1), v1(ix+1, iy, iz, 1))
              !allocate(u1(ix, iy, iz, 1), v1(ix, iy, iz, 1))
              !u1=u
              !v1=v
              deallocate(u,v)
           else
              nm=max(1,int((iz+nprocs-1)/nprocs))
              if ( my_proc_id /= io_proc ) then
                 call mpi_send(u(1,1,1,1),size(u),mpi_real, io_proc, 400*1+my_proc_id, comm, ierr)
                 call mpi_send(v(1,1,1,1),size(v),mpi_real, io_proc, 400*2+my_proc_id, comm, ierr)
                 deallocate(u,v)
              else
                 ! KGao
                 !allocate(u1(ix, iy+1, iz, 1), v1(ix+1, iy, iz, 1))
                 !allocate(u1(ix, iy, iz, 1), v1(ix, iy, iz, 1))

                 do k = 0, nprocs-1
                    ks=k*nm+1             !k-start
                    ke=k*nm+nm            !k-end
                    if ( ke > iz ) ke=iz
                    if ( ks >= 1 .and. ks <= iz .and. ke >= 1 .and. ke <= iz ) then
                       if ( k /= io_proc ) then

                          ! KGao
                          !allocate(dat41(ix, iy+1, ke-ks+1, 1), dat42(ix+1, iy, ke-ks+1, 1))
                          allocate(dat41(ix, iy, ke-ks+1, 1), dat42(ix, iy, ke-ks+1, 1))

                          call mpi_recv(dat41(1,1,1,1), size(dat41), mpi_real, k, 400*1+k, comm, status, ierr)
                          call mpi_recv(dat42(1,1,1,1), size(dat42), mpi_real, k, 400*2+k, comm, status, ierr)
                          u1(:,:,ks:ke,1)=dat41(:,:,1:ke-ks+1,1)
                          v1(:,:,ks:ke,1)=dat42(:,:,1:ke-ks+1,1)
                          deallocate(dat41,dat42)
                       else
                          u1(:,:,ks:ke,1)=u(:,:,1:ke-ks+1,1)
                          v1(:,:,ks:ke,1)=v(:,:,1:ke-ks+1,1)
                          deallocate(u,v)
                       endif
                    endif
                 enddo   !do k = 0, nprocs-1
              endif
           endif   !if ( nprocs > 1 ) then  !need recv

           !---output
           if ( my_proc_id == io_proc ) then
              ! KGao
              call update_hafs_restart(trim(ncfile_uv), 'u', ix, iy, iz, -1, u1)
              call update_hafs_restart(trim(ncfile_uv), 'v', ix, iy, iz, -1, v1)
              deallocate(u1, v1)
           endif
        elseif ( nrecord == 4 .or. nrecord == 5 .or. nrecord == 8 .or. &
                 nrecord == 9 .or. nrecord == 11 ) then
           iz=nz
           nm=max(1,int((iz+nprocs-1)/nprocs))
           if ( my_proc_id == io_proc ) then
              !---get restart data
              allocate(dat4(ix, iy, iz, 1))
              ! KGao: 'T' and 'W' to 't' and 'w')
              if ( nrecord == 4 ) call get_var_data(trim(ncfile_core), 't', ix, iy, iz,1, dat4)
              if ( nrecord == 5 ) call get_var_data(trim(ncfile_tracer), 'sphum', ix, iy, iz,1, dat4)
              if ( nrecord == 8 ) call get_var_data(trim(ncfile_core), 'w', ix, iy, iz,1, dat4)
              ! KGao
              !if ( nrecord == 9 ) call get_var_data(trim(ncfile_core), 'DZ', ix, iy, iz,1, dat4)
              if ( nrecord == 9 ) call get_var_data(trim(ncfile_core), 'delp', ix, iy, iz,1, dat4) ! use a random var 

              if ( nrecord == 11 ) call get_var_data(trim(ncfile_core), 'delp', ix, iy, iz,1, dat4)

              !---send data to other cores
              if ( nprocs > 1 ) then
                 nm=max(1,int((iz+nprocs-1)/nprocs))  !devide iz to each processor
                 do k = 0, nprocs-1
                    ks=k*nm+1              !k-start
                    ke=k*nm+nm            !k-end
                    if ( ke > iz ) ke=iz
                    if ( ks >= 1 .and. ks <= iz .and. ke >= 1 .and. ke <= iz ) then
                       allocate(dat41(ix, iy, ke-ks+1,1))
                       dat41(:,:,1:ke-ks+1,1)=dat4(:,:,ks:ke,1)
                       if ( k /= io_proc ) then
                          call mpi_send(dat41(1,1,1,1), size(dat41), mpi_real, k, 3000+ks, comm, ierr)
                       else
                          allocate(dat42(ix, iy, ke-ks+1,1))
                          dat42=dat41
                       endif
                       deallocate(dat41)
                    endif
                 enddo
              else  !if ( nprocs > 1 ) then
                 allocate(dat42(ix, iy, iz,1))
                 dat42=dat4
              endif
              deallocate(dat4)
           else ! if ( my_proc_id == io_proc ) then
              !---receive dat43
              ks=min(iz,my_proc_id*nm+1)
              ke=min(iz,my_proc_id*nm+nm)
              if ( ks >= 1 .and. ks <= iz .and. ke >= 1 .and. ke <= iz ) then
                 allocate(dat42(ix, iy, ke-ks+1,1))
                 call mpi_recv(dat42(1,1,1,1), size(dat42), mpi_real, io_proc, 3000+ks, comm, status, ierr)
              endif
           endif

           !---merge
           ks=min(iz,my_proc_id*nm+1)
           ke=min(iz,my_proc_id*nm+nm)
           allocate(dat4(ix, iy, ke-ks+1, 1))
           call combine_grids_for_remap(nx,ny,ke-ks+1,1,dat43,ix,iy,ke-ks+1,1,dat42,gwt%gwt_t,dat4)
           deallocate(dat43, dat42)

           !---collect data to io_proc
           if ( nprocs == 1 ) then
              allocate(dat41(ix, iy, iz, 1))
              dat41=dat4
           else
              nm=max(1,int((iz+nprocs-1)/nprocs))  !devide iz to each processor
              if ( my_proc_id /= io_proc ) then
                 call mpi_send(dat4(1,1,1,1),size(dat4),mpi_real, io_proc, 600+ks, comm, ierr)
              else
                 allocate(dat41(ix, iy, iz, 1))
                 do k = 0, nprocs-1
                    ks=k*nm+1
                    ke=k*nm+nm
                    if ( ke > iz ) ke=iz
                    if ( ks >= 1 .and. ks <= iz .and. ke >= 1 .and. ke <= iz ) then
                       if ( k /= io_proc ) then
                          allocate(dat42(ix, iy, ke-ks+1,1))
                          call mpi_recv(dat42(1,1,1,1), size(dat42), mpi_real, k, 600+ks, comm, status, ierr)
                          dat41(:,:,ks:ke,1)=dat42(:,:,1:ke-ks+1,1)
                          deallocate(dat42)
                       else
                          dat41(:,:,ks:ke,1)=dat4(:,:,1:ke-ks+1,1)
                       endif
                    endif
                 enddo
              endif   !if ( my_proc_id /= io_proc ) then
           endif  !if ( nprocs == 1 ) then

           !---output
           if ( my_proc_id == io_proc ) then

              ! KGao: tx -> -1
              !---update restart
              if ( nrecord == 4 ) call update_hafs_restart(trim(ncfile_core), 't', ix, iy, iz, -1, dat41)
              if ( nrecord == 5 ) call update_hafs_restart(trim(ncfile_tracer), 'sphum', ix, iy, iz, -1, dat41)
              ! KGao - do not update w and DZ
              !if ( nrecord == 8 ) call update_hafs_restart(trim(ncfile_core), 'w', ix, iy, iz, 1, dat41)
              !if ( nrecord == 9 ) call update_hafs_restart(trim(ncfile_core), 'DZ', ix, iy, iz, 1, dat41)
              if ( nrecord ==11 ) call update_hafs_restart(trim(ncfile_core), 'delp', ix, iy, iz, -1, dat41)
              deallocate(dat41)

              !---2d phis
              ! KGao: do not update 'phis'
              !if ( nrecord == 9 .and. my_proc_id == io_proc ) then  !phis
              !   allocate(phis2(ix, iy, 1, 1), dat41(ix, iy, 1, 1))
              !   call get_var_data(trim(ncfile_core), 'phis', ix, iy, 1, 1, phis2)
              !   call combine_grids_for_remap(nx,ny,1,1,phis1,ix,iy,1,1,phis2,gwt%gwt_t,dat41)
              !   call update_hafs_restart(trim(ncfile_core), 'phis', ix, iy, 1, 1, dat41)
              !   deallocate(phis1, phis2, dat41)
              !endif

           endif  !if ( my_proc_id == io_proc ) then
           deallocate(dat4)
        endif

     enddo do_record_loop

     !-----------------------------
     ! 6 --- clean
     deallocate( dstgrid%grid_lon, dstgrid%grid_lat, dstgrid%grid_lont, dstgrid%grid_latt)
     deallocate( gwt%gwt_t, gwt%gwt_u, gwt%gwt_v )

  enddo do_nestdom_loop !: do nd = 1, ndom
  close(iunit)
  write(*,*)'--- hafsvi_postproc completed ---'

  return
  end subroutine hafsvi_postproc_ic
