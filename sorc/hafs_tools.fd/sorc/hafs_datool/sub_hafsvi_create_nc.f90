  subroutine hafsvi_create_nc(in_file,ix,iy,iz, out_file)

  use constants
  use netcdf
  use module_mpi
  use var_type

  implicit none

  character (len=*), intent(in) :: out_file
  character (len=2500)   :: infile

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

  !real, allocatable, dimension(:)       :: pfull, phalf
  real, allocatable, dimension(:,:)     :: cangu, sangu, cangv, sangv
  real    :: cputime1, cputime2, cputime3
  integer :: io_proc, nm, ks, ke, nv

  io_proc=nprocs-1

  if ( my_proc_id == io_proc ) then
     do nv = 1, 2
       if (nv==1) then
          allocate(dat4(ix, iy+1, iz,1))
          call get_var_data(trim(infile), 'u_s', ix, iy+1, iz, 1, dat4)
       else if (nv==2) then
          allocate(dat4(ix+1, iy, iz,1))
          call get_var_data(trim(infile), 'v_s', ix+1, iy, iz, 1, dat4)
          dat43=dat4
       endif
       deallocate(dat4)
     enddo  !do nv = 1, 2

          if ( my_proc_id == 0 ) then
        !write(flid_out) nx, ny, nz
        call write_nc_dim(trim(fl_out), 'nx', nx)
        call write_nc_dim(trim(fl_out), 'ny', ny)
        call write_nc_dim(trim(fl_out), 'nz', nz)
        call write_nc_dim(trim(fl_out), 'nz1', nz+1)
        call write_nc_real0d(trim(fl_out), 'lon1', lon1, 'degree', 'longtitude 1')
        call write_nc_real0d(trim(fl_out), 'lat1', lat1, 'degree', 'latitude 1')
        call write_nc_real0d(trim(fl_out), 'lon2', lon2, 'degree', 'longtitude 2')
        call write_nc_real0d(trim(fl_out), 'lat2', lat2, 'degree', 'latitude 1')
        call write_nc_real0d(trim(fl_out), 'cen_lon', cen_lon, 'degree', 'center of longtitude')
        call write_nc_real0d(trim(fl_out), 'cen_lat', cen_lat, 'degree', 'center of latitude')

        call write_nc_real(trim(fl_out), trim(varname), nx, ny, kz, -1, 'nx', 'ny', trim(nzc), '-', dat42, trim(units), trim(varname_long))
