  subroutine hafsvi_create_nc(in_dir, zind_str, out_file)

  use constants
  use netcdf
  use module_mpi
  use var_type

  implicit none

  character (len=*), intent(in) :: in_dir, out_file
  character (len=2500) :: in_file, in_file_ctrl

  integer  :: io_proc, nrecord, ndims, nx, ny, nz, nz1, zind_str
  integer, dimension(nf90_max_var_dims) :: dims

  real, allocatable, dimension(:,:,:,:) :: var_in 
  real, allocatable, dimension(:,:,:) :: var_out 
  real, allocatable, dimension(:,:) :: var_out2
  real, allocatable, dimension(:) :: ak 
  real, allocatable, dimension(:) :: bk

  io_proc=nprocs-1

  if ( my_proc_id == io_proc ) then

     in_file=trim(in_dir)//'/gfs_data.tile7.nc'
     in_file_ctrl=trim(in_dir)//'/gfs_ctrl.nc'

     ! get nx,ny,nz from ic file
     call get_var_dim(trim(in_file), 'delp', ndims, dims)
     nz=dims(3)
     ny=dims(2)
     nx=dims(1)
     nz1 = nz-zind_str+1 

     ! create dims for out nc file
     call write_nc_dim(trim(out_file), 'lon', nx)
     call write_nc_dim(trim(out_file), 'lat', ny)
     call write_nc_dim(trim(out_file), 'lev', nz1)
     call write_nc_dim(trim(out_file), 'levp', nz1+1)
     call write_nc_dim(trim(out_file), 'nvcoord', 2)

     ! loops over 'u', 'v', 'zh', 'delp', 't', 'sphum', 'w', 'vcoord'
     do nrecord = 1, 8 

        ! -- 'u'
        if ( nrecord == 1 ) then
           allocate(var_in (nx, ny+1, nz, 1))
           allocate(var_out(nx, ny, nz1))
           call get_var_data(trim(in_file), 'u_s', nx, ny+1, nz, 1, var_in)
           var_out = 0.5*( var_in(:,1:ny,zind_str:nz,1) + var_in(:,2:ny+1,zind_str:nz,1) )
           call write_nc_real(trim(out_file), 'u', nx, ny, nz1, -1, 'lon', 'lat', 'lev', '-', var_out, 'm/s', 'u_on_a_grid')
           deallocate(var_in, var_out)

        ! -- 'v'
        else if ( nrecord == 2 ) then
           allocate(var_in (nx, ny+1, nz, 1))
           allocate(var_out(nx, ny, nz1))
           call get_var_data(trim(in_file), 'v_s', nx, ny+1, nz, 1, var_in)
           var_out = 0.5*( var_in(:,1:ny,zind_str:nz,1) + var_in(:,2:ny+1,zind_str:nz,1) )
           call write_nc_real(trim(out_file), 'v', nx, ny, nz1, -1, 'lon', 'lat', 'lev', '-', var_out, 'm/s', 'v_on_a_grid')
           deallocate(var_in, var_out)

        ! -- 'zh'
        else if ( nrecord == 3 ) then
           allocate(var_in (nx, ny, nz+1, 1))
           allocate(var_out(nx, ny, nz1+1))
           call get_var_data(trim(in_file), 'zh', nx, ny, nz+1, 1, var_in)
           var_out = var_in(:, :, zind_str:nz+1, 1)
           call write_nc_real(trim(out_file), 'zh', nx, ny, nz1+1, -1, 'lon', 'lat', 'levp', '-', var_out, 'm', 'zh')
           deallocate(var_in, var_out)

        ! -- 'delp'
        else if ( nrecord == 4 ) then
           allocate(var_in (nx, ny, nz, 1))
           allocate(var_out(nx, ny, nz1)) 
           call get_var_data(trim(in_file), 'delp', nx, ny, nz, 1, var_in)
           var_out = var_in(:, :, zind_str:nz, 1)
           call write_nc_real(trim(out_file), 'delp', nx, ny, nz1, -1, 'lon', 'lat', 'lev', '-', var_out, 'pa', 'delp')
           deallocate(var_in, var_out)

        ! -- 't'
        else if ( nrecord == 5 ) then
           allocate(var_in (nx, ny, nz, 1))
           allocate(var_out(nx, ny, nz1))
           call get_var_data(trim(in_file), 't', nx, ny, nz, 1, var_in)
           var_out = var_in(:, :, zind_str:nz, 1)
           call write_nc_real(trim(out_file), 't', nx, ny, nz1, -1, 'lon', 'lat', 'lev', '-', var_out, 'k', 'temp')
           deallocate(var_in, var_out)

        ! -- 'sphum'
        else if ( nrecord == 6 ) then
           allocate(var_in (nx, ny, nz, 1))
           allocate(var_out(nx, ny, nz1))
           call get_var_data(trim(in_file), 'sphum', nx, ny, nz, 1, var_in)
           var_out = var_in(:, :, zind_str:nz, 1)
           call write_nc_real(trim(out_file), 'sphum', nx, ny, nz1, -1, 'lon', 'lat', 'lev', '-', var_out, 'g/g', 'sphum')
           deallocate(var_in, var_out)

        ! -- 'w'
        else if ( nrecord == 7 ) then
           allocate(var_in (nx, ny, nz, 1))
           allocate(var_out(nx, ny, nz1))
           call get_var_data(trim(in_file), 'w', nx, ny, nz, 1, var_in)
           var_out = var_in(:, :, zind_str:nz, 1)
           call write_nc_real(trim(out_file), 'w', nx, ny, nz1, -1, 'lon', 'lat', 'lev', '-', var_out, 'm/s', 'w')
           deallocate(var_in, var_out)

        ! -- 'vcoord'
        else if ( nrecord == 8 ) then
           allocate(var_in(nz+1,  2, 1, 1))
           call get_var_data(trim(in_file_ctrl), 'vcoord', nz+1, 2, 1, 1, var_in)

           allocate(var_out2(nz1+1, 2))
           var_out2 = var_in(zind_str:nz+1, :, 1, 1)
           call write_nc_real(trim(out_file), 'vcoord', -1, -1, nz1+1, 2, '-', '-', 'levp', 'nvcoord', var_out2, 'scalar', 'ak_bk')
           deallocate(var_in, var_out2)

           !allocate(ak(nz+1))
           !allocate(bk(nz+1))
           !ak=var_in(zind_str:nz+1, 1, 1, 1)
           !bk=var_in(zind_str:nz+1, 2, 1, 1)
           !call write_nc_real(trim(out_file), 'ak', -1, -1, nz1+1, -1, '-', '-', 'levp', '-', ak, 'na', 'ak')
           !call write_nc_real(trim(out_file), 'bk', -1, -1, nz1+1, -1, '-', '-', 'levp', '-', bk, 'na', 'bk')
           !deallocate(var_in, ak, bk)

        endif

     !write(*,*) 'selected var shape is', shape(u_a)

     enddo

  endif
  return
  end subroutine 
