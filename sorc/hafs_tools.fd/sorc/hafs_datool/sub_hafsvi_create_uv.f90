  subroutine hafsvi_create_uv_on_agrid(in_file, ix, iy, iz, out_file)

  use constants
  use netcdf
  use module_mpi
  use var_type

  implicit none

  character (len=*), intent(in) :: out_file
  character (len=2500)   :: in_file
  integer  :: ix, iy, iz
  real, allocatable, dimension(:,:,:,:) :: u_s, v_s
  real, allocatable, dimension(:,:,:) u_a, v_a 

  io_proc=nprocs-1

  if ( my_proc_id == io_proc ) then
     allocate(u_s(ix,   iy+1, iz, 1))
     allocate(v_s(ix+1, iy,   iz, 1))

     allocate(u_a(ix, iy, iz))
     allocate(v_a(ix, iy, iz))

     call get_var_data(trim(infile), 'u_s', ix, iy+1, iz, 1, u_s)
     call get_var_data(trim(infile), 'v_s', ix, iy+1, iz, 1, v_s)

     ! in python
     !ua = 0.5*(us[:,:-1,:]+us[:,1:,:])
     !va = 0.5*(vs[:,:-1,:]+vs[:,1:,:])

     u_a = 0.5*( u_s(:,1:iy,:,1) + u_s(:,2:iy+1,:,1) )
     v_a = 0.5*( v_s(:,1:iy,:,1) + v_s(:,2:iy+1,:,1) )

     call write_nc_dim(trim(outfile), 'lon', ix)
     call write_nc_dim(trim(outfile), 'lat', iy)
     call write_nc_dim(trim(outfile), 'lev', iz)
     !call write_nc_real0d(trim(outfile), 'lat1', lat1, 'degree', 'latitude 1')
     call write_nc_real(trim(outfile), 'u', ix, iy, iz, -1, 'lon', 'lat', 'lev', '-', u_a, 'm/s', 'u_on_a_grid')
     call write_nc_real(trim(outfile), 'v', ix, iy, iz, -1, 'lon', 'lat', 'lev', '-', v_a, 'm/s', 'v_on_a_grid')

  endif
  return
  end subroutine hafsvi_preproc_ic
