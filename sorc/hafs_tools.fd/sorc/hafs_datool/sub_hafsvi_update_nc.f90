  subroutine hafsvi_update_nc(in_dir, zind_str, dst_file)

  use constants
  use netcdf
  use module_mpi
  use var_type

  implicit none

  character (len=*), intent(in) :: in_dir, dst_file
  character (len=2500) :: in_file_before, in_file_after
  character (len=50) :: var_name 

  integer  :: io_proc, nrecord, nv, ndims, nx, ny, nz, nz1, zind_str
  integer, dimension(nf90_max_var_dims) :: dims

  real, allocatable, dimension(:,:,:,:) :: var_w, var_s
  real, allocatable, dimension(:,:,:)   :: var_w_new, var_s_new
  real, allocatable, dimension(:,:,:,:) :: var_before_vi, var_after_vi
  real, allocatable, dimension(:,:,:)   :: dvar, dvar_w, dvar_s 

  io_proc=nprocs-1

  if ( my_proc_id == io_proc ) then

     in_file_before=trim(in_dir)//'/gfs_data_for_vi.nc'
     in_file_after=trim(in_dir)//'/gfs_data_after_vi.nc'

     !write (*,*) zind_str
     !write (*,*) in_file_before
     !write (*,*) in_file_after

     ! get nx,ny,nz
     call get_var_dim(trim(in_file_before), 'delp', ndims, dims)
     nz1=dims(3) ! nz1 = nz - zind_str + 1
     ny=dims(2)
     nx=dims(1)
     nz = nz1+zind_str-1 

     ! loops over 'u'/'v', 'delp', 't', 'sphum'
     allocate(var_before_vi (nx, ny, nz1, 1))
     allocate(var_after_vi (nx, ny, nz1, 1))
     allocate(dvar (nx, ny, nz1))

     do nrecord = 1, 4

        ! -- 'u' or 'v'
        if ( nrecord == 1 ) then

           write(*,*) ' --- updating wind'
           allocate(var_w (nx+1, ny, nz, 1))
           allocate(var_s (nx, ny+1, nz, 1))
           allocate(var_w_new (nx+1, ny, nz))
           allocate(var_s_new (nx, ny+1, nz))
           allocate(dvar_w (nx+1, ny, nz1))
           allocate(dvar_s (nx, ny+1, nz1))

           do nv = 1, 2
             if (nv == 1) then
               call get_var_data(trim(in_file_before), 'u', nx, ny, nz1, 1, var_before_vi)
               call get_var_data(trim(in_file_after),  'u', nx, ny, nz1, 1, var_after_vi)
               call get_var_data(trim(dst_file), 'u_w', nx+1, ny,   nz,  1, var_w)
               call get_var_data(trim(dst_file), 'u_s', nx,   ny+1, nz,  1, var_s)
             else if (nv == 2) then
               call get_var_data(trim(in_file_before), 'v', nx, ny, nz1, 1, var_before_vi)
               call get_var_data(trim(in_file_after),  'v', nx, ny, nz1, 1, var_after_vi)
               call get_var_data(trim(dst_file), 'v_w', nx+1, ny,   nz,  1, var_w)
               call get_var_data(trim(dst_file), 'v_s', nx,   ny+1, nz,  1, var_s)
             endif

             dvar = var_after_vi(:,:,:,1) - var_before_vi(:,:,:,1)
             call map_dvar_to_cgrid(dvar, dvar_w, dvar_s, nx, ny, nz1)

             var_w_new = var_w(:,:,:,1)
             var_s_new = var_s(:,:,:,1)
             var_w_new(:,:,zind_str:nz) = var_w_new(:,:,zind_str:nz) + dvar_w
             var_s_new(:,:,zind_str:nz) = var_s_new(:,:,zind_str:nz) + dvar_s

             if (nv == 1) then
               call update_hafs_restart(trim(dst_file), 'u_w', nx+1, ny, nz, -1, var_w_new)
               call update_hafs_restart(trim(dst_file), 'u_s', nx, ny+1, nz, -1, var_s_new)
             else if (nv == 2) then
               call update_hafs_restart(trim(dst_file), 'v_w', nx+1, ny, nz, -1, var_w_new)
               call update_hafs_restart(trim(dst_file), 'v_s', nx, ny+1, nz, -1, var_s_new)
             endif
           enddo ! nv == 1 or 2

           deallocate(var_w, var_s, var_w_new, var_s_new)
           deallocate(dvar_w, dvar_s)

       ! 'delp', 't', 'sphum'    
       else ! nrecord == 2, 3, 4 

           if (nrecord == 2) then
             allocate(var_w (nx, ny, nz, 1))
             allocate(var_w_new (nx, ny, nz))
           endif 
           if (nrecord == 2) then
              var_name = 'delp'
           else if (nrecord == 3) then
              var_name = 't'
           else if (nrecord == 4) then
              var_name = 'sphum'
           endif
           write(*,*) ' --- updating ', var_name

           call get_var_data(trim(in_file_before), trim(var_name), nx, ny, nz1, 1, var_before_vi)
           call get_var_data(trim(in_file_after),  trim(var_name), nx, ny, nz1, 1, var_after_vi)
           call get_var_data(trim(dst_file),       trim(var_name), nx, ny, nz,  1, var_w)

           dvar = var_after_vi(:,:,:,1) - var_before_vi(:,:,:,1)
           var_w_new = var_w(:,:,:,1)
           var_w_new(:,:,zind_str:nz) = var_w_new(:,:,zind_str:nz) + dvar

           call update_hafs_restart(trim(dst_file), trim(var_name), nx, ny, nz, -1, var_w_new)

       endif ! if nrecord

     enddo ! do nrecord

  endif
  return
  end subroutine hafsvi_update_nc

  subroutine map_dvar_to_cgrid(dvar, dvar_w, dvar_s, ix, iy, iz)
  implicit none

  integer :: ix, iy, iz
  real, dimension(ix,   iy,   iz) :: dvar
  real, dimension(ix+1, iy,   iz) :: dvar_w
  real, dimension(ix,   iy+1, iz) :: dvar_s

  dvar_w(2:ix, :, :) = 0.5 * ( dvar(1:ix-1, :, :) + dvar(2:ix, :, :) )
  dvar_w(1,    :, :) = 0.5 *   dvar(1,      :, :) 
  dvar_w(ix+1, :, :) = 0.5 *   dvar(ix,     :, :) 

  dvar_s(:, 2:iy, :) = 0.5 * ( dvar(:, 1:iy-1, :) + dvar(:, 2:iy, :) )
  dvar_s(:, 1,    :) = 0.5 *   dvar(:, 1,      :)
  dvar_s(:, iy+1, :) = 0.5 *   dvar(:, iy,     :)

  return

  end subroutine map_dvar_to_cgrid 

