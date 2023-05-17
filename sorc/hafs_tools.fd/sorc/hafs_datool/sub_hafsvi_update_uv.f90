  subroutine hafsvi_update_uv_on_cgrid(afile_new, afile_old, cfile_new, cfile_old, ix, iy, iz)

  use constants
  use netcdf
  use module_mpi
  use var_type

  implicit none

  character (len=2500)   :: afile_new, afile_old, cfile_new, cfile_old 

  integer  :: ix, iy, iz, nv
  real, allocatable, dimension(:,:,:,:) a_new, a_old, cs_old, cw_old
  real, allocatable, dimension(:,:,:)   da, dcw, dcs, cs_new, cw_new 

  io_proc=nprocs-1

  if ( my_proc_id == io_proc ) then

     allocate(a_new(ix, iy, iz, 1))
     allocate(a_old(ix, iy, iz, 1))
     allocate(cw_old(ix+1, iy,   iz, 1))
     allocate(cs_old(ix,   iy+1, iz, 1))
     allocate(cw_old(ix+1, iy,   iz))
     allocate(cs_old(ix,   iy+1, iz))
     allocate(da(ix,    iy, iz))
     allocate(dcw(ix+1, iy, iz))
     allocate(dcs(ix,   iy+1, iz))

     do nv = 1, 2

        if (nv==1) then ! 'u' 
           call get_var_data(trim(afile_new), 'u', ix, iy, iz, 1, a_new)
           call get_var_data(trim(afile_old), 'u', ix, iy, iz, 1, a_old)
           call get_var_data(trim(cfile_old), 'u_w', ix+1, iy, iz, 1, cw_old)
           call get_var_data(trim(cfile_old), 'u_s', ix, iy+1, iz, 1, cs_old)
           da(:,:,:) = a_new(:,:,:,1) - a_old(:,:,:,1)
           dcw(2:ix,:,:) = 0.5 * ( da(1:ix-1,:,:) + da(2:ix,:,:) ) 
           dcw(1,   :,:) = 0.5 *   da(1,:,:)
           dcw(ix+1,:,:) = 0.5 *   da(ix,:,:)
           dcs(:,2:iy,:) = 0.5 * ( da(:,1:iy-1,:) + da(:,2:iy,:) )
           dcw(:,1,   :) = 0.5 *   da(:,1,     :)
           dcw(:,iy,  :) = 0.5 *   da(:,iy,    :)
           cs_new = cs_old + dcs
           cw_new = cw_old + dcw
           update_hafs_restart(cfile_new, 'u_w', ix+1, iy,   iz, -1, cw_new)
           update_hafs_restart(cfile_new, 'u_s', ix,   iy+1, iz, -1, cs_new)
        else if ( nv == 2 ) then
           call get_var_data(trim(afile_new), 'v', ix, iy, iz, 1, a_new)
           call get_var_data(trim(afile_old), 'v', ix, iy, iz, 1, a_old)
           call get_var_data(trim(cfile_old), 'v_w', ix+1, iy, iz, 1, cw_old)
           call get_var_data(trim(cfile_old), 'v_s', ix, iy+1, iz, 1, cs_old)
           da(:,:,:) = a_new(:,:,:,1) - a_old(:,:,:,1)
           dcw(2:ix,:,:) = 0.5 * ( da(1:ix-1,:,:) + da(2:ix,:,:) )
           dcw(1,   :,:) = 0.5 *   da(1,:,:)
           dcw(ix+1,:,:) = 0.5 *   da(ix,:,:)
           dcs(:,2:iy,:) = 0.5 * ( da(:,1:iy-1,:) + da(:,2:iy,:) )
           dcw(:,1,   :) = 0.5 *   da(:,1,     :)
           dcw(:,iy,  :) = 0.5 *   da(:,iy,    :)
           cs_new = cs_old + dcs
           cw_new = cw_old + dcw
           update_hafs_restart(cfile_new, 'v_w', ix+1, iy,   iz, -1, cw_new)
           update_hafs_restart(cfile_new, 'v_s', ix,   iy+1, iz, -1, cs_new)
        endif

     enddo

     deallocate(a_new)
     deallocate(a_old)
     deallocate(cs_new)
     deallocate(cs_old)
     deallocate(cw_new)
     deallocate(cw_old)
     deallocate(da)
     deallocate(dcs)
     deallocate(dcw)

  endif
  return
  end subroutine hafsvi_update_uv_on_cgrid
