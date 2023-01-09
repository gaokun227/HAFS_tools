#!/bin/bash
#SBATCH --output=/ncrc/home1/Kun.Gao/hafs_tools/scripts/stdout/%x.o%j
#SBATCH --job-name=hafs_atm_vi_test
#SBATCH --partition=batch
#SBATCH --account=gfdl_W
##SBATCH --qos=batch
#SBATCH --qos=debug
#SBATCH --ntasks=4
##SBATCH --tasks-per-node=4
##SBATCH --cpus-per-task=1
#SBATCH --time=03:00:00
#SBATCH --cluster=c4

##SBATCH --nodes=48

export TOTAL_TASKS='24'
export NCTSK='24'
export OMP_THREADS='1'
export envir='prod'
export WHERE_AM_I='gaea'

#export CDATE=2021092400
#export STORMID=18L

export CDATE=2022092000
export STORMID=07L

# what exhafs_atm_vi needs
export HOMEhafs=/ncrc/home1/Kun.Gao/hafs_tools/
export USHhafs=${HOMEhafs}/ush
export EXEChafs=${HOMEhafs}/sorc/hafs_tools.fd/exec/
export FIXhafs=${HOMEhafs}/fix

export WORKhafs=/lustre/f2/scratch/gfdl/Kun.Gao/tshield_tile6/${CDATE}/${STORMID}
export COMhafsprior=${WORKhafs}/../../2021092318/18L # random dir; not needed at this point

# env
#export OMP_NUM_THREADS=1
#export setenv OMP_STACKSIZE=1024M
#unlimit
source ${USHhafs}/hafs_pre_job.sh.inc
module use ${HOMEhafs}/modulefiles
module load modulefile.hafs.gaea
module list
#source ${USHhafs}/hafs_runcmd.sh.inc

# Execute ex-script
${HOMEhafs}/scripts/exhafs_atm_vi.sh
