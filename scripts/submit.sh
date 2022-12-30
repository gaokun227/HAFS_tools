#! /bin/sh
#SBATCH --job-name=hafs_atm_vi_test
#SBATCH --output=/mnt/lfs1/HFIP/hfip-gfdl/Kun.Gao/hafs_tools/scripts/stdout/%x.o%j
#SBATCH --account=hfip-gfdl
##SBATCH --qos=batch
#SBATCH --qos=urgent
##SBATCH --nodes=1-1
#SBATCH --tasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --partition=xjet

export TOTAL_TASKS='24'
export NCTSK='24'
export OMP_THREADS='1'
export envir='prod'
export WHERE_AM_I='jet'

#export CDATE=2021092400
#export STORMID=18L

export CDATE=2022092000
export STORMID=07L

# what exhafs_atm_vi needs
export version=hafs_tools 
export HOMEhafs=/mnt/lfs1/HFIP/hfip-gfdl/Kun.Gao/${version}
export USHhafs=/mnt/lfs1/HFIP/hfip-gfdl/Kun.Gao/${version}/ush
export EXEChafs=/mnt/lfs1/HFIP/hfip-gfdl/Kun.Gao/${version}/sorc/hafs_tools.fd/exec/
#export PARMhafs=/mnt/lfs1/HFIP/hfip-gfdl/Kun.Gao/${version}/parm
export FIXhafs=/mnt/lfs1/HFIP/hfip-gfdl/Kun.Gao/${version}/fix

# not important
export COMhafsprior=/mnt/lfs1/HFIP/hfip-gfdl/Kun.Gao/hafstmp_test/hafsv0p3_bstest_h3db_vida/com/2021092400/18L/../../2021092318/18L

#export WORKhafs=/mnt/lfs1/HFIP/hfip-gfdl/Kun.Gao/hafstmp_test/hafsv0p3_bstest_h3db_vida/${CDATE}/${STORMID}
export WORKhafs=/mnt/lfs1/HFIP/hfip-gfdl/Kun.Gao/tshield_tile6/hafsv0p3_bstest_h3db_vida/${CDATE}/${STORMID}
export WORKhafs=/mnt/lfs1/HFIP/hfip-gfdl/Kun.Gao/tshield_atl_vi/hafsv0p3_bstest_h3db_vida/${CDATE}/${STORMID}

# env
source ${USHhafs}/hafs_pre_job.sh.inc
module use ${HOMEhafs}/modulefiles
module load modulefile.hafs.jet
#module load hafs.jet
module list
#source ${USHhafs}/hafs_runcmd.sh.inc

# Execute ex-script
${HOMEhafs}/scripts/exhafs_atm_vi.sh
