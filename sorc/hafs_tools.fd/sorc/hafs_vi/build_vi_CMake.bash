#!/bin/bash
#set -xeu
set -x # By Kai for C5

HOMEhafs="../../../../"
source ${HOMEhafs}/sorc/machine-setup.sh > /dev/null 2>&1

module use ${HOMEhafs}/modulefiles
module load modulefile.hafs.${target}_c5
module list

BuildDir=${HOMEhafs}/sorc/hafs_tools.fd/sorc/build
if [ -d ${BuildDir} ]; then 
   rm -rf ${BuildDir}
fi
mkdir ${BuildDir}
cd ${BuildDir}

#cmake ../hafs_vi  -DCMAKE_Fortran_COMPILER=ifort -DCMAKE_C_COMPILER=icc -DBUILD_TYPE=RELEASE
cmake ../hafs_vi

make VERBOSE=3
make install

# ./build_vi_CMake.bash > build.log 2>&1 &
