#!/bin/bash

cwd=`pwd`
TOOLS_PATH=${cwd}/hafs_tools.fd

# Build hafs_datool
cd ${TOOLS_PATH}/sorc/hafs_datool/
./build_datool_CMake.bash

# Build hafs_vi
cd ${TOOLS_PATH}/sorc/hafs_vi/
./build_vi_CMake.bash

echo "DONE"
exit
