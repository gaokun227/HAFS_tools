set(CMAKE_Fortran_COMPILER "/apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/bin/intel64/ifort")
set(CMAKE_Fortran_COMPILER_ARG1 "")
set(CMAKE_Fortran_COMPILER_ID "Intel")
set(CMAKE_Fortran_COMPILER_VERSION "18.0.5.20180823")
set(CMAKE_Fortran_COMPILER_WRAPPER "")
set(CMAKE_Fortran_PLATFORM_ID "Linux")
set(CMAKE_Fortran_SIMULATE_ID "")
set(CMAKE_Fortran_SIMULATE_VERSION "")




set(CMAKE_AR "/bin/ar")
set(CMAKE_Fortran_COMPILER_AR "")
set(CMAKE_RANLIB "/bin/ranlib")
set(CMAKE_Fortran_COMPILER_RANLIB "")
set(CMAKE_COMPILER_IS_GNUG77 )
set(CMAKE_Fortran_COMPILER_LOADED 1)
set(CMAKE_Fortran_COMPILER_WORKS TRUE)
set(CMAKE_Fortran_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

set(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_Fortran_COMPILER_ID_RUN 1)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;fpp;FPP;f77;F77;f90;F90;for;For;FOR;f95;F95)
set(CMAKE_Fortran_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_Fortran_LINKER_PREFERENCE 20)
if(UNIX)
  set(CMAKE_Fortran_OUTPUT_EXTENSION .o)
else()
  set(CMAKE_Fortran_OUTPUT_EXTENSION .obj)
endif()

# Save compiler ABI information.
set(CMAKE_Fortran_SIZEOF_DATA_PTR "8")
set(CMAKE_Fortran_COMPILER_ABI "ELF")
set(CMAKE_Fortran_LIBRARY_ARCHITECTURE "")

if(CMAKE_Fortran_SIZEOF_DATA_PTR AND NOT CMAKE_SIZEOF_VOID_P)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_Fortran_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_Fortran_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_Fortran_COMPILER_ABI}")
endif()

if(CMAKE_Fortran_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()





set(CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES "/lfs4/HFIP/hfv3gfs/nwprod/hpc-stack/libs/intel-18.0.5.274/impi-2018.4.274/netcdf/4.7.4/include;/lfs4/HFIP/hfv3gfs/nwprod/hpc-stack/libs/intel-18.0.5.274/impi-2018.4.274/hdf5/1.10.6/include;/lfs4/HFIP/hfv3gfs/nwprod/hpc-stack/libs/intel-18.0.5.274/szip/2.1.1/include;/lfs4/HFIP/hfv3gfs/nwprod/hpc-stack/libs/intel-18.0.5.274/impi-2018.4.274/esmf/8.2.1b04/include;/lfs4/HFIP/hfv3gfs/nwprod/hpc-stack/libs/intel-18.0.5.274/impi-2018.4.274/pio/2.5.2/include;/lfs4/HFIP/hfv3gfs/nwprod/hpc-stack/libs/intel-18.0.5.274/png/1.6.35/include;/lfs4/HFIP/hfv3gfs/nwprod/hpc-stack/libs/intel-18.0.5.274/zlib/1.2.11/include;/lfs4/HFIP/hfv3gfs/nwprod/hpc-stack/libs/intel-18.0.5.274/jasper/2.0.22/include;/apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/ipp/include;/apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/mkl/include;/apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/pstl/include;/apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/tbb/include;/apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/daal/include;/misc/apps/intel/compilers_and_libraries_2018.5.274/linux/compiler/include/intel64;/misc/apps/intel/compilers_and_libraries_2018.5.274/linux/compiler/include/icc;/misc/apps/intel/compilers_and_libraries_2018.5.274/linux/compiler/include;/usr/local/include;/usr/lib/gcc/x86_64-redhat-linux/4.8.5/include;/usr/include")
set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "ifport;ifcoremt;imf;svml;m;ipgo;irc;pthread;svml;c;gcc;gcc_s;irc_s;dl;c")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/lfs4/HFIP/hfv3gfs/nwprod/hpc-stack/libs/intel-18.0.5.274/impi-2018.4.274/esmf/8.2.1b04/lib;/lfs4/HFIP/hfv3gfs/nwprod/hpc-stack/libs/intel-18.0.5.274/impi-2018.4.274/pio/2.5.2/lib;/apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/compiler/lib/intel64_lin;/apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/mkl/lib/intel64_lin;/apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/tbb/lib/intel64/gcc4.7;/apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/daal/lib/intel64_lin;/apps/intel/parallel_studio_xe_2018.4.057/compilers_and_libraries_2018/linux/tbb/lib/intel64_lin/gcc4.4;/misc/apps/intel/compilers_and_libraries_2018.5.274/linux/compiler/lib/intel64_lin;/usr/lib/gcc/x86_64-redhat-linux/4.8.5;/usr/lib64;/lib64;/usr/lib;/lib")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
