set(CMAKE_Fortran_COMPILER "/opt/cray/pe/craype/2.6.3/bin/ftn")
set(CMAKE_Fortran_COMPILER_ARG1 "")
set(CMAKE_Fortran_COMPILER_ID "Intel")
set(CMAKE_Fortran_COMPILER_VERSION "19.0.5.20190815")
set(CMAKE_Fortran_COMPILER_WRAPPER "CrayPrgEnv")
set(CMAKE_Fortran_PLATFORM_ID "Linux")
set(CMAKE_Fortran_SIMULATE_ID "")
set(CMAKE_Fortran_SIMULATE_VERSION "")




set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_Fortran_COMPILER_AR "")
set(CMAKE_RANLIB "/usr/bin/ranlib")
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





set(CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES "/opt/cray/pe/libsci/19.06.1/INTEL/16.0/x86_64/include;/opt/cray/pe/mpt/7.7.11/gni/mpich-intel/16.0/include;/opt/cray/rca/2.2.20-7.0.2.1_2.93__g8e3fb5b.ari/include;/opt/cray/alps/6.6.59-7.0.2.1_3.85__g872a8d62.ari/include;/opt/cray/xpmem/2.2.20-7.0.2.1_2.72__g87eb960.ari/include;/opt/cray/gni-headers/5.0.12.0-7.0.2.1_2.34__g3b1768f.ari/include;/opt/cray/pe/pmi/5.0.15/include;/opt/cray/ugni/6.0.14.0-7.0.2.1_3.77__ge78e5b0.ari/include;/opt/cray/udreg/2.3.2-7.0.2.1_2.52__g8175d3d.ari/include;/opt/cray/wlm_detect/1.3.3-7.0.2.1_2.23__g7109084.ari/include;/opt/cray/krca/2.2.7-7.0.2.1_2.81__ge897ee1.ari/include;/opt/cray-hss-devel/9.0.0/include;/opt/intel/compilers_and_libraries_2019.5.281/linux/ipp/include;/opt/intel/compilers_and_libraries_2019.5.281/linux/mkl/include;/opt/intel/compilers_and_libraries_2019.5.281/linux/pstl/include;/opt/intel/compilers_and_libraries_2019.5.281/linux/tbb/include;/opt/intel/compilers_and_libraries_2019.5.281/linux/daal/include;/opt/intel/compilers_and_libraries_2019.5.281/linux/compiler/include/intel64;/opt/intel/compilers_and_libraries_2019.5.281/linux/compiler/include/icc;/opt/intel/compilers_and_libraries_2019.5.281/linux/compiler/include;/usr/local/include;/usr/lib64/gcc/x86_64-suse-linux/7/include;/usr/lib64/gcc/x86_64-suse-linux/7/include-fixed;/usr/include")
set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "AtpSigHandler;AtpSigHCommData;pthread;darshan;fmpich;mpichcxx;darshan;z;mpichf90_intel;rt;ugni;pthread;pmi;imf;m;dl;sci_intel_mpi;sci_intel;imf;m;dl;mpich_intel;rt;ugni;pthread;pmi;imf;m;dl;pmi;pthread;alpslli;pthread;wlm_detect;alpsutil;pthread;rca;xpmem;ugni;pthread;udreg;sci_intel;imf;m;pthread;dl;hugetlbfs;imf;m;pthread;ifport;ifcore;imf;svml;m;ipgo;irc;svml;c;gcc;irc_s;dl;c")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/opt/cray/pe/libsci/19.06.1/INTEL/16.0/x86_64/lib;/opt/cray/dmapp/default/lib64;/opt/cray/pe/mpt/7.7.11/gni/mpich-intel/16.0/lib;/opt/cray/rca/2.2.20-7.0.2.1_2.93__g8e3fb5b.ari/lib64;/opt/cray/alps/6.6.59-7.0.2.1_3.85__g872a8d62.ari/lib64;/opt/cray/xpmem/2.2.20-7.0.2.1_2.72__g87eb960.ari/lib64;/opt/cray/pe/pmi/5.0.15/lib64;/opt/cray/ugni/6.0.14.0-7.0.2.1_3.77__ge78e5b0.ari/lib64;/opt/cray/udreg/2.3.2-7.0.2.1_2.52__g8175d3d.ari/lib64;/sw/gaea-cle7/darshan/3.2.1-1/runtime/lib;/opt/cray/pe/atp/2.1.3/libApp;/opt/cray/wlm_detect/1.3.3-7.0.2.1_2.23__g7109084.ari/lib64;/opt/intel/compilers_and_libraries_2019.5.281/linux/mpi/intel64/libfabric/lib;/opt/intel/compilers_and_libraries_2019.5.281/linux/ipp/lib/intel64;/opt/intel/compilers_and_libraries_2019.5.281/linux/compiler/lib/intel64_lin;/opt/intel/compilers_and_libraries_2019.5.281/linux/mkl/lib/intel64_lin;/opt/intel/compilers_and_libraries_2019.5.281/linux/tbb/lib/intel64/gcc4.7;/opt/intel/compilers_and_libraries_2019.5.281/linux/daal/lib/intel64_lin;/sw/gaea-cle7/uasw/ncrc/envs/20200417/opt/linux-sles15-x86_64/gcc-7.5.0/globus-toolkit-6.0.17-klqyvmmhxqsf77ita7vvlw3wgyire7df/lib;/usr/lib64/gcc/x86_64-suse-linux/7;/usr/lib64;/lib64;/usr/x86_64-suse-linux/lib;/lib;/usr/lib")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
