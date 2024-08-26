module load PE-gnu/3.0
module load mkl/2018.1.163
module load cmake/3.20.3

# PFLOTRAN build
git clone https://github.com/bsulman/pflotran-elm-interface.git
pushd pflotran-elm-interface/src/pflotran
git checkout pflotran-elm-interface
export PETSC_DIR=/software/user_tools/current/cades-ccsi/packages_PE-gnu3.0_ext/petsc-x-noopt/openmpi-3.1.5-gcc-8.1.0
make pflotran pflotran_rxn
cp libpflotranchem.a /nfs/data/ccsi/proj-shared/b0u/ELM-PFLOTRAN/alquimia_PE3/

# Alquimia build
export PETSC_DIR=/software/user_tools/current/cades-ccsi/packages_PE-gnu3.0_ext/petsc-x-noopt/openmpi-3.1.5-gcc-8.1.0
export PFLOTRAN_DIR=/home/b0u/models/PFLOTRAN/pflotran_PE3/pflotran-elm-interface/src/pflotran
PETSC_ARCH='' cmake .. \
-DCMAKE_INSTALL_PREFIX=/nfs/data/ccsi/proj-shared/b0u/ELM-PFLOTRAN/alquimia_PE3 \
-DCMAKE_C_COMPILER=$OPENMPI_DIR/bin/mpicc \
-DCMAKE_CXX_COMPILER=$OPENMPI_DIR/bin/mpicxx \
-DCMAKE_Fortran_COMPILER=$OPENMPI_DIR/bin/mpif90 \
-DCMAKE_BUILD_TYPE=Debug \
-DXSDK_WITH_PFLOTRAN=ON \
-DTPL_PFLOTRAN_LIBRARIES=/nfs/data/ccsi/proj-shared/b0u/ELM-PFLOTRAN/alquimia_PE3/libpflotranchem.a \
-DTPL_PFLOTRAN_INCLUDE_DIRS=$PFLOTRAN_DIR

# E3SM is set up assuming these are all in the same directory. But we could change that
ln -s /nfs/data/ccsi/proj-shared/b0u/ELM-PFLOTRAN/alquimia_PE3/lib/*.so /nfs/data/ccsi/proj-shared/b0u/ELM-PFLOTRAN/alquimia_PE3/
# We should pull request to alquimia to include the .mod files in installation include
cp alquimia/*.mod /nfs/data/ccsi/proj-shared/b0u/ELM-PFLOTRAN/alquimia_PE3/