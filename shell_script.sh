if [ ! -d ~/tmp ]; then
  mkdir ~/tmp
  chmod 777 ~/tmp
fi
export TMP=~/tmp
export TMPDIR=~/tmp
source /etc/profile.d/modules.sh
source /etc/profile.d/pbs.sh
cp ~/GTD2D_lib/F/gtd2d_libinter_cfunvec.a .
chmod u=rwx num_core.sh
module load intel-suite mpi fftw/3.3.3-double netcdf/4.4.4-fortran netcdf/4.4.1-c matlab/R2019a
make all install

qsub -q pqaero CX1.pbs
sleep 10
