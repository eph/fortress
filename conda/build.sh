#export PREFIX=/home/eherbst/anaconda3/
# make clean
which $PYTHON
$PYTHON --version
$PYTHON setup.py build
$PYTHON setup.py install
export CONDA_BUILD=1
export FC=gfortran
make clean
make libfortress.so
cp libfortress.so $PREFIX/lib
mkdir -p $PREFIX/include/fortress
cp *.mod $PREFIX/include/fortress
#make test_driver
#cp test_driver $PREFIX/bin 
