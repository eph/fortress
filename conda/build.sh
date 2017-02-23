export FC=gfortran
make clean
make libfortress.so
cp libfortress.so $PREFIX/lib
mkdir -p $PREFIX/include/fortress
cp *.mod $PREFIX/include/fortress
#make test_driver
#cp test_driver $PREFIX/bin 
$PYTHON setup.py clean
$PYTHON setup.py install
