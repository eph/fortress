export FC=ifort
#source /opt/intel/compilers_and_libraries_2017/linux/bin/compilervars.sh intel64
#source ~/.bashrc
make clean
make libfortress.so
cp libfortress.so $PREFIX/lib
mkdir -p $PREFIX/include/fortress
cp *.mod $PREFIX/include/fortress
make test_driver
cp test_driver $PREFIX/bin 
./test_driver



