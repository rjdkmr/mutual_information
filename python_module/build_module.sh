#!/bin/sh

# Change this environment variable according to the system
export PYTHON_INCLUDE=-I/usr/include/python2.7
export PYTHON_LIB=-I/usr/lib64/python2.7

rm MI_depend.o digamma.o mi_swig_wrap.o mi_swig_wrap.cxx mi.py Mutual_information.o _mi.so

echo "========================"
echo "Making python module..."
echo "========================"

swig -python -c++ mi_swig.i

g++ -fPIC -pthread -c ../src/Mutual_information.cpp
g++ -fPIC -pthread -c ../src/MI_depend.cpp
g++ -fPIC -pthread -c ../src/digamma.cpp
g++ -fPIC -pthread -c mi_swig_wrap.cxx ${PYTHON_INCLUDE} ${PYTHON_LIB}

g++ -shared -lpthread Mutual_information.o MI_depend.o digamma.o mi_swig_wrap.o -o _mi.so

rm MI_depend.o digamma.o mi_swig_wrap.o mi_swig_wrap.cxx Mutual_information.o


echo "Finished...."
echo "   "
echo "To use python module in the python script, make sure mi.py and _mi.so should be in PYTHONPATH environment variable"
echo "   "
echo "   "

echo "========================"
echo "TEST"
echo "========================"
python ../test/example_mi.py ../test/proj_1.dat ../test/proj_2.dat
echo " "
echo "========================"
echo "If there was ERROR in TEST, TEST was not successful"
echo "========================"
