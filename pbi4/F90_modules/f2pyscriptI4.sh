#!/bin/sh  

rm *.so

f2py -c powerI4.f90 -m powerI4
f2py -c -L/usr/local/Cellar/fftw/3.3.10_1/lib  -lfftw3  -lfftw3_threads -I/usr/local/Cellar/fftw/3.3.10_1/include bispectrumI4.f90 -m bispectrumI4

mv *.so ../

