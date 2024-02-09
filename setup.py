from setuptools import find_packages

from numpy.distutils.core import setup, Extension

ext1 = Extension(name='powerI4',
                 sources=['pbi4/F90_modules/powerI4.f90'],
                )

ext2 = Extension(name='bispectrumI4',
                 sources=['pbi4/F90_modules/bispectrumI4.f90'],
                library_dirs=['/opt/ohpc/pub/libs/gnu12/openmpi4/fftw/3.3.10/lib'],
                libraries = ['fftw3','fftw3_threads'],
                include_dirs= ['/opt/ohpc/pub/libs/gnu12/openmpi4/fftw/3.3.10/include'],
                )

setup(name="PBI4",
      version="1.1",
      package_dir={"": "pbi4"},
      packages=find_packages(),
      ext_modules=[ext1, ext2])