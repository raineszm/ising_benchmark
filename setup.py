from setuptools import setup
setup(name='ising',
      install_requires=['cffi>=1.0'],
      cffi_modules=['ising/random_build.py:ffi'])
