from setuptools import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize("tests/cython_tests/cpython_test.pyx")
)