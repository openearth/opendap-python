from setuptools import setup, find_packages

setup(
    name='OpenDAP',
    version='0.0',
    author='Bas Hoonhout',
    author_email='bas.hoonhout@deltares.nl',
    packages=find_packages(),
    description='OpenDAP toolbox',
    long_description=open('README.txt').read(),
    install_requires=['netCDF4',
                      'pandas',
                      'datetime',
                      'logging'],
)
