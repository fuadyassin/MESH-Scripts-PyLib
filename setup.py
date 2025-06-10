from setuptools import setup, find_packages

setup(
    name='MESH-Scripts-PyLib',
    version='0.1',
    package_dir={'': 'src'},  # Tells setuptools to look inside /src
    packages=find_packages(where='src'),
    install_requires=[
        'pandas',
        'numpy',
        'geopandas',
        'netCDF4',
        'requests',
        'xarray',
        'shapely',
        'owslib'
    ],
    description='Python package for processing and analyzing hydrological data.',
    author='Fuad Yassin',
    author_email='fuad.yassin@usask.ca',
    url='https://github.com/MESH-Model/MESH-Scripts-PyLib',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)
