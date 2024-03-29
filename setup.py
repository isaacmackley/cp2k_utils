# setup.py

from setuptools import setup, find_packages

setup(
    name='cp2k_utils',
    version='0.1',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'xyz-to-cif = cp2k_utils.cli.structures:main',
            'pdos-plot = cp2k_utils.cli.pdos:main'
        ],
    },
    install_requires=[
        'numpy',
        'pymatgen'
    ],
)
