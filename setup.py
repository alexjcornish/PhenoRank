# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('LICENSE') as f:
    license = f.read()

setup(
    name='phenorank',
    version='0.1.5',
    description='Reducing study bias in gene prioritisation through simulation',
    author='Alex J. Cornish',
    license=license,
    packages=['phenorank'],
    package_dir={'phenorank': 'phenorank'},
    package_data={'phenorank': ['data_phenorank/*', 'data_prince/*', 'test/data_phenorank/*', 'test/data_prince/*']}
)
