#!/usr/bin/env python
from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

version='0.0.2'

requirements = [
                "ete3",
                "datetime",
                "statistics",
                "svgwrite",
                "pycairo",
                "cairosvg",
                "argparse",

                "six"] #could not use ete3 without six

classifiers = [
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        "Intended Audience :: Healthcare Industry"
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        "Programming Language :: Python :: 3.9",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization"]

setup(
    author="Jennifer L Havens",
    author_email='jhavens@ucsd.edu',
    python_requires='>=3.7',
    classifiers=classifiers,
    description="Visualize lineages overtime, with phylogentic context, based on viral genomes",
    packages=find_packages(),#['00_scripts'],#find_packages(include=['epimuller', 'epimuller.*']),
    entry_points={
        'console_scripts': [
            'epimuller=00_scripts.mutationLinages_report:main',
            'epimuller-draw=00_scripts.drawMuller:main',
            'epimuller-define=00_scripts.defineAndCountClades:main',
            'epimuller-parse=00_scripts.parseFastaNames:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme,
    long_description_content_type="text/markdown",
    include_package_data=True,
    keywords='epimuller',
    name='epimuller',
    url='https://github.com/jennifer-bio/epimuller',
    version=version,
)