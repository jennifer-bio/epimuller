#!/usr/bin/env python
from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()


requirements = ["ete3",
                "datetime",
                "statistics",
                "svgwrite",
                "random",
                "pycairo",
                "cairosvg",
                "argparse" ]
version='0.0.2'

setup(
    author="Jennifer L Havens",
    author_email='jhavens@ucsd.edu',
    python_requires='>=3.7',
    classifiers=[
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
        "Topic :: Scientific/Engineering :: Visualization"
    ],
    description="Visualize lineages with Muller plot based on viral genomes ",
    entry_points={
        'console_scripts': [
            'epimuller=epimuller.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme,
    long_description_content_type="text/markdown",
    include_package_data=True,
    keywords='epimuller',
    name='epimuller',
    packages=find_packages(include=['epimuller', 'epimuller.*']),
    url='https://github.com/jennifer-bio/epimuller',
    version=version,
)
