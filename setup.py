from setuptools import setup, find_packages

setup(
    name='archdbmap',
    version='v0.1',

    description='Map Smotifs from ArchDB14 into sequences by sequence similarity',
    long_description='',

    # The project's main homepage.
    url='https://github.com/structuralbioinformatics/archdbmap',

    # Author details
    author='Jaume Bonet',
    author_email='jaume.bonet@gmail.com',

    # Choose your license
    license='MIT',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Environment :: Console',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Operating System :: MacOS',
        'Operating System :: Unix',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],

    project_urls={
        'Source': 'https://github.com/structuralbioinformatics/archdbmap',
        'Tracker': 'https://github.com/structuralbioinformatics/archdbmap/issues',
    },

    platforms='UNIX',
    keywords='development',

    install_requires=[x.strip() for x in open('REQUIREMENTS').readlines()],

    packages=find_packages(exclude=['docs', 'demo', 'sphinx-docs']),
    include_package_data=True,
    package_data={
        'archdbmap': ['REQUIREMENTS'],
    },
    entry_points={
        'console_scripts':
            ['archdbmap=archdbmap.__main__:main',
             ]
    },

    zip_safe=False,
)
