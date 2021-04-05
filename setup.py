import setuptools

with open("README.rst", "r") as fh:
    long_description = fh.read()

version = {}
with open("lcheapo_obspy/version.py") as fp:
    exec(fp.read(), version)

setuptools.setup(
    name="lcheapo_obspy",
    version=version['__version__'],
    author="Wayne Crawford",
    author_email="crawford@ipgp.fr",
    description="LCHEAPO data reading and plotting",
    long_description=long_description,
    long_description_content_type="text/x-rst; charset=UTF-8",
    url="https://github.com/WayneCrawford/lcheapo_obspy",
    # packages=setuptools.find_packages(),
    packages=['lcheapo_obspy'],
    package_dir={'lcheapo_obspy': 'lcheapo_obspy'},
    package_data={'lcheapo_obspy': ['data/*.xml', 'data/*.json',
                                    '_examples/*.py', '_examples/*.yaml']},
    # include_package_data=True,
    install_requires=[
          'obspy>=1.1',
          'pyyaml>5.0'
          'jsonschema>=2.6',
          'jsonref>=0.2',
          'progress>=1.5',
          'lcheapo>=0.73'
      ],
    entry_points={
         'console_scripts': [
             'lcplot=lcheapo_obspy.lcread:_plot_command',
             'lc2SDS_weak=lcheapo_obspy.lc2SDS:lc2SDS',
             'lctest=lcheapo_obspy.lctest:main',
             'lc_examples=lcheapo_obspy.lcputexamples:main'
         ]
    },
    python_requires='>=3.6',
    classifiers=(
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics"
    ),
    keywords='oceanography, marine, OBS'
)
