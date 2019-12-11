import setuptools

with open("README.rst", "r") as fh:
    long_description = fh.read()
    
version={}
with open("lcheapo_obspy/version.py") as fp:
    exec(fp.read(),version)

setuptools.setup(
    name="lcheapo_obspy",
    version=version['__version__'],
    author="Wayne Crawford",
    author_email="crawford@ipgp.fr",
    description="LCHEAPO data reading and plotting",
    long_description=long_description,
    long_description_content_type="text/x-rst; charset=UTF-8",
    url="https://github.com/WayneCrawford/lcheapo_obspy",
    packages=setuptools.find_packages(),
    include_package_data=True,
    install_requires=[
          'obspy>=1.1',
          'lcheapo>=0.7',
          'numpy>=1.17',
          'scipy>=1.3',
          'matplotlib>3.0',
          'PyYAML>3.0'
      ],
    entry_points={
         'console_scripts': [
             'lcplot=lcheapo_obspy.lc_read:_plot_command',
             'lc2ms_w=lcheapo_obspy.lc_read:_to_mseed_command',
             'lctest=lcheapo_obspy.lctest:main'
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