import sys 
if sys.version_info[0] < 3:
    sys.exit('Sorry, Python < 3.x is not supported')

from setuptools import setup, find_packages, Extension
import numpy
from Cython.Distutils import build_ext

setup(name="fast_dbscan",
      packages=find_packages(),
      cmdclass={'build_ext': build_ext},
      version='0.1.1',
      description="fast dbscan clustering on peptide strings",
      #long_description=open("README.md").read(),
      author='Michael J. Harms',
      author_email='harmsm@gmail.com',
      url='https://github.com/harmslab/fast_dbscan',
      download_url='https://github.com/harmslab/fast_dbscan/archive/0.1.0.tar.gz',
      ext_modules=[Extension("dbscan",
                   sources=["src/dbscan_module.pyx", "src/dbscan.c"],
                   include_dirs=["src/",numpy.get_include()])],
      install_requires=["numpy"],
      zip_safe=False,
      package_data={"":["*.h","src/*.h"]},
      classifiers=['Programming Language :: Python'],
      entry_points = {
            'console_scripts': [
                'fast_dbscan = fast_dbscan.fast_dbscan:main',
            ]
      })
