#    Copyright (C) 2012 Daniel Gamermann <gamermann@gmail.com>
#
#    This file is part of ExGUtils
#
#    ExGUtils is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    ExGUtils is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with ExGUtils.  If not, see <http://www.gnu.org/licenses/>.
#
#    
#    Please, cite us in your reasearch!
#


from distutils.core import setup, Extension

ld = """This package contains modules with functions for statistical work
        with the ex-gaussian function and some routines for numerical 
        calculations."""

bla = raw_input("Compile C [Y/n]? ")

if bla.strip() == "Y":
    setup(name = "ExGUtils",
          version = "3.0",
          description = "Tools for statistical analyses using the ex-gaussian function (C sources compiled).",
          long_description = ld,
          url = "http://pypi.python.org/pypi/ExGUtils",
          author = "Daniel Gamermann",
          author_email = "gamermann@gmail.com",
          license = "GNU GPLv3",
          package_dir={'ExGUtils' : 'ExGUtils'},
          packages=['ExGUtils'],
          ext_modules=[Extension("ExGUtils.uts", sources=["ExGUtils/csources/_Uts.c", "ExGUtils/csources/routs.c", "ExGUtils/csources/matrices.c", "ExGUtils/csources/randoms.c", "ExGUtils/csources/integral.c"])],
          requires = ["math", "scipy", "numpy"],
          classifiers=[
              'Programming Language :: Python',
              'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
              'Operating System :: OS Independent',
              'Intended Audience :: Science/Research',
              'Intended Audience :: End Users/Desktop',
              'Intended Audience :: Developers',
              'Topic :: Scientific/Engineering :: Mathematics',
              'Topic :: Scientific/Engineering :: Bio-Informatics'
              ])
else:
    setup(name = "ExGUtils",
          version = "3.0",
          description = "Tools for statistical analyses using the ex-gaussian function (pure python sources).",
          long_description = ld,
          url = "http://pypi.python.org/pypi/ExGUtils",
          author = "Daniel Gamermann",
          author_email = "gamermann@gmail.com",
          license = "GNU GPLv3",
          package_dir={'ExGUtils' : 'PYExGUtils'},
          packages=['ExGUtils'],
          requires = ["random", "math", "scipy", "numpy"],
          classifiers=[
              'Programming Language :: Python',
              'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
              'Operating System :: OS Independent',
              'Intended Audience :: Science/Research',
              'Intended Audience :: End Users/Desktop',
              'Intended Audience :: Developers',
              'Topic :: Scientific/Engineering :: Mathematics',
              'Topic :: Scientific/Engineering :: Bio-Informatics'
              ])


