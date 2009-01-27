SPAGeDi -- a program for Spatial Pattern Analysis of Genetic Diversity
  Copyright (c) 2002-2009 Olivier Hardy and Xavier Vekemans

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.


NOTE

SPAGeDi has been tested on several data sets and results were checked for con-
sistency with alternative softwares whenever possible. It may nevertheless still
contain bugs (corrected bugs are listed at the end of this manual). Some of
these bugs are probably easy to detect by causing the program to crash or lead-
ing to obvious erroneous results for particular data sets and analyses. But
others, more critical, may just cause biased results that appear plausible.
Hence, it is advised to take much care checking the consistency of the infor-
mation from the results file. The authors would appreciate being informed of any
detected bug. The authors claim no responsibility if or whenever a bug causes a
misinterpretation of the results given by SPAGeDi.

INSTALLATION

SPAGeDi requires CMake (http://www.cmake.org/) in order for it to be build from
sources.  Many Unix-like operating systems can install CMake through their
package systems.  Extract SPAGeDi and issue the following commands in the ex-
tracted directory:

cmake -g "Unix Makefiles" .
gmake
gmake install

The '-g' option to cmake can be changed provide different build system options.
See the cmake manual for information.

MANUAL

See manual.pdf for detailed instructions.

CITATION

Hardy, O. J. & X. Vekemans (2002). SPAGeDi: a versatile computer program to
analyse spatial genetic structure at the individual or population levels.
Molecular Ecology Notes 2: 618-620.

PORTABLE/UNIX PORT

Portable source code (Unix port) created by Reed A. Cartwright <reed@scit.us>

