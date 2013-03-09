SPAGeDi
=======

SPAGeDi: Spatial Pattern Analysis of Genetic Diversity  
Copyright &copy; 2002&ndash;2009 Olivier Hardy (<ohardy@ulb.ac.be>) and Xavier Vekemans

Portable/Unix source code maintained by Reed A. Cartwright (<cartwright@asu.edu>).

Website: <http://ebe.ulb.ac.be/ebe/SPAGeDi.html>  
GitHub: <https://github.com/reedacartwright/spagedi>  
Downloads: <http://scit.us/spagedi/>  

License
-------

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but **without any warranty**; without even the implied warranty of **merchantability** or **fitness for a particular purpose**.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

Citation
--------

[Hardy OJ & Vekemans X (2002)](http://iee.ulb.ac.be/sciences/lagev/fichiers/Spagedi_MENotes2002.pdf) SPAGeDi: a versatile computer program to analyse spatial genetic structure at the individual or population levels. Molecular Ecology Notes 2: 618&ndash;620.

Nota Bene
---------

SPAGeDi has been tested on several data sets and results were checked for consistency with alternative softwares whenever possible. It may nevertheless still contain bugs (corrected bugs are listed at the end of the manual). Some of these bugs are probably easy to detect by causing the program to crash or leading to obvious erroneous results for particular data sets and analyses. But others, more critical, may just cause biased results that appear plausible. Hence, it is advised to take much care checking the consistency of the information from the results file. The authors would appreciate being informed of any detected bug. The authors claim no responsibility if or whenever a bug causes a misinterpretation of the results given by SPAGeDi.

Using SPAGeDi
-------------

See manual.pdf for detailed instructions on using SPAGeDi.

Downloading and Installing
--------------------------

Binary packages for SPAGeDi can be downloaded from <http://scit.us/spagedi/>.

Compiling from Source Code
--------------------------

SPAGeDi requires CMake 2.8 (<http://www.cmake.org/>) to build it from sources.  Many Unix-like operating systems can install CMake through their package systems.

Download the SPAGeDi source code (<http://scit.us/spagedi/>) and issue the following commands:

	tar xvzf SPAGeDi-.*.tar.bz2
	cd SPAGeDi-.*/build
	cmake ..
	make
	make install

You may have to adjust the first two commands to fit the archive that you download.

The '-G' option to cmake is used to specify different build systems, e.g. Unix Makefiles versus KDevelop3 project.  The '-D' option to cmake can be used to set different cmake variables from the command line:

	cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr .
	make
	make install

This will build an optimized version of SPAGeDi and install it to '/usr/bin'. To specify your own build flags you need to set the environment variables CFLAGS and LDFLAGS as neccessary.  Then specify

	cmake -DCMAKE_BUILD_TYPE= .

See CMake's manual for additional information.

