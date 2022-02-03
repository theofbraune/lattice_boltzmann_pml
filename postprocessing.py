#!/usr/bin/env python
# -*- coding: utf-8 -*-
############################################################################
# EduLB - The educational lattice Boltzmann solver
# 
# Copyright (C) 2013 Andreas Hantsch (edulb@gmx-topmail.de)
# Version 0.4
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
############################################################################

import sys, math, glob, os

timefactor=1.
velocityfactor=1000.
pngheight=600

print("Running postprocessing...")
print("Reading ./system/controlDict...")
controlDict = open("./system/controlDict")
for line in controlDict:
	if line.find("GEOFILE")>-1:
		item = line.split()
		geofile= item[1]
print("\bdone.")
controlDict.close()

print ("Reading %s..." %(geofile))
controlDict = open(geofile)
for line in controlDict:
	if line.find("LX")>-1:
		item = line.split()
		lx=float(item[1])-1
	if line.find("LY")>-1:
		item = line.split()
		ly=float(item[1])-1
print("\bdone.")
controlDict.close()

os.chdir("./results/")
print("Removing old image files...")
for file in glob.glob('*.png')+glob.glob('*.eps'):
	os.remove(file)
print("\bdone.")

number_files=len(glob.glob('*_m.ssv'))
count=0.
dc=0.0
pngwidth=int(float(lx)/float(ly)*pngheight)
print("Running gnuplot...")
for file in glob.glob('*m.ssv'):
	file_number=file.replace("_m.ssv","")
	if ((float(count)/float(number_files))>=dc):
		print ("%3.0f%%" %(float(count)/float(number_files)*100))
		dc=dc+0.1
	tempfile=open("temp.plt", 'w')
	tempfile.write("""# temporary plot file
set title \"time step %.0f\"
set xrange[0:%d]
set yrange[0:%d]
set cbrange[0:1]
unset colorbox
set palette gray
set view map
set terminal png nocrop enhanced 10 size %d,%d
set out \"%s_m.png\"
splot \"%s_m.ssv\" using 1:2:6 with pm3d,\
	\"%s_m.ssv\" using 1:2:(0):(%.15f*$3):(%.15f*$4):(0) every 4:4 ti \"\" with vectors head filled lt 1  lc rgb \"blue\"
"""
%((float(file_number)*timefactor),lx,ly,pngwidth,pngheight,file_number,file_number,file_number,velocityfactor,velocityfactor))
	tempfile.close()
	os.system("gnuplot temp.plt")
	os.remove("temp.plt")
	count=count+1.

print ("100%")
print ("\bdone.")
print ("done.")
