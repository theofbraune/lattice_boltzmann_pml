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

import os,sys
from PIL import Image
import pdb

try:
	if len(sys.argv) != 2:
		sys.stderr.write("Usage: %s <inputfile>%s"%(sys.argv[0],os.linesep))
	else:
		print("Converting PNG to monochrome BMP file.")
		pngfile=sys.argv[1]
		
		bmpfile=pngfile.replace("png","bmp")
		geofile=pngfile.replace(".png","")
		print('png: ', pngfile)
		print('bmp: ',bmpfile)
		#pdb.set_trace()
		os.system("convert -monochrome %s %s" %(pngfile, bmpfile))
		print("Reading BMP file.")
		im=Image.open(bmpfile)
		pix=im.load()
	
		out=[]
		
		(width,height)=im.size
		for x in range(0,width):
			for y in range(0,height):
				dot=pix[x,y]
				#if(dot>0):
				#	print('dot value: ', dot)
				dot//=255
				if (dot == 0) or (dot == 1):
					"""
					flip vertically: obstacle coordinates start
					at bottom, image coordinates at top
					"""
					if(dot==0):
						dot = 1
						#pass
					else:
						dot=0
						#pass
					out.append("%s %s %s"%(x,height-y-1,dot))
				else:
					print('dot: ', dot)
					print("wrong data")
					raise Exception
		print("Writing numerical data to geometry file.")
		outfile=open(geofile, 'w')
		outfile.write("LX %d\n" %width)
		outfile.write("LY %d\n" %height)
		for dot in out:
			outfile.write(dot)
			outfile.write("\n")
		outfile.close()
		
		print("Deleting temporary BMP file.")
		os.system("rm %s" %bmpfile)
		
		print("Done.")
except:
	pass
"""
except Exception,e:
	print "Error: %s"%e
"""