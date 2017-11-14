# -*- coding: utf-8 -*-
#
#  SwathProfile.py
#  
#  Copyright (C) 2016  J. Vicente Perez, Universidad de Granada
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.

#  For additional information, contact to:
#  Jose Vicente Perez Pena
#  Dpto. Geodinamica-Universidad de Granada
#  18071 Granada, Spain
#  vperez@ugr.es // geolovic@gmail.com

#  Version: 1.1
#  November 14, 2017

#  Last modified November 14, 2017

import swath as sw
import matplotlib.pyplot as plt
import praster as p
import argparse


# ARGUMENT PARSER
# ===============
parser = argparse.ArgumentParser()
parser.add_argument("line_shapefile", help="Shapefile with profile center lines")
parser.add_argument("dem", help="Digital Elevation Model")
parser.add_argument("width", help="Half width of the swath", type=float)
parser.add_argument("-n", "--name_field", help="Field with names", default="")
parser.add_argument("-l", "--nlines", help="Number of profile lines at each side", type=int, default=0)
parser.add_argument("-s", "--stepsize", help="Step size to take elevations along profiles", type=float, default=0.) 
parser.add_argument("-f", "--fullres", help="Full resolution (much slower!!)", action="store_true")
args = parser.parse_args()


# PARAMETERS
# ==========
line_shp = args.line_shapefile
dem = args.dem
name_field = args.name_field
width = args.width
nlines = args.nlines
step = args.stepsize
fullres = args.fullres


def main(line_shp, dem, name_field, width, nlines, step, fullres):   
    
    # Read line shapefile and get names
    centerlines, names = sw.read_line_shapefile(line_shp, name_field)

    # Open dem
    dem_raster = p.open_raster(dem)
    if not dem_raster:
        print("It was a problem opening the DEM")
        return
    
    # If nlines and step are 0, get default values
    if nlines == 0:
        nlines = None
    if step == 0:
        step = None
    
    # If fullresolution, take maximum resolution and discard previous values
    if fullres:
        nlines = int(width/dem.cellsize)
        step = dem.cellsize

    # Create the swath profiles list
    swath_profiles = []
    
    for idx, line in enumerate(centerlines):
        sline = sw.SwathProfile(line, dem_raster, width=width, n_lines=nlines, 
                                           step_size=step, name=names[idx])      
        if len(sline.data) > 0:
            swath_profiles.append(sline)

    # Create graphic and draw swath profiles
    if len(swath_profiles) > 0:
        fig = plt.figure()
        return sw.SwathGraph(swath_profiles, fig)

if __name__ == "__main__":
    # execute only if run as a script
    my_graph = main(line_shp, dem, name_field, width, nlines, step, fullres)   
    plt.show()
