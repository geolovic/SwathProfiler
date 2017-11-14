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

## DEBUG PARAMETERS
## ================
line_shp = "../SampleData/gisdata/sample_lines.shp"
dem = "../SampleData/gisdata/dem40.tif"
name_field = "name"
width = 2500
nlines = 0
step = 0
fullres = False

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
