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

#  Version: 1.0 
#  November 23, 2016

#  Last modified March 23, 2017

import ogr
import swath
from shapely.geometry import LineString
import matplotlib.pyplot as plt
import praster as p


def main(line, dem, width, nlines, step, names, drawlegend, drawdata):   
    # Read line shapefile and raster
    centerlines, names = read_line_shapefile(line, names)

    dem_raster = p.open_raster(dem)
    if nlines == 0:
        nlines = None
    if step == 0:
        step = None

    # Create the swath profiles list
    swath_profiles = []
    n = 0
    for centerline in centerlines:
        sline = swath.SwathProfile(centerline, dem_raster, width, nlines, step, name=names[n])
        if len(sline.swaths) > 0:
            swath_profiles.append(sline)
        n += 1

    # Create graphic and draw swath profiles
    if len(swath_profiles) > 0:
        fig = plt.figure()
        swath_graph = swath.SwathGraph(swath_profiles, fig, drawlegend, drawdata)
        plt.show()
        return swath_graph


def read_line_shapefile(shapefile, names_field=""):
    """
    This function reads a line shapefile and returns a tuple (out_lines, out_names)
    
    Parameters:
    ================
    shapefile :: *str*
      Path to the line shapefile with profile centerlines
    names_field :: *str (Default: ""*
      Name of the field with the profile names. If skipped, profiles will be named sequentially
    
    Returns:
    ==============
    (out_lines, out_names) :: *tuple*
        out_lines :: List with shapely.geometry.LineString objects representing centerlines
        out_names :: List of string with profile names
    """
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataset = driver.Open(shapefile)
    layer = dataset.GetLayer()
    
    # layerdef = layer.GetLayerDefn()
    if not layer.GetGeomType() == 2:
        return

    out_names = []
    out_lines = []
    perfil_id = 0
    for feat in layer:
        out_names.append(feat.GetField(names_field))
        geom = feat.GetGeometryRef()
        coords = []
        for n in range(geom.GetPointCount()):
            pt = geom.GetPoint(n)
            coords.append((pt[0], pt[1]))
        out_lines.append(LineString(coords))
        perfil_id += 1
        
    return out_lines, out_names
