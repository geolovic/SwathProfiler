#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test for SwathGraph class

Created on Tue Nov 14 2017
Last Modified on Tue Nov 14 2017
@author: J Vicente Perez
"""

from swath import read_line_shapefile
from swath import SwathProfile, SwathGraph
import matplotlib.pyplot as plt
import praster as p


def test_01():
    # Test for class SwathProfile
    print("Test 01 for SwathGraph class")
    print("Creating swath Graph")
    print("Test started")
    print("...")
    
    in_dem = "../SampleData/gisdata/dem40.tif"
    in_shp = "../SampleData/gisdata/sample_lines.shp"
    name_fld = "name"
    width = 2500
    
    lines, names = read_line_shapefile(in_shp, name_fld)
    dem = p.open_raster(in_dem)
    swath_profiles = []
    n_lines = 50
    
    for idx, line in enumerate(lines):
        ssize = line.length / 1024.
        print("Processing {0} of {1} lines".format(idx+1, len(lines)))
        swath_profiles.append(SwathProfile(line, dem, width=width, n_lines=n_lines, 
                                           step_size=ssize, name=names[idx]))
   
    fig = plt.figure()
    return SwathGraph(swath_profiles, fig)


my_graph = test_01()
