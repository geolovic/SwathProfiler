#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test for SwathGraph class

Testing different number of lines

Created on Tue Nov 14 2017
Last Modified on Tue Nov 14 2017
@author: J Vicente Perez
"""

from swath import read_line_shapefile
from swath import SwathProfile
import matplotlib.pyplot as plt
import praster as p
import time


def test_01():
    # Test for class SwathProfile
    print("Test 01 for SwathGraph class")
    print("Examining different number of lines ")
    print("Test started")
    print("...")
    
    in_dem = "../SampleData/gisdata/dem40.tif"
    in_shp = "../SampleData/gisdata/sample_lines.shp"
    name_fld = "name"
    width = 2500
    
    lines, names = read_line_shapefile(in_shp, name_fld)
    line = lines[0]
    name = names[0]
    
    dem = p.open_raster(in_dem)
    ssize = line.length / 1024.
    
    n_lines = list(range(20, 141, 20))
    n_lines.append(None)
    fig = plt.figure()
    
    for idx, nl in enumerate(n_lines):
        inicio = time.time()
        ax = fig.add_subplot(4, 2, idx+1)
        sw = SwathProfile(line, dem, width=width, n_lines=nl, step_size=ssize, name=name)
        sw.draw_swath(ax, drawbg=True, q=True)
        fin = time.time()
        if not nl:
            nl = int(sw.data.shape[1] / 2)

        ax.set_title("{0} lines ({1:.1f} seconds)".format(nl, fin-inicio), fontsize=10)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_xticks(())
        ax.set_yticks(())
    
    plt.tight_layout()

    
test_01()
        
        