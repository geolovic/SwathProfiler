#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test for SwathProfile class

Created on Sun Nov 12 2017
Last Modified on Tue Nov 14 2017
@author: J Vicente Perez
"""
from swath import read_line_shapefile
from swath import SwathProfile
import matplotlib.pyplot as plt
import time
import praster as p


def test_01():
    # Test for class SwathProfile
    print("Test 01 for SwathProfile class")
    print("Creating a simple swath profile")
    print("Using 50 lines left and right and 1024 elev values.")
    inicio = time.time()
    print("Test started")
    print("...")
    
    in_dem = "../SampleData/gisdata/dem40.tif"
    in_shp = "../SampleData/gisdata/sample_lines.shp"
    names = "name"
    width = 2500
    
    lines, names = read_line_shapefile(in_shp, names)
    line = lines[0]
    name = names[0]
    dem = p.open_raster(in_dem)
    
    ssize = line.length / 1024.
    
    sw = SwathProfile(line, dem, width=width, n_lines=50, step_size=ssize, name=name)
    fin = time.time()
    print("Profile created in {0:.3f} seconds".format(fin-inicio))
    
    print("Plotting results")
    
    fig, ax = plt.subplots()
    sw.draw_swath(ax, drawdata=False, q=True, r=True, legend=True)


def test_02():
    # Test for class SwathProfile
    print("="*20)
    print("Test 02 for SwathProfile class")
    print("Creating a simple swath profile")
    print("Using full resolution")
    inicio = time.time()
    print("Test started")
    print("...")
    
    in_dem = "../SampleData/gisdata/dem40.tif"
    in_shp = "../SampleData/gisdata/sample_lines.shp"
    names = "name"
    width = 2500
    
    lines, names = read_line_shapefile(in_shp, names)
    line = lines[0]
    name = names[0]
    dem = p.open_raster(in_dem)
    
    sw = SwathProfile(line, dem, width=width, n_lines=None, step_size=None, name=name)
    fin = time.time()
    print("Profile created in {0:.3f} seconds".format(fin-inicio))
    
    print("Plotting results")
    
    fig, ax = plt.subplots()
    sw.draw_swath(ax, drawdata=False, q=True, r=True, legend=True)
    
    
test_01()
#test_02()
