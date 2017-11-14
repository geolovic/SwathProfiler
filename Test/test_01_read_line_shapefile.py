#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test for read_line_shapefile() function

Created on Sun Nov 12 2017
Last Modified on Sun Nov 12 2017
@author: J Vicente Perez
"""
from swath import read_line_shapefile
import matplotlib.pyplot as plt
import numpy as np
import time


def test_01():
    # Test for read_line_shapefile
    # Reads a shapefile and plots its lines
    print("Test 01 for read_line_shapefile")
    print("Testing a line shapefile with a valid name field")
    inicio = time.time()
    print("Test started")
    print("...")
    
    in_shp = "../SampleData/gisdata/sample_lines.shp"
    names = "name"
    
    lines, names = read_line_shapefile(in_shp, names)
    fin = time.time()
    print("Test finished in {0:.3f} seconds".format(fin-inicio))
    print("Plotting results")
    
    fig, ax = plt.subplots()
    for idx, line in enumerate(lines):
        lbl = names[idx]
        xy_coords = np.array(line.coords)
        ax.plot(xy_coords[:,0], xy_coords[:,1], label=lbl)
    ax.set_title("name field for names")      
    ax.legend(loc="lower right")


def test_02():
    # Test for read_line_shapefile
    # Reads a shapefile and plots its lines
    print("Test 02 for read_line_shapefile")
    print("Testing a line shapefile without a name field")
    inicio = time.time()
    print("Test started")
    print("...")
    
    in_shp = "../SampleData/gisdata/sample_lines.shp"
    
    lines, names = read_line_shapefile(in_shp)
    fin = time.time()
    print("Test finished in {0:.3f} seconds".format(fin-inicio))
    print("Plotting results")
    
    fig, ax = plt.subplots()
    for idx, line in enumerate(lines):
        lbl = names[idx]
        xy_coords = np.array(line.coords)
        ax.plot(xy_coords[:,0], xy_coords[:,1], label=lbl)
    ax.set_title("Without names field")  
    ax.legend(loc="lower right")
  
    
def test_03():
    # Test for read_line_shapefile
    # Reads a shapefile and plots its lines
    print("Test 03 for read_line_shapefile")
    print("Testing a line shapefile with a integer names field")
    inicio = time.time()
    print("Test started")
    print("...")
    
    in_shp = "../SampleData/gisdata/sample_lines.shp"
    name_field = "id"
    
    lines, names = read_line_shapefile(in_shp, name_field)
    fin = time.time()
    print("Test finished in {0:.3f} seconds".format(fin-inicio))
    print("Plotting results")
    
    fig, ax = plt.subplots()
    for idx, line in enumerate(lines):
        lbl = names[idx]
        xy_coords = np.array(line.coords)
        ax.plot(xy_coords[:,0], xy_coords[:,1], label=lbl)
        
    ax.set_title("Id as names field")    
    ax.legend(loc="lower right")    

test_01()
test_02()
test_03()
