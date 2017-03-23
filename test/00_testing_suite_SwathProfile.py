# -*- coding: utf-8 -*-
"""

José Vicente Pérez
Granada University (Spain)
November 13, 2016 

Testing suite for the pRaster Python class
"""

import time
from SwathProfile import main


def test01():
    """
    Test for main function of SwathProfile.py
    """
    print "=" * 70
    print "Test 01 para SwathProfile.py"
    inicio = time.time()
    # Input parameters
    line = "data/perfiles.shp"
    dem = "data/coskie.tif"
    names = "Name"
    width = 500
    
    drawdata = False
    drawlegend = False
    nlines = 0
    step = 0

    ret = main(line, dem, width, nlines, step, names, drawlegend, drawdata)
    
    fin = time.time()
    print "Test finalizado en " + str(fin - inicio) + " segundos"
    print "=" * 70
    
    return ret


res = test01()
