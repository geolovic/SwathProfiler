# -*- coding: utf-8 -*-
"""

José Vicente Pérez
Granada University (Spain)
November 13, 2016 

Testing suite for the pRaster Python class
"""

import os, sys
import time

#change to test alternative praster modules
from SwathProfile import main


#Cambiamos el working directoy
working_folder = os.path.dirname(sys.argv[0])
os.chdir(working_folder)

def test01():
    """
    Test for main function of SwathProfile.py
    """
    print "=" * 70
    print "Test 01 para main() function of SwathProfile.py"
    inicio = time.time()
    #Input parameters
    line = "data/in/perfiles.shp"
    dem = "data/in/coskie.tif"
    names = "Name"
    width = 500
    
    drawdata = True
    drawlegend = False
    nlines = 0
    step = 0

    ret = main(line, dem, width, nlines, step, names, drawlegend, drawdata)
    
    fin = time.time()
    print "Test finalizado en " + str(fin - inicio) + " segundos"
    print "=" * 70
    
    return ret
    


res = test01()