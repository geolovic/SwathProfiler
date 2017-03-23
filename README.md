# SwathProfiler
This tool calculates the swath profiles for a set of lines for a given swath witdh. After Swath Profiles are calculated, it will show a graphic window with the output swaths. In that window use left and right mouse button to navigate between the different profiles. 
It also calculates THI* index (see Pérez-Peña et al., 2016).

## Instalation instructions:
1. Open QGIS and go to: Scripts >> Tools >> Add script from file (To activate the Processing Toolbox go to Menu >> Processing >> Toolbox or press Ctrl + Alt + T)
2. Navigate to the folder qgis_release and select SwathProfile.py
3. The new script will be now accessible in Scripts >> User scripts >> SwathProfiler
4. To add the help file, copy SwathProfile.py.help to UserFolder/.qgis2/processing/script
5. Double clic to launch the script. In the display window use left and right buttons to move to next / previous profile (if input line shapefile contained more than a linea). 


References:
Pérez-Peña, J.V., Al-Awabdeh, M., Azañón, J.M., Galve, J.P., Booth-Rea, G., Notti, D., 2016. SwathProfiler and NProfiler: Two new ArcGIS Add-ins for the automatic extraction of swath and normalized river profiles. Comput. Geosci. doi:10.1016/j.cageo.2016.08.008
