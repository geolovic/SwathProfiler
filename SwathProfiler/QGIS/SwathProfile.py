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


import matplotlib.pyplot as plt
import gdal
import numpy as np
import ogr
from shapely.geometry import LineString
import matplotlib.patches as mpatches


# QGIS ARGUMENTS
# ==============
##Line_shapefile=vector
##DEM=raster
##Width=number 1000
##Name_field=field Line_shapefile
##Number_of_lines=number 0
##Step_size=number 0
##Full_resolution=boolean False

line_shp = str(Line_shapefile)
dem = str(DEM)
name_field = str(Name_field)
width = float(Width)
nlines = int(Number_of_lines)
step = float(Step_size)
fullres = bool(Full_resolution)


# IMPORTED MODULES [November 14, 2017]
# ====================================
# swath.py
# ========
# Disable some keymap characters that interfere with graph key events
plt.rcParams["keymap.xscale"] = [""]
plt.rcParams["keymap.yscale"] = [""]
plt.rcParams["keymap.save"] = [u'ctrl+s']
DIRECTIONS = {"LEFT": -1, "RIGHT": 1}

def read_line_shapefile(shapefile, names_field=""):
    """
    This function reads a line shapefile and returns a tuple (out_lines, out_names)
    
    Parameters:
    ================
    shapefile : *str*
      Path to the line shapefile with profile centerline
    names_field : *str*
      Name of the field with the profile names. If skipped, profiles will be named sequentially
    
    Returns:
    ==============
    (out_lines, out_names) : *tuple*
        out_lines : List with shapely.geometry.LineString objects representing shapefile lines
        out_names : List of string with profile names
    """
    # Open the dataset and get the layer
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataset = driver.Open(shapefile)
    layer = dataset.GetLayer()
    
    # Check if layer has the right geometry
    if not layer.GetGeomType() == 2:
        return

    # Get a list of layer fields
    layerdef = layer.GetLayerDefn()
    fields = [layerdef.GetFieldDefn(idx).GetName() for idx in range(layerdef.GetFieldCount())]
    take_field = True
    if not names_field in fields:
        take_field = False

    out_names = []
    out_lines = []
    perfil_id = 0
    
    for feat in layer:
        if take_field:
            out_names.append(str(feat.GetField(names_field)))
        else:
            out_names.append(str(perfil_id))
        
        geom = feat.GetGeometryRef()
        # If the feature is multipart, only the first part is considered
        if geom.GetGeometryCount() > 0:
            geom = geom.GetGeometryRef(0)
        
        coords = []
        for n in range(geom.GetPointCount()):
            pt = geom.GetPoint(n)
            coords.append((pt[0], pt[1]))
        out_lines.append(LineString(coords))
        perfil_id += 1
        
    return out_lines, out_names


class SwathProfile:
    def __init__(self, line=None, dem=None, width=0, n_lines=None, step_size=None, name=""):
        """
        Class to create a swath profile object and related parameters

        :param line: shapely.geometry.LineString - LineString the swath profile centerline
        :param dem: praster.pRaster - pRaster with the Digital Elevation Model
        :param width: float - Half-width of the swath profile (in data units)
        :param n_lines: int - number of lines at each side of the swath centerline
        :param step_size: float - Step-size to get elevation points along the profile
        :param name: str - Name of the profile
        """

        self.name = str(name)
        
        # Creates an empty SwathProfile Object
        if line is None:
            return
        
        # Get step size (By default dem.cellsize if was not specified)
        if step_size is None or step_size < dem.cellsize:
            step_size = dem.cellsize
        
        # Get number of lines (By default width/dem.cellsize)
        if n_lines is None or n_lines > int(width/dem.cellsize):
            n_lines = int(width/dem.cellsize)
        
        # Get distance between lines
        line_distance = float(width) / n_lines
        
        # Get profile distances
        self.li = np.arange(0., line.length, step_size)
        
        # Get the number of points for each swath line
        npoints = self.li.shape[0]

        # Create the elevation data array with the first line (baseline)
        self.data = self._get_zi(line, dem, npoints)

        # Simplify baseline
        sline = line.simplify(tolerance=dem.cellsize*5)
        
        # Create the elevation data for the Swath
        for n in range(n_lines):
            dist = line_distance * (n+1)
            left_line = sline.parallel_offset(dist, side="left")
            right_line = sline.parallel_offset(dist, side="right") 
            # Sometimes parallel_offset produces MultiLineStrings ÃÂ¿??
            if left_line.type == "MultiLineString":
                left_line = self._combine_multilines(left_line)
            if right_line.type == "MultiLineString":
                right_line = self._combine_multilines(right_line)

            right_line = self._flip(right_line)
            l_elev = self._get_zi(left_line, dem, npoints)
            r_elev = self._get_zi(right_line, dem, npoints)
            self.data = np.append(self.data, r_elev, axis=1)
            self.data = np.append(self.data, l_elev, axis=1)

        # Get parameters (max, min, mean, q1, q3, HI, relief)
        self.maxz = np.nanmax(self.data, axis=1)
        self.minz = np.nanmin(self.data, axis=1)
        self.meanz = np.nanmean(self.data, axis=1)
        self.q1 = np.nanpercentile(self.data, q=25, axis=1)
        self.q3 = np.nanpercentile(self.data, q=75, axis=1)
        self.HI = (self.meanz - self.minz) / (self.maxz - self.minz)
        self.relief = self.maxz - self.minz
        
        # Get a background polygon for the data
        xi = np.append(self.li, self.li[::-1])
        yi = np.append(self.maxz, self.minz[::-1])
        xi = xi.reshape((xi.size, 1))
        yi = yi.reshape((yi.size, 1))
        self.bg_dat = np.append(xi, yi, axis=1)
        
        # Length of the swath
        self.length = self.li[-1]
        
    def _get_zi(self, line, dem, npoints):
        """
        Get elevations along a line in npoints equally spaced. If any point of the line falls
        outside the DEM or in a NoData cell, a np.nan value will be asigned.
        :param line : Shapely.LineString object. Input LineString
        :param dem : pRaster object. DEM with elevatations.
        :param npoints : int. Number of points along the line to get elevations
        :return zi : Numpy.ndarray. Array with size (npoints, 1) with elevations
        """
        step_size = 1.0/npoints
        zi = []
        for idx in range(npoints):
            pt = line.interpolate(step_size * idx, normalized=True)
            xy = list(pt.coords)[0]
            z = dem.get_xy_value(xy)
            if z == dem.nodata or not z:
                z = np.nan
            zi.append(z)

        return np.array(zi, dtype="float").reshape((len(zi), 1))
 
    def _flip(self, line):
        """
        Flips a LineString object. Returns the new line flipped
        :param line : Shapely.LineString object. Input LineString
        :return line : Shapely.LineString object. Fliped LineString
        """
        coords = list(line.coords)
        coords = np.array(coords)[::-1]
        return LineString(coords)
       
    def _combine_multilines(self, line):
        """
        Combines all the parts of a MultiLineString in a single LineString
        :param line : Shapely.LineString object. Input MultiLineString
        :return line : Shapely.LineString object. Ouput LineString
        """
        xyarr = np.array([], dtype="float32").reshape((0, 2))
        for n in range(len(line.geoms)):
            xyarr = np.append(xyarr, np.array(line.geoms[n].coords), axis=0)
        return LineString(xyarr)
            
    def draw_swath(self, ax, legend=False, drawdata=False, drawbg=False, q=False, **kwargs):
        """
        Draw the swat profile in an Axe
        :param ax : Axe where the profile will be painted. Its cleared before drawing
        :param legend : boolean. Draw the legend
        :param drawdata: boolean. Draw the data (all the profiles)
        :param drawbg: boolean. Draw a background instead of the data
        :param p : boolean. Draw the q1, q3 quartiles
        :kwargs : Colors and line widths for the profiles. See colors and linew dictionaries
        """
        ax.clear()
        colors = {"data": (0.75, 0.75, 0.75),
                  "min": (202./255, 111./255, 30./255),
                  "max": (169./255, 50./255, 38./255),
                  "mean": (22./255, 160./255, 133./255),
                  "q1": (0./255, 191./255, 255./255),
                  "q3": (0./255, 191./255, 255./255)}
        
        colors.update(kwargs)
        linew = {"dataw":0.65, "linesw":1.}
        linew.update(kwargs)

        if drawbg:
            poly = mpatches.Polygon(self.bg_dat, facecolor="0.85")
            ax.add_patch(poly) 
            
        if drawdata:
            for n in range(self.data.shape[1]):
                ax.plot(self.li, self.data[:, n], lw=linew["dataw"], color=colors["data"])

        ax.plot(self.li, self.maxz, lw=linew["linesw"], color=colors["max"], label="max")
        ax.plot(self.li, self.minz, lw=linew["linesw"], color=colors["min"], label="min")
        ax.plot(self.li, self.meanz, lw=linew["linesw"], color=colors["mean"], label="mean")
        
        if q:
            ax.plot(self.li, self.q1, lw=linew["linesw"], color=colors["q1"], label="q1")
            ax.plot(self.li, self.q3, lw=linew["linesw"], color=colors["q3"], label="q3")

        ax.set_xlabel("Distance [m]")
        ax.set_ylabel("Elevation [m]")
        
        ax.set_title(self.name)
        
        # QGIS Adjustment (to make the graphic nicer)
        dz = (self.maxz.max() - self.minz.min()) * 0.05
        ax.set_xlim(0, self.length)
        ax.set_ylim(self.minz.min() - dz, self.maxz.max()  + dz)
        
        if legend:
            legend = ax.legend()
            for tx in legend.texts:
                tx.set_fontsize(12)

    def draw_thi(self, ax, enhanced=False):
        """
        Draws the THI profile in an input Axe
        
        :param ax : matplotlib.Axe object to draw the THI profile
        :param enhanced : boolean. Specify if the enhanced THI (THI*) is calculated
        """
        if enhanced:
            hi = (self.HI - 0.2) / 0.6
        else:
            hi = self.HI

        max_relief = float(np.nanmax(self.relief))
        wi = 0.2 * np.log(self.relief / max_relief) + 1
        thi = (hi - 0.5) * wi + 0.5
        
        ax.plot(self.li, thi, c="k", linewidth=1.2)

        ax.plot([0, self.length], [0.5, 0.5], linestyle="--",
                linewidth=0.75, color=(0.4, 0.4, 0.4))

        ax.set_ylim((0.0, 1.0))
        ax.set_xlim((0.0, self.length))
        ax.set_xlabel("Distance [m]")
        
        if enhanced:
            label = "THI*"
        else:
            label = "THI"
        
        ax.set_ylabel(label)
        ax.set_yticks((0.0, 0.5, 1.0))


class SwathGraph:
    def __init__(self, swaths, fig):
        """
        Class to represent SwathProfiles
        
        :param swaths : list. List with SwathProfile objects
        :param fig : matplotlib.Figure. Figure to represent the swath profiles
        """
        # Create graphic swaths and select the active one
        self.swaths = swaths
        self.nswaths = len(swaths)
        self.id = 0
        self.active_swath = self.swaths[self.id]
        
        # Configure figure and axes
        self.figure = fig
        self.ax1 = None
        self.ax2 = None
        
        # Options
        self.drawTHI = False
        self.legend = True
        self.drawdata = False
        self.drawbg = True
        self.drawq = False
        self.enhanced = False
        
        # Connect to keyboard events
        self.cid = self.figure.canvas.mpl_connect("key_press_event", self.key_press)
        
        # Call the draw function
        self.draw()

    def activate(self):
        """
        Activates the Graphic by connecting to events
        """
        self.cid = self.figure.canvas.mpl_connect("key_press_event", self.key_press)

    def draw(self):
        """
        Function that draws the active SwathProfile with selected parameters
        """
        # Clear the figure and select the active profile
        self.figure.clear()
        self.active_swath = self.swaths[self.id]
        
        # Show one Axe or two if drawTHI is set to True
        if self.drawTHI:
            ax1 = self.figure.add_axes((0.1, 0.25, 0.85, 0.65))
            ax2 = self.figure.add_axes((0.1, 0.10, 0.85, 0.15))
            ax1.set_xticks(())
            ax1.set_xticklabels(())
            ax2.set_yticks((0.0, 0.5, 1.0))
            ax2.yaxis.tick_right()
        else:
            ax1 = self.figure.add_axes((0.1, 0.1, 0.85, 0.80))
            
        # Draw the swath profile with selected options
        self.active_swath.draw_swath(ax1, legend=self.legend, drawdata=self.drawdata,
                                     drawbg = self.drawbg, q=self.drawq)
        
        # Draw the THI profile if drawTHI is set to True
        if self.drawTHI:
            self.active_swath.draw_thi(ax2, self.enhanced)
        
        # Force canvas redraw
        self.figure.canvas.draw()

    def key_press(self, event):
        """
        Interactive key events for hypsometric curves:
            LEFT / RIGHT >> Draws the previous/next swath profile
            L >> Displays/hides the legend
            T >> Displays/hides the THI Graphic below swath graphic
            Q >> Displays/hides the q1 and q3 profiles
            E >> Displays/hides the enhanced THI index
            D >> Draws all the profile data (it can be slow)
            B >> Draws a gray background representing the profile data
    
        """
        if event.key.upper() == "LEFT" or event.key.upper() == "RIGHT":
            self.id += DIRECTIONS[event.key.upper()]
            self.id = self.id % self.nswaths
            self.draw()
            
        elif event.key.upper() == "L":
            if self.legend == False:
                self.legend = True
            else:
                self.legend = False
            self.draw()
            
        elif event.key.upper() == "T":
            if self.drawTHI == False:
                self.drawTHI = True
            else:
                self.drawTHI = False
            self.draw()
        
        elif event.key.upper() == "Q":
            if self.drawq == False:
                self.drawq = True
            else:
                self.drawq = False
            self.draw()
        
        elif event.key.upper() == "E" and self.drawTHI:
            if self.enhanced == False:
                self.enhanced = True
            else:
                self.enhanced = False
            self.draw()
            
        elif event.key.upper() == "D":
            if self.drawdata == False:
                self.drawdata = True
            else:
                self.drawdata = False
            self.draw()
        
        elif event.key.upper() == "B":
            self.drawdata = False
            if self.drawbg == False:
                self.drawbg = True
            else:
                self.drawbg = False
            self.draw()

# praster.py
# ==========
NTYPES = {'int8': 3, 'int16': 3, 'int32': 5, 'int64': 5, 'uint8': 1, 'uint16': 2,
          'uint32': 4, 'uint64': 4, 'float16': 6, 'float32': 6, 'float64': 7}

GTYPES = {1: 'uint8', 2: 'uint16', 3: 'int16', 4: 'uint32', 5: 'int32', 6: 'float32', 7: 'float64'}

def open_raster(raster_path):
    """
    This function open a raster and returns a pRaster instance
    :param raster_path: [str] Path to the raster to open
    :return: pRaster class instance
    """
    raster = gdal.Open(raster_path)
    if not raster:
        return
    array = raster.GetRasterBand(1).ReadAsArray()
    geot = raster.GetGeoTransform()
    proj = raster.GetProjection()
    nodata = raster.GetRasterBand(1).GetNoDataValue()

    return PRaster(array, geot, proj, nodata)


def create_raster(xsize, ysize, dtype=gdal.GDT_Int16, proj="", geot=(0.0, 1.0, 0.0, 0.0, 0.0, 1.0), nodata=0.0):
    """
    This function creates a pRaster object "In Memory", to save it use the method Save(path)
    :param xsize: *int* -- Number of columns of the raster
    :param ysize: *int* -- Number of rows of the raster
    :param dtype: *gdal.GDT type* -- Data type of the new raster (Default = gdal.GDT_Int16)
    :param proj:  *str* -- Projection of the new raster in wkt (Default = "")
    :param geot:  *tuple* -- Geotransform matrix for the new raster (Default = (0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
    :param nodata: *float* -- Nodata value for the new raster (Default = 0.0)
    :return: pRaster instance
    """

    # Creates an empty array and fill up with nodata values
    arrdata = np.empty((ysize, xsize)).astype(GTYPES[dtype])
    arrdata.fill(nodata)

    return PRaster(arrdata, proj, geot, nodata)


def create_from_template(template, dtype=None, nodata=None):
    """
    This function creates a raster with the same parameters than the template. The created raster is "In Memory", to
    save it use the method Save(path)
    :param template: *str* -- Path to the raster template
    :param dtype:  *gdal.GDT type* -- Data type for the new raster (Default = None -> Takes dtype from template)
    :param nodata: *float* -- NoData value for the new raster (Default = None -> Takes nodata from template)
    :return: pRaster instance
    """
    temp_raster = gdal.Open(template)
    temp_banda = temp_raster.GetRasterBand(1)
    geot = temp_raster.GetGeoTransform()
    proj = temp_raster.GetProjection()
    xsize = temp_banda.XSize
    ysize = temp_banda.YSize

    if dtype is None:
        dtype = temp_banda.DataType

    if nodata is None:
        nodata = temp_banda.GetNoDataValue()

    arrdata = np.empty((ysize, xsize)).astype(GTYPES[dtype])
    arrdata.fill(nodata)

    return PRaster(arrdata, geot, proj, nodata)


class PRaster:

    def __init__(self, array, geot=(0.0, 1.0, 0.0, 0.0, 0.0, -1.0), proj="", nodata=None):
        """
        Class to manipulate Raster objects
        :param array: *numpy.ndarray* -- Numpy array with raster data
        :param geot:  *tuple* -- Geotramsformation Matrix (upX, xcell, 0, upY, 0, ycell) (Default = (0,1,0,0,0,-1))
        :param proj:  *str* -- Projection of the new raster in wkt (Default = "")
        :param nodata: *float* -- Nodata value for the new raster (Default = None)
        """
        self.geot = geot
        self.proj = proj        
        self.cellsize = geot[1]
        self.array = array
        self.nodata = nodata
        self.YSize = array.shape[0]
        self.XSize = array.shape[1]       
        self.XMin = geot[0]
        self.YMin = geot[3] - self.cellsize * self.YSize
        self.XMax = geot[0] + self.cellsize * self.XSize
        self.YMax = geot[3]
        
    def get_cell_value(self, cell):
        """
        Get the raster value at a cell location
        :param cell: *tuple* -- (row, col) Cell location
        :return: Value of raster in the cell location
        """
        if cell[0] < 0 or cell[1] < 0:
            return None
        elif cell[0] >= self.YSize or cell[1] >= self.XSize:
            return None
        else:
            return self.array[cell[0], cell[1]]
        
    def get_xy_value(self, point):
        """
        Get the raster value at a point location
        :param point: *tuple* -- (X, Y) Point location
        :return: Value of raster in the point location
        """
        cell = self.xy_2_cell(point)
        return self.get_cell_value(cell)
 
    def xy_2_cell(self, point):
        """
        Get the cell position (row, col) of a point
        :param point: *tuple* -- (X, Y) Point location
        :return: Cell position (row, col)
        """
        row = int((self.YMax - point[1]) / self.cellsize)
        col = int((point[0] - self.XMin) / self.cellsize)
        return row, col
    
    def cell_2_xy(self, cell):
        """
        Get the XY position (X, Y) of a raster cell
        :param cell: *tuple* -- (row, col) Cell location
        :return: XY position (X, Y)
        """
        x = self.XMin + self.cellsize * cell[1] + self.cellsize / 2.
        y = self.YMax - self.cellsize * cell[0] - self.cellsize / 2.
        return x, y
        
    def set_cell_value(self, cell, value):
        """
        Set the raster value at a cell location
        :param cell: *tuple* -- (row, col) Cell location
        :param value: *number* -- New value for the pRaster at cell location
        """
        self.array[cell[0], cell[1]] = value

    def get_window(self, cell, ncells):
        """
        Get a ncells x ncells window around an specific cell. Function includes edge treatment.
        :param cell: *tuple* -- (row, col) Cell location
        :param ncells: *int* -- Number of cells around cell (Full window size will be 2*ncells + 1)
        :return: numpy.ndarray with the data of the raster within the window
        """
        row0 = cell[0] - ncells
        col0 = cell[1] - ncells
        
        nrows = ncells * 2 + 1
        ncols = ncells * 2 + 1
        
        # Check boundary conditions
        if row0 < 0:
            nrows += row0
            row0 = 0   
        if col0 < 0:
            ncols += col0
            col0 = 0
            
        if row0 + nrows >= self.array.shape[0]:
            nrows = self.array.shape[0] - row0
        if col0 + ncols >= self.array.shape[1]:
            ncols = self.array.shape[1] - col0 
        
        # Get raster data in the window and return an array
        window = self.array[row0:row0+nrows, col0:col0+ncols]
        
        if window.size > 0:
            return window
        else:
            return None
        
    def get_flow(self, cell):
        """
        Get the next cell where the flow goes (in case pRaster is a Flow Accumulation raster)
        :param cell: *tuple* -- (row, col) Cell location
        :return: cell (row, col) where the flow goes. It returns None at the edges of the raster
        """

        vec_adyacentes = [(-1, 0), (0, -1), (0, 1), (1, 0)]
        vec_diagonales = [(-1, -1), (-1, 1), (1, -1), (1, 1)]
        
        # Suponemos que el valor mÃÂ¡ximo es el mismo
        cell_value = self.get_cell_value(cell)
        
        max_value = cell_value
        max_pos = cell
    
        # La celda a la que va el flujo no tiene porque ser la de mayor valor de flow accumulation
        # En el caso de que el flujo haga una L la mÃÂ¡xima es la diagonal, pero el flujo va a la adyacente
        # Por ello primero se comprueban las celdas adyacentes y luego las diagonales
    
        for n in vec_adyacentes:
            row = cell[0] + n[0]
            col = cell[1] + n[1]
            if (row < self.YSize and row >= 0) and (col < self.XSize and col >= 0):
                # A veces en raster UInt16 o UInt32 nodata puede ser un valor muy alto
                if self.get_cell_value((row, col)) > max_value:
                    max_value = self.get_cell_value((row, col))
                    max_pos = (row, col)
        
        if cell_value == max_value:
            # Si no hay ninguna celda adyacente con un valor mayor de f
            for n in vec_diagonales:
                row = cell[0] + n[0]
                col = cell[1] + n[1]
                if (row < self.YSize and row >= 0) and (col < self.XSize and col >= 0):
                    # A veces en raster UInt16 o UInt32 nodata puede ser un valor muy alto
                    if self.get_cell_value((row, col)) > max_value and max_value != self.nodata:
                        max_value = self.get_cell_value((row, col))
                        max_pos = (row, col)            
        
        if max_value == cell_value or max_value == self.nodata:
            return None
        else:
            return max_pos
            
    def save_raster(self, path):
        """
        Saves the pRaster in the disk
        :param path: *str* -- Path where new raster will be saved
        """

        if str(self.array.dtype) not in NTYPES.keys():
            return
        else:
            tipo = NTYPES[str(self.array.dtype)]

        driver = gdal.GetDriverByName("GTiff")
        raster = driver.Create(path, self.XSize, self.YSize, 1, tipo)
        raster.SetGeoTransform(self.geot)
        raster.SetProjection(self.proj)
        if self.nodata:
            raster.GetRasterBand(1).SetNoDataValue(self.nodata)
        raster.GetRasterBand(1).WriteArray(self.array)


# MAIN PROGRAM
# =============
def main(line_shp, dem, name_field, width, nlines, step, fullres):   
    # Read line shapefile and get names
    centerlines, names = read_line_shapefile(line_shp, name_field)

    # Open dem
    dem_raster = open_raster(dem)
    if not dem_raster:
        print("It was a problem opening the DEM")
        return
    
    # If nlines and step are 0, get default values
    if nlines == 0:
        nlines = 50
    if step == 0:
        step = dem_raster.cellsize
    
    # If fullresolution, take maximum resolution and discard previous values
    if fullres:
        nlines = int(width/dem_raster.cellsize)
        step = dem_raster.cellsize
        
    # Log Messages
    progress.setText("SwathProfiler started")
    progress.setText("Found " + str(len(centerlines)) + " lines in the shapefile")
    progress.setText("")
    progress.setText("Parameters:")
    progress.setText("Width: " + str(width) + " m" )
    progress.setText("Number of lines: " + str(nlines))
    progress.setText("Step size: " + str(step))
    progress.setText("")
    
    # Create the swath profiles list
    swath_profiles = []
    
    for idx, line in enumerate(centerlines):
        
        # Log Messages
        progress.setText("Processing {0} of {1} lines".format(idx+1, len(centerlines)))
        
        sline = SwathProfile(line, dem_raster, width=width, n_lines=nlines, 
                                           step_size=step, name=names[idx])      
        if len(sline.data) > 0:
            swath_profiles.append(sline)

    # Create graphic and draw swath profiles
    if len(swath_profiles) > 0:
        fig = plt.figure()
        return SwathGraph(swath_profiles, fig)


my_graph = main(line_shp, dem, name_field, width, nlines, step, fullres)   
plt.show()
