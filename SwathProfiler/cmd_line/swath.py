# -*- coding: utf-8 -*-
#
#  swath.py
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
#  November 12, 2017

#  Last modified November 14, 2017

import ogr
from shapely.geometry import LineString
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

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
        
        # Get number of lines (By default 50)
        if n_lines is None:
            n_lines = 50
        elif n_lines > int(width/dem.cellsize):
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
            # Sometimes parallel_offset produces MultiLineStrings ¿??
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
        ax.set_xlim((0.0, self.length))
        
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