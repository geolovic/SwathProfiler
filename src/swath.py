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

#  Version: 1.0
#  November 23, 2016

#  Last modified March 23, 2017

from shapely.geometry import LineString
from shapely.geometry import Polygon
import numpy as np


class SwathProfile:
    def __init__(self, line=None, dem=None, width=0, n_lines=None, step_size=None, name=""):
        """
        Class to create a swath profile and related parameters

        :param line: shapely.geometry.LineString - LineString the swath profile centerline
        :param dem: praster.pRaster - pRaster with the Digital Elevation Model
        :param width: float or int - Half-width of the swath profile (in data units)
        :param n_lines: int - number of lines at each side of the swath centerline
        :param step_size: float - Step-size to get elevation points along the profile
        :param name: str - Name of the profile
        """

        self.swaths = []
        self.name = str(name)

        # Creates an empty SwathProfile Object
        if line is None:
            return
        
        # Check if all the profiles will lie within the DEM raster
        if not self._is_inside(line, width, dem):
            return
        
        # If step_size is None, get dem resolution as step_size
        if step_size is None:
            step_size = dem.cellsize
        
        # Get number of points for each swath line
        npoints = int(line.length/step_size)

        # Get number of lines. By default width/dem.cellsize
        if n_lines is None:
            n_lines = int(width/dem.cellsize)

        # Get distances and center line (columns 0 and 1)
        self.li = self._get_li(npoints, step_size)
        self.swaths.append(self._get_zi(line, dem, npoints))

        # Get right and left swath lines
        sline = line.simplify(tolerance=dem.cellsize)

        line_distance = float(width) / n_lines
        for n in range(n_lines):
            dist = line_distance * (n+1)
            left_line = sline.parallel_offset(dist, side="left")
            right_line = sline.parallel_offset(dist, side="right")
            right_line = self._flip(right_line)
            self.swaths.append(self._get_zi(left_line, dem, npoints))
            self.swaths.append(self._get_zi(right_line, dem, npoints))
        
        self.swaths = np.array(self.swaths).T

        # Get max and min data arrays
        self.maxz = np.nanmax(self.swaths, axis=1)
        self.minz = np.nanmin(self.swaths, axis=1)
        self.meanz = np.nanmean(self.swaths, axis=1)
        self.q1 = np.nanpercentile(self.swaths, q=25, axis=1)
        self.q3 = np.nanpercentile(self.swaths, q=75, axis=1)
        self.HI = (self.meanz - self.minz) / (self.maxz - self.minz)
        self.relief = self.maxz - self.minz
        
        # Length of the swath
        self.length = line.length

    def draw_swath(self, ax, legend=False, drawdata=True, **kwargs):
        
        colors = {"data": (0.8, 0.8, 0.8),
                  "min": (202./255, 111./255, 30./255),
                  "max": (169./255, 50./255, 38./255),
                  "mean": (247./255, 220./255, 111./255),
                  "q1": (22./255, 160./255, 133./255),
                  "q3": (22./255, 160./255, 133./255)}
        
        colors.update(kwargs)

        if drawdata:
            for n in range(self.swaths.shape[1]):
                ax.plot(self.li, self.swaths[:, n], linewidth=0.75, color=colors["data"])

        ax.plot(self.li, self.maxz, linewidth=1.5, color=colors["max"], label="max")
        ax.plot(self.li, self.minz, linewidth=1.5, color=colors["min"], label="min")
        ax.plot(self.li, self.meanz, linewidth=1.5, color=colors["mean"], label="mean")
        ax.plot(self.li, self.q1, linewidth=1.5, color=colors["q1"], label="q1")
        ax.plot(self.li, self.q3, linewidth=1.5, color=colors["q3"], label="q3")

        ax.set_xlabel("Distance [m]")
        ax.set_ylabel("Elevation [m]")
        
        ax.set_title(self.name)

        ax.set_xlim((0.0, self.length))
        ax.set_xticks(())
        
        if legend:
            legend = ax.legend()
            for tx in legend.texts:
                tx.set_fontsize(12)

    def draw_thi(self, ax, enhanced=False):
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

    def _get_zi(self, line, dem, npoints):
        step_size = 1.0/npoints
        zi = []
        for idx in range(npoints):
            pt = line.interpolate(step_size * idx, normalized=True)
            xy = list(pt.coords)[0]
            z = dem.get_xy_value(xy)
            if z == dem.nodata:
                z = np.nan
            zi.append(z)

        return zi

    def _get_li(self, npoints, step_size):
        """
        Get distances of the swath profile
        """
        li = []
        for n in range(npoints):
            li.append(n * step_size)
        return li

    def _flip(self, line):
        """
        Flip a LineString object. Returns a new line flipped
        """
        coords = list(line.coords)
        coords = np.array(coords)[::-1]
        return LineString(coords)
    
    def _is_inside(self, line, width, raster):
        buff_pol = line.buffer(width)
        raster_extent = Polygon([(raster.XMin, raster.YMin), (raster.XMin, raster.YMax), (raster.XMax, raster.YMax),
                                 (raster.XMax, raster.YMin), (raster.XMin, raster.YMin)])
        
        return raster_extent.covers(buff_pol) 


class SwathGraph:
    def __init__(self, swaths, fig, legend=True, drawdata=True):
        # Select active swath
        self.swaths = swaths
        self.nswaths = len(swaths)
        self.id = 0
        self.active_swath = self.swaths[self.id]
        self.cid = None
        # Configure chart
        self.ax1 = fig.add_axes((0.1, 0.25, 0.85, 0.65))
        self.ax2 = fig.add_axes((0.1, 0.10, 0.85, 0.15))
        self.ax1.set_xticks(())
        self.ax2.set_yticks((0.0, 0.5, 1.0))
        self.ax2.yaxis.tick_right()
        # Options
        self.legend = legend
        self.drawdata = drawdata
        # Activate graphic
        self.activate()
        self.draw()

    def activate(self):
        self.cid = self.ax1.figure.canvas.mpl_connect("button_press_event", self.button_press)

    def draw(self):
        self.ax1.clear()
        self.ax2.clear()
        self.active_swath.draw_swath(self.ax1, legend=self.legend, drawdata=self.drawdata)
        self.active_swath.draw_thi(self.ax2, True)
        self.ax1.figure.canvas.draw()

    def button_press(self, event):
        if event.button == 1:
            self.id -= 1
        elif event.button > 1:
            self.id += 1

        if self.id < 0:
            self.id = 0
        if self.id == self.nswaths:
            self.id = self.nswaths - 1

        self.active_swath = self.swaths[self.id]
        self.draw()
