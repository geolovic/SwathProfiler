# -*- coding: utf-8 -*-
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

#  Version: 1.0 
#  November 23, 2016

#  Last modified March 23, 2017

import ogr
import gdal
from shapely.geometry import LineString
from shapely.geometry import Polygon
import numpy as np
import matplotlib.pyplot as plt

# QGIS Toolbox input parameters
# ##Line_Shapefile=vector
# ##Digital_Elevation_Model=raster
# ##Names_field=field Line_Shapefile
# ##Swath_width=number 1000
# ##Number_of_lines=number 0
# ##Step_size=number 0
# ##Draw_data=boolean True
# ##Draw_legend=boolean True

# Debug code (mimic parameters from QGIS toolbox)
Line_Shapefile = "data/perfiles.shp"
Digital_Elevation_Model = "data/coskie.tif"
Names_field = "Name"
Swath_width = 1000
Number_of_lines = 0
Step_size = 0
Draw_data = True
Draw_legend = True


def open_raster(raster_path):
    """
    This function open a raster and returns a pRaster instance

    raster_path :: *str*
      Path to the raster
    """
    raster = gdal.Open(raster_path)
    if not raster:
        return
    array = raster.GetRasterBand(1).ReadAsArray()
    geot = raster.GetGeoTransform()
    proj = raster.GetProjection()
    nodata = raster.GetRasterBand(1).GetNoDataValue()

    return PRaster(array, geot, proj, nodata)


def create(path, xsize, ysize, datatype=gdal.GDT_Int16, projection="",
           geotransform=(0.0, 1.0, 0.0, 0.0, 0.0, 1.0), nodata=None):
    """
    This function creates a pRaster object (in Memory)

    Parameters:
    ==================
    xsize :: *int*
        Number of columns of the raster
    ysize :: *int*
        Number of rows of the raster
    datatype :: *gdal.GDT type (Default=gdal.GDT_Int16)*
        Type of the data for the new raster (default gdal.GDT_Int16)
    projection :: *str (Default="")*
        Projection for the new raster in wkt
    geotransform :: *tuple (Default=(0.0, 1.0, 0.0, 0.0, 0.0, 1.0))*
        Geotransform matrix for the new raster
    nodata :: *value (Default=None)*
        Nodata value for the new raster (default None)
    """
    driver = gdal.GetDriverByName("GTiff")
    raster = driver.Create(path, xsize, ysize, 1, datatype)
    raster.SetGeoTransform(geotransform)
    raster.SetProjection(projection)
    if nodata is not None:
        banda = raster.GetRasterBand(1)
        banda.SetNoDataValue(nodata)

    return open_raster(path)


def create_from_template(path, template, datatype=gdal.GDT_Int16, nodata=None):
    """
    This function creates a raster with the same parameters than the template
    and returns a pRaster object

    Parameters:
    ==================
    path :: *str*
        Full path to the raster to be created
    template :: *str*
        Full path to the raster template
    datatype :: *gdal.GDT type (Default=gdal.GDT_Int16)*
        Pixel datatype
    nodata :: *value (Default=None)*
        Value for nodata
    """
    driver = gdal.GetDriverByName("GTiff")
    temp_raster = gdal.Open(template)
    temp_banda = temp_raster.GetRasterBand(1)
    xsize = temp_banda.XSize
    ysize = temp_banda.YSize
    raster = driver.Create(path, xsize, ysize, 1, datatype)
    if nodata is not None:
        banda = raster.GetRasterBand(1)
        banda.SetNoDataValue(nodata)
    raster.SetGeoTransform(temp_raster.GetGeoTransform())
    raster.SetProjection(temp_raster.GetProjection())

    return open_raster(path)


class PRaster:
    """
    Class that defines by combining a numpy.array with some Raster parameters
    as cellsize, extension, nodata values, etc.

    Parameters:
    ==================
    array :: *numpy.array*
        Array with the data for the raster
    geot :: *tuple (Default: (0.0, 1.0 ,0.0 ,0.0 ,0.0 ,1.0)*
        Tuple with GeoTransform values
    proj :: *str (Default: "")
        String with coordinate system in Well Known Text format
    nodata :: *value (Default: None)*
        Value (int or float) for nodata
    """

    def __init__(self, array, geot=(0.0, 1.0, 0.0, 0.0, 0.0, 1.0), proj="", nodata=None):
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
        This function returns the value of the raster at a cell location.
        cell :: *tuple*
            Tuple of ints (row, col)

        """
        if cell[0] < 0 or cell[1] < 0:
            return None
        elif cell[0] >= self.YSize or cell[1] >= self.XSize:
            return None
        else:
            return self.array[cell[0], cell[1]]

    def get_xy_value(self, point):
        """
        This function returns the value of the raster at a point location.
        point :: *tuple*
            Tuple of floats (x, y)
        """
        cell = self.xy_2_cell(point)
        return self.get_cell_value(cell)

    def xy_2_cell(self, point):
        """
        This function returns the position of a point in row-col coordinates
        point :: *tuple*
            Tuple of floats (x, y)
        """
        row = int((self.YMax - point[1]) / self.cellsize)
        col = int((point[0] - self.XMin) / self.cellsize)
        return row, col

    def cell_2_xy(self, cell):
        """
        This function returns the position of a cell in XY coordinates
        cell :: *tuple*
            Tuple of ints (row, col)
        """
        x = self.XMin + (self.cellsize * cell[1]) + (self.cellsize / 2)
        y = self.YMax - (self.cellsize * cell[0]) - (self.cellsize / 2)
        return x, y

    def set_cell_value(self, cell, value):
        """
        This function writes a value at a cell location. Values are wrote in
        internal raster numpy array.
        value :: *number*
            Value to write in the raster
        cell :: *tuple*
            Tuple of ints (row, col) where the value will write
        """
        self.array[cell[0], cell[1]] = value

    def get_window(self, cell, ncells):
        """
        Esta función obtiene una ventana de NxN alrededor de un punto del raster.
        La función incluye un tratamiento para obtener los valores en los bordes.

        Parametros:
        ===========
        cell :: *tuple*
            Tuple of ints (row, col) que representa el punto central
        ncells :: *int*
          Numero de celdas de la ventana alrededor del punto central

        Retorna:
        ==========
        arr :: *numpy.Array*
          Numpy array con los datos en la ventana de análisis
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
        window = self.array[row0:row0 + nrows, col0:col0 + ncols]

        if window.size > 0:
            return window
        else:
            return None

    def get_flow(self, cell):
        """
        This function returns the next cell where the flow goes
        It suppouses that the current is a flow accumulation raster

        Parameters:
        ==================
        *cell* : tuple
            Tuple of ints (row, col)

        Output:
        ===================
        *flow_pos* : tuple
            Tuple (row, col) of the next cell where the flow goes down stream
        """
        vec_adyacentes = [(-1, 0), (0, -1), (0, 1), (1, 0)]
        vec_diagonales = [(-1, -1), (-1, 1), (1, -1), (1, 1)]

        # Suponemos que el valor máximo es el mismo
        cell_value = self.get_cell_value(cell)

        max_value = cell_value
        max_pos = cell

        # La celda a la que va el flujo no tiene porque ser la de mayor valor de flow accumulation
        # En el caso de que el flujo haga una L la máxima es la diagonal, pero el flujo va a la adyacente
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

    def save(self, path):
        """
        This function writes the content of the internal array into a new raster.
        """
        types = {'int8': 3, 'int16': 3, 'int32': 5, 'int64': 5, 'uint8': 1, 'uint16': 2,
                 'uint32': 4, 'uint64': 4, 'float16': 6, 'float32': 6, 'float64': 7}
        if self.array.dtype not in types.keys():
            return
        else:
            tipo = types[str(self.array.dtype)]

        driver = gdal.GetDriverByName("GTiff")
        raster = driver.Create(path, self.XSize, self.YSize, 1, tipo)
        raster.SetGeoTransform(self.geot)
        raster.SetProjection(self.proj)
        if self.nodata:
            raster.GetRasterBand(1).SetNoDataValue(self.nodata)
        raster.GetRasterBand(1).WriteArray(self.array)


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
        npoints = int(line.length / step_size)

        # Get number of lines. By default width/dem.cellsize
        if n_lines is None:
            n_lines = int(width / dem.cellsize)

        # Get distances and center line (columns 0 and 1)
        self.li = self._get_li(npoints, step_size)
        self.swaths.append(self._get_zi(line, dem, npoints))

        # Get right and left swath lines
        sline = line.simplify(tolerance=dem.cellsize)

        line_distance = float(width) / n_lines
        for n in range(n_lines):
            dist = line_distance * (n + 1)
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
                  "min": (202. / 255, 111. / 255, 30. / 255),
                  "max": (169. / 255, 50. / 255, 38. / 255),
                  "mean": (247. / 255, 220. / 255, 111. / 255),
                  "q1": (22. / 255, 160. / 255, 133. / 255),
                  "q3": (22. / 255, 160. / 255, 133. / 255)}

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
        step_size = 1.0 / npoints
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


def main(line, dem, width, nlines, step, names, drawlegend, drawdata):
    # Read line shapefile and raster
    centerlines, names = read_line_shapefile(line, names)

    dem_raster = open_raster(dem)
    if nlines == 0:
        nlines = None
    if step == 0:
        step = None

    # Create the swath profiles list
    swath_profiles = []
    n = 0
    for centerline in centerlines:
        sline = SwathProfile(centerline, dem_raster, width, nlines, step, name=names[n])
        if len(sline.swaths) > 0:
            swath_profiles.append(sline)
        n += 1

    # Create graphic and draw swath profiles
    if len(swath_profiles) > 0:
        fig = plt.figure()
        swath_graph = SwathGraph(swath_profiles, fig, drawlegend, drawdata)
        plt.show()
        return swath_graph


def read_line_shapefile(shapefile, names_field=""):
    """
    This function reads a line shapefile and returns a tuple (out_lines, out_names)

    Parameters:
    ================
    shapefile :: *str*
      Path to the line shapefile with profile centerlines
    names_field :: *str (Default: ""*
      Name of the field with the profile names. If skipped, profiles will be named sequentially

    Returns:
    ==============
    (out_lines, out_names) :: *tuple*
        out_lines :: List with shapely.geometry.LineString objects representing centerlines
        out_names :: List of string with profile names
    """
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataset = driver.Open(shapefile)
    layer = dataset.GetLayer()

    # layerdef = layer.GetLayerDefn()
    if not layer.GetGeomType() == 2:
        return

    out_names = []
    out_lines = []
    perfil_id = 0
    for feat in layer:
        out_names.append(feat.GetField(names_field))
        geom = feat.GetGeometryRef()
        coords = []
        for n in range(geom.GetPointCount()):
            pt = geom.GetPoint(n)
            coords.append((pt[0], pt[1]))
        out_lines.append(LineString(coords))
        perfil_id += 1

    return out_lines, out_names


sg = main(Line_Shapefile, Digital_Elevation_Model, Swath_width, Number_of_lines,
          Step_size, str(Names_field), Draw_legend, Draw_data)
