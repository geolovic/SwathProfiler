# -*- coding: iso-8859-15 -*-
"""
José Vicente Pérez
Granada University (Spain)
November 13, 2016 
  
Copyright (C) 2016  J. Vicente Pérez, Universidad de Granada
Unpublished material - All rights reserved under the Copyright Laws of Spain

For additional information, contact to:
José Vicente Pérez Peña
Dpto. Geodinámica-Universidad de Granada
18071 Granada, Spain
vperez@ugr.es // geolovic@gmail.com

Version: 3.0 
Last modified December 03, 2016
"""

import gdal


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
        window = self.array[row0:row0+nrows, col0:col0+ncols]
        
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
