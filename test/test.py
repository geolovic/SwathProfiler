from shapely.geometry import LineString
from shapely.geometry import Polygon
import matplotlib.pyplot as plt

def paint_polygon(polygon, ax):
    xs = []
    ys = []
    for x, y in list(polygon.exterior.coords):
        xs.append(x)
        ys.append(y)
    ax.plot(xs, ys)



width = 1000
XMax = 394291.222183
XMin = 387918.803394
YMax = 4457348.45963
YMin = 4452007.67884

line = LineString([(388724.03, 4455466.51), (393151.81,  4455466.51)])
buff = line.buffer(width,cap_style=2)
raster_extent = Polygon([(XMin, YMin), (XMin, YMax), (XMax, YMax), (XMax, YMin), (XMin, YMin)])

print raster_extent.covers(line)
print raster_extent.covers(buff)

fig = plt.figure()
ax = fig.add_subplot(111)
paint_polygon(raster_extent, ax)
paint_polygon(buff, ax)

