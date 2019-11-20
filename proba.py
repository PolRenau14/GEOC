
from matplotlib import pyplot as plt
from shapely.geometry.polygon import LinearRing, Polygon



# 1: valid ring

#poly = Polygon([(0, 0), (0, 2), (1, 1),
#                (2, 2), (2, 0), (1, 0.8), (0, 0)])
#x,y = poly.exterior.xy

points = [(0, 0), (0.5, 7), (1, 7.5),(2.1,6.8), (3.5, 7.25), (5,3),(3.5,3),(3.7,3.3),(4.1,3.21),(2,4.5),(1.6,3.8),(2.1,4),(1.9,2),(5,2),(5,1.5),(1,1.5),(0.6,1.7),(0.5,1),(6,1),(6,3),(7,3), (7,0),(0, 0)]
ring = LinearRing(points)
x, y = ring.xy

fig = plt.figure(1, figsize=(15,15), dpi=90)
ax = fig.add_subplot(111)
ax.plot(x, y)
ax.scatter(x,y,s=20) #draw all vertex.
#select reflex vertex and drawthem.
reflex = [points[3],points[8],points[9],points[10],points[12],points[16],points[17],points[18]]
reflexX = []
reflexY =[]
for (w,s) in reflex:
    reflexX.append(w)
    reflexY.append(s)
ax.scatter(reflexX,reflexY,s=20,color="red")
ax.set_title('Polygon Edges')
plt.savefig('1.png')
plt.close()

#Pintem el nou poligon.
aux = []
auy = []
list = [3,8,9,10,12,16,17,18]
for i in range(0,len(x)):
    if i not in list:
        aux.append(x[i])
        auy.append(y[i])
fig = plt.figure(1, figsize=(15,15), dpi=90)
ax = fig.add_subplot(111)
ax.plot(aux, auy)
ax.scatter(aux,auy,s=20) #draw all vertex.
reflexX = [aux[6:8]]
reflexY = [auy[6:8]]
ax.scatter(reflexX,reflexY,s=20,color="red")
plt.savefig('2.png')
plt.close()

x = []
y = []
list = [5,6,7]
for i in range(0,len(aux)):
    if i not in list:
        x.append(aux[i])
        y.append(auy[i])
fig = plt.figure(1, figsize=(15,15), dpi=90)
ax = fig.add_subplot(111)
ax.plot(x, y)
ax.scatter(x,y,s=20) #draw all vertex.
plt.show()
