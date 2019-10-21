import matplotlib.pyplot as plt

def getlists(p):
    xl = []
    yl = []
    for (x,y) in p:
        xl.append(x)
        yl.append(y)
    return xl,yl

def convex_hull(points):
    """Computes the convex hull of a set of 2D points.

    Input: an iterable sequence of (x, y) pairs representing the points.
    Output: a list of vertices of the convex hull in counter-clockwise order,
      starting from the vertex with the lexicographically smallest coordinates.
    Implements Andrew's monotone chain algorithm. O(n log n) complexity.
    """

    # Sort the points lexicographically (tuples are compared lexicographically).
    # Remove duplicates to detect the case we have just one unique point.
    points = sorted(set(points))

    # Boring case: no points or a single point, possibly repeated multiple times.
    if len(points) <= 1:
        return points

    # 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
    # Returns a positive value, if OAB makes a counter-clockwise turn,
    # negative for clockwise turn, and zero if the points are collinear.
    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    # Build lower hull
    lower = []
    for p in points:
        cont = 1
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            print("antes   "), print(cont), print(lower)
            lower.pop()
            print("despues     "),print(lower)
            cont += 1
        lower.append(p)
    xlower ,ylower = getlists(lower)
    plt.plot(xlower,ylower,color="yellow")
    # Build upper hull
    upper = []
    for p in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    print(upper)
    print("hello2         ")
    print(cross((2,0),(2,4),(2.5,3)))

    xupper ,yupper = getlists(upper)
    plt.plot(xupper,yupper,color="blue")


    return lower[:-1] + upper[:-1]

# The worst case, cause we add all points, but when arribe tot te prev last, we see that there is no
points = [(0,0),(1,1),(2,2.1),(3,3.4),(4,5),(5,3),(6,10),(0,0)]
xpoint, ypoint = getlists(points)

plt.plot(xpoint,ypoint, color='black',marker='o',markersize=5)
p = convex_hull(points)
print(p)

plt.show()
