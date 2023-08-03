#!/usr/bin/env python
# coding: utf-8

# In[1353]:


import matplotlib.pyplot as plt
import pulp
import math
import random
import numpy as np
import simpy


# Points and Distances

# In[1354]:


resolution = 30000  # city of about 10km diameter


# In[1355]:


def dist(p1, p2):
    (x1, y1) = p1
    (x2, y2) = p2
    return int(math.sqrt((x1-x2)**2+(y1-y2)**2))


# In[1356]:


def nearest(X, P):
    minD = math.inf
    minP = None
    for p in P:
        d=dist(X, p)
        if d<minD:
            minP, minD = p, d
    return minP    


# In[1357]:


def generatePoints(n):
    
    def gen():
        x = int(random.gauss(mu=resolution/2, sigma=resolution/8))
        return min(resolution-50, max(50, x))
    
    mindist = resolution/(2*math.sqrt(n)) 
    # avoid points to neart to each other
    P = []
    while len(P)<n:
        i=len(P)
    
        x0=gen()
        y0=gen()
        
        # don't place new points near existing points
        P.append((x0,y0))
        for j in range(0,i):
            if dist(P[i],P[j])<mindist:
                P=P[:-1]
                break
    return P


# Lines
# 
# Lines are represented as pairs of Points, but there is no preferred order of points

# In[1358]:


def equalLines(l1, l2):
    if l1==l2:
        return True
    else:
        return l1[0]==l2[1] and l1[1]==l2[0]


# In[1359]:


def rev(L):
    return L[1],L[0]


# In[1360]:


def solve(x11, x12, y1, x21, x22, y2):
    
    def Det(x11, x12, x21, x22):
        return x11*x22-x21*x12

    D = Det(x11, x12, x21, x22)
    Ds = Det(y1, x12, y2, x22)
    Dt = Det(x11, y1, x21, y2)
    if D==0:
        return False
    s=Ds/D
    t=Dt/D
    return 0 <= s and s <= 1 and 0 <= t and t <= 1
    
def intersecting(l1, l2):
    p1, p2 = l1
    q1, q2 = l2
    if p1==q1 or p1==q2 or p2==q1 or p2==q2:
        return False
    xp1, yp1 = p1
    xp2, yp2 = p2
    xq1, yq1 = q1
    xq2, yq2 = q2
    return solve(xp2-xp1, xq1-xq2, xq1-xp1,
                 yp2-yp1, yq1-yq2, yq1-yp1)


# Triangles

# In[1361]:


def equalTriangles(t1, t2):
    P1, P2, P3 = t1
    Q1, Q2, Q3 = t2
    if P1==Q1:
        if P2==Q2:
            return P3==Q3
        elif P2==Q3:
            return P3==Q2
        else:
            return False
    elif P1==Q2:
        if P2==Q1:
            return P3==Q3
        elif P2==Q3:
            return P3==Q1
        else:
            return False
    elif P1==Q3:
        if P2==Q1:
            return P3==Q2
        elif P2==Q2:
            return P3==Q1
        else:
            return False
    else:
        return False


# In[1362]:


def removeTriangle(t, T):
    for tt in T:
        if equalTriangles(t, tt):
            T.remove(tt)
            return True
    return False    


# In[1363]:


def addTriangle(t, T):
    for tt in T:
        if equalTriangles(t, tt):
            return
    T.append(t)  


# In[1364]:


def area(p1, p2, p3):
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3
    return abs(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))/2


# In[1365]:


def insideTriangle(x, t):
    p, q, r = t
    if x==p or x==q or x==r:
        return False
    return abs(area(p,q,x) + area(p,r,x)+ area(q,r,x) - area(p,q,r))<0.00001 


# In[1366]:


def sides(t):
    A, B, C = t
    return [(B, C), (C, A), (A, B)]


# In[1367]:


def longestSide(t):
    return sorted(sides(t), key=lambda s: dist(s[0], s[1]), reverse=True)[0]


# In[1368]:


def commonSide(t1, t2):
    S1 = sides(t1)
    S2 = sides(t2)
    for s1 in S1:
        for s2 in S2:
            if equalLines(s1, s2):
                return s1
    return None


# In[1369]:


def oppositePoint(t, l):
    A, B, C = t
    P, Q = l
    if A==P:
        return C if B==Q else B
    elif A==Q:
        return C if B==P else B
    else:
        return A


# In[1370]:


def intersectingTriangles(c1, c2):
    c= commonSide(c1, c2)
    if c is None:
        return False
    A, B = c
    C1 = oppositePoint(c1, c)
    C2 = oppositePoint(c2, c)
    return intersecting((A, C1), (B,C2)) or             intersecting((A, C2), (B, C1))

def intersectingLines(c1, c2):
    c= commonSide(c1, c2)
    if c is None:
        return None
    A, B = c
    C1 = oppositePoint(c1, c)
    C2 = oppositePoint(c2, c)
    if intersecting((A, C1), (B,C2)):
        return (A, C1), (B,C2)
    elif intersecting((A, C2), (B, C1)):
        return (A, C2), (B, C1)
    else:
        return None


# In[1371]:


def defuse(t1, t2):
    # is only called if intersectingTriangles(t1, t2)
    l1, l2 = intersectingLines(t1, t2)
    c = commonSide(t1, t2)
    if c is None:
        return None
    A, B = c
    C1 = oppositePoint(t1, c)
    C2 = oppositePoint(t2, c)
    if intersecting((A, C1), (B,C2)):
        return [C1, C2, A]
    elif intersecting((A, C2), (B, C1)):
        return [C1, C2, B]


# In[1372]:


def slimTriangle(t):
    A, B, C = t
    a = dist(B,C)
    b = dist(A,C)
    c = dist(A,B)
    [a, b, c] = sorted([a, b, c])
    return a+b<1.07*c


# In[1373]:


def plotTriangle(t, style='r-o', lw=1, ms=3):
    p1, p2, p3 = t
    plt.plot( [ p1[0], p2[0], p3[0], p1[0] ],
              [ p1[1], p2[1], p3[1], p1[1] ], 
              style, lw=lw, ms=ms)

def plotTriangles(T, style='r-o', lw=1, ms=3):
    plt.xlim(0,resolution)
    plt.ylim(0,resolution)
    plt.axis('off')
    for t in T:
        plotTriangle(t, style, lw=lw, ms=ms)
    plt.show()


# Triangulation
# 
# The triangulations algorithm is an only slighty modified version as documented by: 
# 
# Agryzkov T., Oliver J.L., Tortosa L., Vicent J.F. (2014) A Method to Triangulate a Set of Points in the Plane. In: Murgante B. et al. (eds) Computational Science and Its Applications – ICCSA 2014. ICCSA 2014. Lecture Notes in Computer Science, vol 8580. Springer, Cham DOI 10.1007/978-3-319-09129-7_25

# In[1374]:


def triangulation(P, show=False):
    x = [ p[0] for p in P ]
    y = [ p[1] for p in P ]
    minx = min(x)
    maxx = max(x)
    miny = min(y)
    maxy = max(y)
    nmaxx = 1.1 * maxx
    nmaxy = 1.1 * maxy
    nminx = max(0.5 * minx, 0)
    nminy = max(0.5 * miny, 0)
    distx = nmaxx-nminx
    disty = nmaxy-nminy
    d = math.sqrt(distx*disty/(2*len(P)))
    dx = math.ceil(distx/d)
    dy = math.ceil(disty/d)

        
    def kx(i):
        return nminx+i*d
    
    def ky(j):
        return nminy+j*d
    
    def k(i, j):
        return (kx(i), ky(j))
    
    # def n(i, j):
    #     return nearest( (kx(i), ky(j)), P)
    def n(i, j):
        X = (kx(i), ky(j))
        minD = math.inf
        minP = None
        for p in P:
            d=dist(X, p)
            if d<minD:
                minP, minD = p, d
        return minP  

    if show:

        plt.figure(0)
        plt.xlim(0,resolution)
        plt.ylim(0,resolution)
        plt.axis('off')
    
        plt.plot([ p[0] for p in P ], 
                 [ p[1] for p in P ], 'ro', ms=10)    
    
        for i in range(0, dx+1):
            for j in range(0, dy):
                plt.plot( [ kx(i), kx(i) ],
                          [ ky(j), ky(j+1) ], 'b:.')
            
        for j in range(0, dy+1):
            for i in range(0, dx):
                plt.plot( [ kx(i), kx(i+1) ],
                          [ ky(j), ky(j) ], 'b:.')
            
        for i in range(0, dx):
            for j in range(0, dy):
                plt.plot( [ kx(i), kx(i+1) ],
                          [ ky(j), ky(j+1) ], 'b:.')
    
        for i in range(0, dx+1):
            for j in range(0, dy+1):
                nx, ny = n(i,j)
                plt.plot( [ kx(i), nx ],
                          [ ky(j), ny ], 'g:')


    
    # set of triangles
    T = []
    for i in range(0, dx):
        for j in range(0, dy):
            p1=n(i,j)
            p2=n(i+1, j)
            p3=n(i, j+1)
            p4=n(i+1, j+1)
            if p1!=p2 and p2!=p4 and p4!=p1:
                addTriangle([p1, p2, p4], T)
            if p1!=p3 and p3!=p4 and p4!=p1:
                addTriangle([p1, p3, p4], T)
                  
    # Replace a triangle containing an inner point
    # with three triangles formed from this inner point
    C = T.copy()
    for p in P:
        for t in T:
            if insideTriangle(p, t):
                p1, p2, p3 = t
                C.remove(t)
                C.append( (p1, p2, p) )
                C.append( (p2, p3, p) )
                C.append( (p3, p1, p) )
         
    PLOT = 0
    
    def showTriangulation(T, style='r-o', check=False):
        nonlocal PLOT
        PLOT += 1
        plt.figure(PLOT)
        plt.xlim(0,resolution)
        plt.ylim(0,resolution)
        plt.axis('off')
        for t in T:
            plotTriangle(t, 'r-o', lw=0.5)
        
        if check:
            found=False
            for t1 in T:
                for t2 in T:
                    if t1!=t2 and intersectingTriangles(t1, t2):
                        plotTriangle(t1, 'b-o')
                        plotTriangle(t2, 'g-o')
                        found=True
                        break
                if found:
                    break
  
    found = True
    while found:
        
        if show:
            showTriangulation(C)
        
        D=C.copy()
        found = False
        for c1 in C:
            for c2 in C:
                if c1!=c2 and intersectingTriangles(c1, c2):
                    # print(f"replacing")
                    # print(triangle_String(c1))
                    # print(triangle_String(c2))
                    c3=defuse(c1, c2)
                    removeTriangle(c2, D)
                    # print(triangle_String(c3))
                    addTriangle(c3, D)
                    found=True
                    break
            if found:
                break
        C=D
        
    # Check for isolated points
    # This can now probably not happen anymore, 
    # but just leave it in for now
    
    def isolated(p, C):
        for c in C:
            p1, p2, p3 = c
            if p==p1 or p==p2 or p==p3:
                return False
        return True
    
    for p in P:
        if isolated(p, D):
            plt.plot( [ p[0] ], [ p[1] ], 'ko', ms=15)
                    
    return C


# Graphs

# In[1375]:


def element(p, S):
    for s in S:
        if s==p: # s[0]==p[0] and s[1]==p[1]:
            return True
    return False
        
def vertices(T):
    S = []
    for t in T:
        A, B, C = t
        if not element(A, S):
            S.append(A)
        if not element(B, S):
            S.append(B)
        if not element(C, S):
            S.append(C)
    return sorted(S, key=lambda p: p[0])

def edges(T):
    S = []
    for t in T:
        A, B, C = sides(t)
        if not element(A, S) and not element(rev(A), S):
            S.append(A)
        if not element(B, S) and not element(rev(B), S):
            S.append(B)
        if not element(C, S) and not element(rev(C), S):
            S.append(C)
    return sorted(S, key=lambda p: p[0][0])

def graph(T):
    return vertices(T), edges(T)


# Lists and Paths

# In[1376]:


def pathLength(P):
    return 0 if len(P)<=1 else             dist(P[0], P[1])+pathLength(P[1:])


# In[1377]:


def reverse(P):
    return [ P[-i] for i in range(1,len(P)+1) ]


# In[1378]:


def index(x, L):
    for i in range(len(L)):
        if x==L[i]: 
            return i
    return None


# In[1379]:


def addWithoutDuplicates(L, X):
    for i in range(len(X)):
        if X[i] not in L:
            L.append(X[i])
    return L


# Maps

# In[1380]:


def createMap(P, p=0.2):
    T = triangulation(P)
    V, E = graph(T)
    
    for t in T:
        if slimTriangle(t):
            s = longestSide(t)
            if s in E:
                E.remove(s)
            elif (s[1], s[0]) in E:
                E.remove((s[1], s[0]))
        else:
            for tt in T:
                if not equalTriangles(t, tt):
                    s = commonSide(t, tt)
                    if s is None:
                        continue
                    if random.random()<p:
                        if s in E:
                            E.remove(s)
                            break
                        elif (s[1], s[0]) in E:
                            E.remove((s[1], s[0]))
                            break
                    
    return V, E


# Place Warehouse Somewhere Inside 

# In[1381]:


def placeWarehouse(M):
    V, _ = M
    dmin = math.inf
    pos = (random.randint(resolution/4, resolution*3/4),
            random.randint(resolution/4, resolution*3/4))
    for p in V:
        dp = dist(p, pos)
        if dp<dmin:
            w = p
            dmin = dp
    return w


# Generate Delivery Points

# In[1382]:


def splitEdge(V, E, s):
    A, B = s
    p = random.uniform(0.2,0.8)
    x = int(A[0]+p*(B[0]-A[0]))
    y = int(A[1]+p*(B[1]-A[1]))
    t = (x,y)
    E.remove(s) 
    E.append((A, t))
    E.append((t, B))
    V.append(t)
    return (V, E), t


# In[1383]:


def addTargets(M, n=5):
    V, E = M
    V, E = V.copy(), E.copy()
    T = []
    # we want to ensure that the beginning of the 
    # sequence of points generated randomly stays
    # the same
    mindist = 0.05
    while len(T)<n:
        S = random.sample(E, 1)
        s = S[0]
        A, B = s
        if dist(A,B)>resolution/20: # avoid targets placed narrowly
            (V, E), t = splitEdge(V, E, s)
            T.append(t)
    return (V, E), T


# Plot Map with Delivery Route

# In[1384]:


def plotMap(G, T=[], P=[], w=None,
            style='r-o', lw=1, ms=3, 
            styleT='go', msT=5,
            styleP='b-o', lwP=3, msP=1,
            stylePT='go', msPT=7,
            styleW='ro', msW=9,
            text=None, grid=False, size=6):
    fig = plt.gcf()
    fig.set_size_inches(size, size)
    resolution = 10000
    plt.xlim(0,resolution)
    plt.ylim(0,resolution)
    if not grid:
        plt.axis('off')
    V, E = G
    for e in E:
        p1, p2 = e
        plt.plot( [ p1[0], p2[0] ],
                  [ p1[1], p2[1] ], 
                  style, lw=lw, ms=ms)
    for t in T:
        plt.plot( [ t[0] ], [ t[1] ], 
                  styleT, ms=msT)
    plt.plot( [ p[0] for p in P ],
              [ p[1] for p in P ], 
              styleP, lw=lwP, ms=msP)
    for p in P:
        if p in T:
            plt.plot( [ p[0] ], [ p[1] ], 
                      stylePT, ms=msPT)
    if w is not None:
        plt.plot( [ w[0] ], [ w[1] ], 
                      styleW, ms=msW)
    if grid:
        plt.grid()
        if text is not None:
            plt.title(text)
    else:
        if text is not None:
            plt.text(0.8*resolution, 0, text)

    plt.show()


# In[1385]:


def generateData(seed=None, nodes=50, customers=100, 
                 plot=False, log=False):

    if seed is None:

        print("Usage:  M, W, C = generateData(seed=None, plot=False, log=False)")
        print("")
        print("  seed  the seed value to be used for data generation. ")
        print("        To test the application use seed=0, it will create")
        print("        a small map, with a very few customer locations and")
        print("        a small set of delivery data.")
        print("        For the generating proper simulation data, use the last")
        print("        four digits of your student ID as seed value.")
        print("")
        print("  log   Controls print output during data generation.")
        print("")
        print("  plot  Controls graphical output during data generation.")
        print("")
        print("Returns:")
        print("")
        print("  M = (V, E) is the generated map given as a graph")
        print("    where V is a list of vertices, with each vertice ")
        print("    given as a pair (x, y) of integer coordinates, ")
        print("    and E is a list of edges, with each edge given")
        print("    as a pair (A, B) of vertices, with each vertex again")
        print("    given as a pair (x, y) of integer coordinates")
        print("")
        print("  W ∈ V  is the location of the distribution warehouse")
        print("    given as a pair (x, y) of integer coordinates")
        print("")
        print("  C ⊆ V  is a list of customer locations")
        print("    given as pairs (x, y) of integer coordinates")
        print("    len(C) gives the number of customers generated")
        print("")
        
        seed = 0
    
    if seed==0:          # generate very simple test data 
        nodes = 10       # number of points in map
        customers = 5    # number of  customers
        grid = True
            
    else:
        grid = False
        
    random.seed(seed)
    
    P = generatePoints(nodes)
    
    M = createMap(P)
    W = placeWarehouse(M)
    MT, C = addTargets(M, customers)

    if log:
        print(f"Generated map with {nodes:d} nodes and " 
              f"{customers:d} customer locations")
    if plot:
        label="" if seed==0 else f"seed={seed:4d}"
        plotMap(MT, T=C, w=W, text=label, grid=grid)
    
    return MT, W, C


# Data Generation is reproducible

# In[1386]:


D1 = generateData(1234)
D2 = generateData(1234)
D1 == D2


# Generate Data
# This section demonstrates how you can generate the test data for the problem.
# 
# General Help Message
# If you use generateData() without any parameters you will get a general help message.

# In[1387]:


M, W, C = generateData()


# Analysing Simple Test Data
# This section illustrates the data structure generated.

# In[1388]:


M, W, C = generateData(seed=0, log=True, plot=True)


# The Graph
# You can identify the points in the grid above. The vertices of the graph are:

# In[1389]:


V, E = M
V


# The edges of the graph are:

# In[1390]:


E


# Customer Addresses
# The customer addresses (green dots in the map) are:

# In[1391]:


C


# Real Sample Data¶
# This section shows sample data as you you may get them for your required simulation.

# In[1392]:


_ = generateData(1234, plot=True, log=True)


# In[1393]:


data = generateData(9999, plot=True, log=True)


# Save sample data as pickle file:

# In[1394]:


import pickle
with open('data.pickled', 'wb') as f:
    pickle.dump(data, f)


# In[1395]:


sampleData = generateData(seed=0)


# In[1396]:


import pickle
with open('sampleData.pickled', 'wb') as f:
    pickle.dump(sampleData, f)


# In[1397]:


import scipy.stats as stats

def histplot(data, title="", xlabel="",
             width=None, height=None):
    
    minx = min(data)
    maxx = max(data)
    μ = np.mean(data)
    σ = np.std(data)
    
    fig = plt.figure()
    fig.set_figwidth(width if width is not None else 4)
    fig.set_figheight(height if height is not None else 2.5)
    ax = fig.gca()
        
    bins=min(50, min(len(data),maxx-minx)//5+1)
    hist=plt.hist(data, density=True, bins=bins)
    plt.xlabel(xlabel)
    plt.ylabel('Density')
    plt.title(title)
        
    x = np.linspace(minx, maxx, 100)
    y = [ stats.norm(loc=μ, scale=σ).pdf(p) for p in x]
    ax.plot(x, y, lw=1, color='red')
    ax.axvline(x=μ, color='red')
    maxy = max(max(y), max(hist[0]))
    ax.text(maxx, maxy, 
            f'μ={μ:2.2f}\nσ={σ:2.2f}', 
            ha='right', va='top', 
            color='red', fontsize=12)
    ax.grid(True)
    plt.show()


# In[1398]:


def dailyPlot(data,
              title="", ylabel="",
              width=None, height=None):
    
    days = len(data)
    
    fig = plt.figure()
    fig.set_figwidth(width if width is not None else 6)
    fig.set_figheight(height if height is not None else 2)
    
    ax = fig.gca()
    diff = (max(data)-min(data))*0.1
    ymin = int(math.floor(min(data)-diff))
    ymax = int(math.ceil(max(data)+diff))
    ax.set_xlim(-1, days)
    ax.set_ylim(ymin, ymax)
    ax.grid(True)
    
    ms = 2 if len(data)>100 else 5
    lw = 0.5 if len(data)>100 else 1

    x = np.arange(0, len(data))
    y = np.array([ y for y in data ])
    b, m = np.polynomial.polynomial.polyfit(x, y, 1)
    
    plt.plot(x, y, 'bo-', linewidth=lw, markersize=ms)
    plt.plot(x, m*x+b, 'r-')
    
    plt.xlabel('Day')
    plt.ylabel(ylabel)
    plt.title(title)
    plt.show()


# In[1399]:


def barPlot(count, title=''):
    fig = plt.figure()
    fig.set_figwidth(3 if len(count)<=10 else 5)
    fig.set_figheight(2)
    ax = fig.gca()
    xpos = np.arange(0, len(count))
    labels = [ c for c in range(0, len(count)) ]
    ax.set_xticks(xpos, labels)
    plt.title(title)
    plt.bar(range(0, len(count)), count)


# Geometry

# In[1400]:


def dist(p1, p2):
    (x1, y1) = p1
    (x2, y2) = p2
    return int(math.ceil(math.sqrt((x1-x2)**2+(y1-y2)**2)))


# In[1401]:


def pathLength(P):
    return 0 if len(P)<=1 else             dist(P[0], P[1])+pathLength(P[1:])


# Time Handling
# Convention: In this project we measure time in seconds. The simulation will start at 0:00. Time related methods will be added as they are needed.

# In[1402]:


def day(now):
    return int(now//(24*3600))


# timestamp(t) generates a timestamp string in the form [dd] hh:mm:ss.d. It will be used for trace messages.

# In[1403]:


timestamp(24*3600*3+17*3600+615.1)


# `nextHour` generates the simulation time for the next time the hour is equal to the requested hour. 

# In[1404]:


def nextHour(env, hour):
    beginningOfDay = day(env.now)*24*3600
    timeOfDay = env.now-beginningOfDay
    if hour*3600 > timeOfDay:
        return hour*3600 - timeOfDay
    else:
        return hour*3600 + 24*3600 - timeOfDay


# Generation of Input Data
# 
# Generate Parcels
# 
# generateParcels returns for each day a list of customers that should receive parcels on that day. The data are generated for a given number of customers and days under assumption that the demand is independent and identically distributed with a mean value
# 

# In[1405]:


def generateParcels(days=50, customers=100, 
                    p=0.25, plot=False, log=False):

    # daily delivery demand is generated for each customer 
    # independently using an expovariate distribution
    
    parcelsForCustomersPerDay = [ [] for i in range(days) ]
    parcelsPerCustomer = [ 0 for i in range(customers) ]
    parcelsPerDay = [ 0 for i in range(days) ]

    for c in range(customers):
        cp = []
        t = random.expovariate(p/(24*3600)) # arrival time secs
        d = day(t)
        deliveryDays = []
        while d<days:
            if d not in deliveryDays:
                deliveryDays.append(d)
            parcelsForCustomersPerDay[d].append(c)
            parcelsPerDay[d] += 1
            parcelsPerCustomer[c] += 1
            if d not in cp:
                cp.append(d)
            t += random.expovariate(p/(24*3600))
            d = day(t)
        
        if log:
            print(f"customer {c:3d} "
                  f"expects {parcelsPerCustomer[c]:2d} deliveries "
                  f"on days {str(cp):s}")
    
    if plot:
        
        histplot(parcelsPerCustomer,
                 xlabel=f'Number of Parcels (over {days:3,d} days, p={p:2.2f})',
                 title=f'Number of Parcels per Customer (N={customers:3,d})')
        
        histplot(parcelsPerDay,
                 xlabel=f'Number of Parcels (over {customers:3,d} customers, p={p:2.2f})',
                 title=f'Number of Parcels per Day (N={days:3,d})')
        
    return parcelsForCustomersPerDay


# Testing Parcel Data Generation

# In[1406]:


random.seed(0)
generateParcels(days=4, customers=5, p=0.5, log=True)


# Statistical Analysis of Parcel Data

# In[1407]:


random.seed(0)
_ = generateParcels(days=1000, customers=1000, p=0.4, plot=True)


# Shortest Path Finding
# A* Algorithm

# In[1408]:


def shortestPath(M, A, B):

    def h(p): 
        return pathLength(p)+dist(p[-1],B)
    
    # candidates C are pairs of the path so far and 
    # the heuristic function of that path, 
    # sorted by the heuristic function, as maintained by
    # insert function
    def insert(C, p):
        hp = h(p)
        c = (p, hp)
        for i in range(len(C)):
            if C[i][1]>hp:
                return C[:i]+[c]+C[i:]
        return C+[c]   
        
    V, E = M
    assert(A in V and B in V)    
    C = insert([], [A])

    while len(C)>0:
        # take the first candidate out of the list of candidates
        path, _ = C[0]
        C = C[1:]
        if path[-1]==B:
            return path
        else:
            for (x, y) in E:
                if path[-1]==x and y not in path:
                    C = insert(C, path+[y])
                elif path[-1]==y and x not in path:
                    C = insert(C, path+[x])
    return None


# Floyd-Warshall Algorithm
# The following is a modified Floyd-Warshall Algorithm calculating in parallel the distances and the shortest path between vertices in a graph M.

# In[1409]:


def FW(M):
    
    V, E = M

    n = len(V)
    d = [ [ math.inf for j in range(n) ] for i in range(n) ]
    p = [ [ None for j in range(n) ] for i in range(n) ]

    for (A, B) in E:
        a = V.index(A)
        b = V.index(B)
        d[a][b] = d[b][a] = dist(A, B)
        p[a][b] = [A, B]
        p[b][a] = [B, A]
    
    for i in range(n):
        d[i][i] = 0
        p[i][i] = [V[i]]
    
    for k in range(n):
        for i in range(n):
            for j in range(n):
                dk = d[i][k] + d[k][j]
                if d[i][j] > dk:
                    d[i][j] = dk
                    p[i][j] = p[i][k][:-1] + p[k][j]
                    
    return d, p


# Testing Shortest Path

# In[1410]:


import pickle
with open('data.pickled', 'rb') as f:
    M, W, C = pickle.load(f)


# In[1411]:


V, E = M
d, p = FW(M)
random.seed(0)
[A, B] = random.sample(V, k=2)
a = V.index(A)
b = V.index(B)
P1 = p[a][b]
L1 = d[a][b]
plotMap(M, T=[A, B], P=P1, text=f"FW path: {L1:d}m")
P2 = shortestPath(M, A, B)
L2 = pathLength(P2)
plotMap(M, T=[A, B], P=P2, text=f"A* path: {L2:d}m")
P1==P2


# Create Roundtrip
# Using Iterative Integer Programming

# In[1412]:


import pulp


# In[1413]:


def roundtrips(x, n):
    
    def isElem(x, l):
        for i in range(len(l)):
            if l[i]==x:
                return True
        return False

    def startpoint(trips):
        for i in range(n):
            for t in trips:
                if isElem(i, t):
                    break
            else:
                return i
    
    def totalLength(trips):
        s=0
        for i in range(0, len(trips)):
            s += len(trips[i])-1
        return s

    trips = []
    while totalLength(trips)<n:
        start = startpoint(trips)
        trip = [ start ]
        i = start
        while len(trip) < n-totalLength(trips):
            for j in range(0, n):
                if pulp.value(x[i][j])==1:
                    trip.append(j)
                    i=j
                    break        
            if pulp.value(x[trip[-1]][start])==1:
                trip.append(start)
                break
        trips.append(trip)
    return sorted(trips, key=lambda t: len(t), reverse=True)


# In[1414]:


def createLoop(M, D, P, T):
    V, E = M
    n = len(T)
    # create variables
    x = pulp.LpVariable.dicts("x", ( range(n), range(n) ),
                            lowBound=0, upBound=1, cat=pulp.LpInteger)
    # create problem
    prob = pulp.LpProblem("Loop",pulp.LpMinimize)
    # add objective function
    prob += pulp.lpSum([ D[V.index(T[i])][V.index(T[j])]*x[i][j] 
                             for i in range(n) for j in range(n) ])
    # add constraints
    constraints=0
    for j in range(n):
        prob += pulp.lpSum([ x[i][j] for i in range(n) if i!=j ]) ==1
    constraints += n
    for i in range(n):
        prob += pulp.lpSum([ x[i][j] for j in range(n) if i!=j ]) ==1
    constraints += n
    for i in range(n):
        for j in range(n):
            if i!=j:
                prob += x[i][j]+x[j][i] <= 1
                constraints += 1
    # initialise solver
    solvers = pulp.listSolvers(onlyAvailable=True)
    solver = pulp.getSolver(solvers[0], msg=0)
    # print(f"{constraints:d} Constraints")
    prob.solve(solver)
    trips = roundtrips(x, n)
    while len(trips)>1:
        for t in trips:
            prob += pulp.lpSum([ x[t[i]][t[i+1]] + x[t[i+1]][t[i]]
                            for i in range(0,len(t)-1) ]) <= len(t)-2
            constraints += 1
        # print(f"{constraints:d} Constraints")
        prob.solve(solver)
        trips = roundtrips(x, n)
    trip = trips[0]
    # print(trip)
    loop = []
    for k in range(len(trip)-1):
        sub = P[V.index(T[trip[k]])][V.index(T[trip[k+1]])]
        loop += sub if len(loop)==0 else sub[1:]
    return loop


# Using Heuristic Method

# In[1415]:


def createLoopH(M, D, P, T, plot=False):
    
    def makeLoop(L):
        loop = []
        for i in range(len(L)-1):
            A = L[i]
            B = L[i+1]
            a = V.index(A)
            b = V.index(B)
            sub = P[a][b]
            loop += sub if len(loop)==0 else sub[1:]
        return loop
    
    V, E = M
    W = T[0]
    customers = T[1:]
    if len(T)==1:
        L = T
    elif len(T)<=3:
        L = T + [T[0]]
    else:
        L = T[:3]+[T[0]]
        T = T[3:]
        while len(T)>0:
            if plot:
                loop = makeLoop(L)
                plotMap(M, T=L, P=loop, w=W, 
                        grid=True, text=f"{pathLength(loop):,d}m")
            minExt = math.inf
            minInd = None
            selInd = None
            for k in range(len(T)):
                C = T[k]
                c = V.index(C)
                for i in range(0, len(L)-1):
                    A = L[i]
                    B = L[i+1]
                    a = V.index(A)
                    b = V.index(B)
                    ext = D[a][c] + D[c][b] - D[a][b]
                    if ext<minExt:
                        minExt, minInd, selInd = ext, i+1, k
            L = L[:minInd]+[T[selInd]]+L[minInd:]
            T = T[:selInd]+T[selInd+1:]
    return makeLoop(L)


# Testing Generating Roundtrip

# In[1416]:


import pickle
with open('data.pickled', 'rb') as f:
    M, W, C = pickle.load(f)


# In[1417]:


D, P = FW(M)


# In[1418]:


random.seed(0)
customers = random.sample(C, k=25)


# In[1419]:


L = createLoop(M, D, P, [W]+customers)


# In[1420]:


plotMap(M, T=customers, P=L, w=W, 
        text=f"IP Solution: {pathLength(L):,d}m")


# In[1421]:


L2 = createLoopH(M, D, P, [W]+customers)


# In[1422]:


plotMap(M, T=customers, P=L, w=W, 
        text=f"Heuristic SolutionL {pathLength(L):,d}m")


# Recorder

# In[1423]:


def combineRecorders(recs): 
    M = recs[0].M
    W = recs[0].W
    C = recs[0].C
    days = sum( [ rec.days for rec in recs ] )
    parcels = sum( [ rec.parcels for rec in recs ])
    # ensure the models are compatible
    for rec in recs[1:]:
        assert(M == rec.M)
        assert(W == rec.W)
        assert(C == rec.C)
        
    r = Recorder(None, M, W, C, days, parcels)
    
    # join the data frames
    r.daily = pd.concat([ rec.daily for rec in recs ], ignore_index=True)
    r.precs = pd.concat([ rec.precs for rec in recs ], ignore_index=True)
    return r


# In[1424]:


class Recorder:
    
    def __init__(self, env, M, W, C, days, parcels,
                 log=False, plot=False):
        self.env = env
        self.M = M
        self.W = W
        self.C = C
        self.days = days
        self.parcels = parcels
        self.log = log
        self.plot = plot
        Customer.REGISTER = []
        Parcel.REGISTER = []
                
        # create a data frame for time records per working day
        self.daily = pd.DataFrame()
        self.daily['begin'] = [None]*days
        self.daily['end'] = [None]*days
        self.daily['dist'] = [None]*days
        self.daily['left'] = [None]*days
        
        # create a data frame for records per parcel
        self.precs = pd.DataFrame()
        self.precs['arrived'] = [None]*parcels
        self.precs['delivered'] = [None]*parcels

    def trace(self, event):
        if self.log:
            print(timestamp(self.env.now), event)

    def recordDriverBeginsWork(self):
        self.trace("Driver arrives for work")
        self.daily.at[day(self.env.now), 'begin'] = int(round(self.env.now))
        
    def recordDriverEndsWork(self):
        self.trace("Driver goes home")
        self.daily.at[day(self.env.now), 'end'] = int(round(self.env.now))

    def recordTourLength(self, length):
        self.daily.at[day(self.env.now), 'dist'] = int(length)

    def recordParcelsLeftOver(self, numberOfParcels):
        self.trace(f"{numberOfParcels:d} left over for next day")
        self.daily.at[day(self.env.now), 'left'] = numberOfParcels

    def recordParcelArrival(self, parcel):
        self.trace(f"{str(parcel):s} arrived in delivery centre")
        self.precs.at[parcel.i, 'arrived'] = day(self.env.now)
        
    def recordParcelDelivery(self, parcel):
        self.trace(f"{str(parcel):s} delivered to customer")
        self.precs.at[parcel.i, 'delivered'] = day(self.env.now)

    def finish(self):
        self.daily['time'] = self.daily['end']-self.daily['begin']
        # if at the end of the simulation there are parcels 
        # not yet delivered, register them for the next day
        for i in range(self.parcels):
            if self.precs.at[i, 'delivered'] is None:
                self.precs.at[i, 'delivered'] = self.days
        self.precs['delay'] = self.precs['delivered']-self.precs['arrived']
        # calculate max delay time of parcels delivered on a day
        self.daily['maxdelay'] = [0] * self.days
        for i in range(self.parcels):
            deliveryDay = self.precs.at[i, 'delivered']
            if deliveryDay < self.days:
                if self.daily.at[deliveryDay, 'maxdelay'] < self.precs.at[i, 'delay']:
                    self.daily.at[deliveryDay, 'maxdelay'] = self.precs.at[i, 'delay']
        # simulation is finished for good
        # by removing the simulation environment we can
        # pickle recorder
        self.env = None
        
    def histWorkingTime(self):
        histplot(self.daily['time']//60,
                 xlabel='Working Time [min]',
                 title='Daily Working Time')
        
    def plotWorkingTime(self):
        dailyPlot(self.daily['time']//60,
                  ylabel='Working Time [min]',
                  title='Daily Working Time')
            
    def histTourLength(self):
        histplot(self.daily['dist'],
                 xlabel='Tour Length [m]',
                 title='Daily Tour Length')
            
    def plotTourLength(self):
        dailyPlot(self.daily['dist'],
                  ylabel='Tour Length [m]',
                  title='Daily Tour Length')
            
    def histLeftOver(self):
        histplot(self.daily['left'],
                 xlabel='Left-Over Parcels',
                 title='Daily Left-Over Parcels')
                
    def plotLeftOver(self):
        dailyPlot(self.daily['left'],
                  ylabel='Number of Parcels',
                  title='Daily Left-Over Parcels')
        
    def plotParcelDelay(self):
        dailyPlot(self.daily['maxdelay'],
                  ylabel='Max Delay [days]',
                  title='Max Delivery Delay')
        
    def barplotParcelDelay(self):
        maxDelay = max(self.precs['delay'])
        count = [ 0 for c in range(maxDelay+1) ]
        for i in range(self.parcels):
            count[self.precs.at[i, 'delay']] += 1
        barPlot(count, title='Parcel Delay [days]')
        
    def tableParcelDelay(self):
        maxDelay = max(self.precs['delay'])
        count = [ 0 for c in range(maxDelay+1) ]
        for i in range(self.parcels):
            count[self.precs.at[i, 'delay']] += 1
        print(f'Delivery Delay ({self.parcels:d} parcels)')
        print(f'{"None":>7s}:   {count[0]:4d}   {count[0]/self.parcels*100:4.1f}%')
        for c in range(1, len(count)):
            print(f'{c:2d} days:   {count[c]:4d}   {count[c]/self.parcels*100:4.1f}%')


# Model
# Class Parcel
# 
# Parcels follow through a sequence of states:
# 
# processing
# in transit (from manufacture to distribution centre)
# arrived in distribution centre
# ready for delivery
# out for delivery
# customer not present
# returned to distribution centre
# delivered

# In[1425]:


class Parcel:
    
    REGISTER = []
    
    def __init__(self, rec, cust):
        self.rec = rec
        self.cust = cust
        self.dest = cust.location
        self.status = [ ] # status record and
        self.timing = [ ] # timing
        self.i = len(Parcel.REGISTER)
        Parcel.REGISTER.append(self)

    def __str__(self):
        return f"Parcel: {self.i:2d} for customer {self.cust.i:d}"
    
    def __reg(self, state):
        self.status += [ state ]
        self.timing += [ self.rec.env.now ]
        self.rec.trace(str(self)+" "+state)
        
    def arrivedAtDeliveryCentre(self):
        self.__reg('arr at delivery centre')
        
    def outForDelivery(self): 
        self.__reg('out for delivery')
        
    def returnFromDelivery(self):
        self.__reg('return from delivery')
        
    def delivered(self):
        self.rec.recordParcelDelivery(parcel)
        self.__reg('delivered to customer')


# Class Customer
# The customer is considered passive. The customer is at home with q given probability q. The time to answer the door is given with an expovariate distribution. When the customer is not answering within a certain period of time, the parcel is returned to the delivery centre.

# In[1426]:


class Customer:
    
    REGISTER = []

    def __init__(self, rec, location, q):
        self.rec = rec
        self.location = location
        self.q = q
        self.parcelsReceived = []
        self.i = len(Customer.REGISTER)
        Customer.REGISTER.append(self)
        
        
    def __str__(self):
        return f"Customer: {self.i:2d} {str(self.location):s}"
    
    # factory method ensures that there is only
    # one customer per location
    def getCustomer(rec, location, q=None):
        for c in Customer.REGISTER:
            if c.location == location:
                return c
        assert(q is not None)
        return Customer(rec, location, q)
    
    def responseTime(self):
        if random.random()<self.q:
            return math.inf # customer is not home
        else:
            return random.expovariate(1/AVERAGE_TIME_ANSWER_DOOR)
            
    def acceptParcel(self, parcel):
        self.parcelsReceived.append(parcel.i)
        self.rec.recordParcelDelivery(parcel)


# Class Driver

# In[1427]:


class Driver:
    
    def __init__(self, rec, DC, patience):
        self.rec = rec
        self.DC = DC
        self.patience = patience
        self.location = None
        self.parcels = None
        self.tour = None
        self.rec.env.process(self.process())
        
    # activity
    def __drive(self, target):
        assert(self.tour[0] == self.location)
        while self.location!=target:
            d = dist(self.location, self.tour[1])
            yield self.rec.env.timeout(d / AVERAGE_SPEED)
            self.location = self.tour[1]
            self.tour = self.tour[1:]
        assert(self.tour[0] == self.location == target)
    
    def arriveForWork(self):
        self.location = self.DC.W
        self.parcels = []
        self.returns = []
        self.tour = [ self.DC.W ]
        self.rec.recordDriverBeginsWork()
        
    def leaveForDelivery(self, tour, parcels):
        for p in parcels:
            yield self.rec.env.timeout(PREP_TIME_PER_PARCEL)
            p.outForDelivery()
        self.tour, self.parcels = tour, parcels
        self.rec.trace(f"Driver leaves for delivery "                        f"of {len(parcels):d} parcels")
                
    def goHome(self):
        self.location = self.DC.W
        self.parcels = None
        self.returns = None
        self.tour = [ self.DC.W ]
        self.rec.recordDriverEndsWork()
        
    def process(self):
        yield self.rec.env.timeout(nextHour(self.rec.env, 18))
        while day(self.rec.env.now)<self.rec.days:
            self.arriveForWork()
            tour, parcels = self.DC.sendForDelivery()
            yield from self.leaveForDelivery(tour, parcels)
            self.rec.recordTourLength(pathLength(tour))
            while len(self.parcels)>0:
                # drive to customer
                custLocation = self.parcels[0].dest
                cust = Customer.getCustomer(self.rec, custLocation)
                self.rec.trace("Driver drives to "+str(cust))
                yield from self.__drive(custLocation)
                self.rec.trace("Driver arrived at "+str(cust))
                wait_time = cust.responseTime()
                
                if wait_time<self.patience: 
                    # customer answered door

                    yield self.rec.env.timeout(wait_time)
                    self.rec.trace(str(cust)+" answers door")
                    while len(self.parcels)>0 and                             custLocation == self.parcels[0].dest:
                        handover_time = random.expovariate(1/AVERAGE_TIME_HANDOVER)
                        yield self.rec.env.timeout(handover_time)
                        cust.acceptParcel(self.parcels[0])
                        self.parcels = self.parcels[1:]
                    signoff_time = random.expovariate(1/AVERAGE_TIME_SIGNOFF)
                    yield self.rec.env.timeout(signoff_time)
                    self.rec.trace(str(cust)+" signed off")
                else:
                    # customer not at home or to slow
                    yield self.rec.env.timeout(self.patience)
                    self.rec.trace(str(cust)+" doesn't answer the door")
                    while len(self.parcels)>0 and                              custLocation == self.parcels[0].dest:
                        self.returns.append(self.parcels[0])
                        self.parcels = self.parcels[1:]

            # return to delivery centre
            self.rec.trace("Driver returns to delivery centre")
            yield from self.__drive(self.DC.W)
            yield from self.DC.returnFromDelivery(self.returns)
            
            self.rec.recordParcelsLeftOver(len(self.DC.parcels)+
                                           len(self.DC.leftOver))

            self.goHome()
            yield self.rec.env.timeout(nextHour(self.rec.env, 18))


# Class Delivery Centre

# In[1428]:


class DeliveryCentre:
    
    def __init__(self, rec, M, W, limit, heuristic):
        self.rec = rec
        self.M = M
        self.W = W
        self.D, self.P = FW(M)
        self.limit = limit
        self.heuristic = heuristic
        
        self.leftOver = []    # list of parcels
        self.parcels = []     # list of parcels scheduled for delivery
        self.dest = []        # list of unique customer destinations
        self.tour = [W]       # tour planned for delivery
        
    # builds tour incrementally:
    # every single parcel is in the sequence of their
    # acceptance either incorporated into the tour or
    # assigned to left overs 
    def __accept(self, parcel):
        custLoc = parcel.dest
        if custLoc in self.dest:
            self.parcels.append(parcel)
        else:
            if self.heuristic:
                S = createLoopH(self.M, self.D, self.P, [self.W] + self.dest + [custLoc])
            else:
                S = createLoop(self.M, self.D, self.P, [self.W] + self.dest + [custLoc])
            if pathLength(S)<self.limit:
                self.parcels.append(parcel)
                self.dest.append(custLoc)
                self.tour = S
            else:
                self.leftOver.append(parcel)
            
    def acceptParcel(self, parcel):
        parcel.arrivedAtDeliveryCentre()
        self.rec.recordParcelArrival(parcel)
        self.__accept(parcel)
            
    def sendForDelivery(self):
        parcels = []
        tour = self.tour
        addresses = self.dest
        
        # pick parcels in sequence to be delivered
        for dest in tour:
            if dest in self.dest:
                for p in self.parcels:
                    if p.dest == dest:
                        parcels += [p]
         
        # rearrange the left overs
        L = self.leftOver
        self.tour = [self.W]
        self.parcels = []
        self.leftOver = []
        self.dest = []
        for p in L:
            self.__accept(p)
        
        if self.rec.plot:
            plotMap(self.rec.M, T=addresses, P=tour, w=tour[0], 
                    grid=True, size=5,
                    text=f"Day {day(self.rec.env.now):2d}, {pathLength(tour):,d}m")

        return tour, parcels
                 
    def returnFromDelivery(self, parcels):
        for p in parcels:
            p.returnFromDelivery()
            yield self.rec.env.timeout(RETURN_TIME_PER_PARCEL)
            self.__accept(p)
        yield self.rec.env.timeout(END_OF_DAY_TIME)


# Simulation
# Parameters
# The time required for driving is based on the distance between way points at an average speed of 15km/h.

# In[1429]:


AVERAGE_SPEED = 15/3.6


# The cumulative preparation time (route planning and sorting of the parcels in the delivery order and packing the cargo-bike) is assumed to be 50 sec per parcel to be delivered.

# In[1430]:


PREP_TIME_PER_PARCEL = 50


# The time to process returned parcels in the delivery centre is 30 sec per parcel.

# In[1431]:


RETURN_TIME_PER_PARCEL = 30


# The average customer time to answer the door, accept a parcel, or signoff.

# In[1432]:


AVERAGE_TIME_ANSWER_DOOR = 40
AVERAGE_TIME_HANDOVER = 10
AVERAGE_TIME_SIGNOFF = 10


# The time for end of day closing procedure

# In[1433]:


END_OF_DAY_TIME = 600


# Simulation Routine

# In[1434]:


def simulation(M, W, C,      # map geometry
               days,         # run simulation for number of days
               p=0.2,        # parcels per day and customer, 
               limit=30000,  # bike range limit
               q=0.1,        # probability that the customer is not at home 
               patience=60,  # max wait time for customer answering door
               heuristic=False, log=False, plot=False, ticks=False):

    P = generateParcels(days=days, p=p, customers=len(C))
        
    parcels = sum([ len(d) for d in P ])
    
    print(f"Simulating the delivery of {parcels:d} parcels "
          f"over {days:d} days to {len(C):d} customers") 
    
    env = simpy.Environment()
    rec = Recorder(env, M, W, C, days, parcels, log=log, plot=plot)

    def generatorProcess(env):
                
        DC = DeliveryCentre(rec, M, W, limit, heuristic)
        D = Driver(rec, DC, patience)
        
        # process the parcels day by day
        for CL in P: 
            
            if log:
                print()
            if ticks:
                print(".", end="")

            yield env.timeout(12*3600) # days
            for ci in CL:
                cust = Customer.getCustomer(rec, C[ci], q)
                p = Parcel(rec, cust)
                DC.acceptParcel(p)
            yield env.timeout(12*3600) # days

    env.process(generatorProcess(env))
    env.run()
    
    rec.finish()
    
    return rec


# Model Verification

# In[1435]:


import pickle
with open('sampleData.pickled', 'rb') as f:
    M, W, C = pickle.load(f)


# In[1436]:


C


# In[1437]:


W


# In[1438]:


M


# In[1439]:


plotMap(M, T=C, w=W, text='Sample Data', grid=True)


# In[1440]:


random.seed(0)
_ = simulation(M, W, C, days=4, p=0.4, limit=30000, log=True, plot=True)


# Run Small Simulation

# In[1441]:


import pickle
with open('data.pickled', 'rb') as f:
    M, W, C = pickle.load(f)


# In[1442]:


random.seed(0)
rec = simulation(M, W, C, days=10, p=0.2, limit=30000, 
                 heuristic=True, plot=True)


# Analysing Working Time

# In[1443]:


rec.histWorkingTime()


# In[1444]:


rec.plotWorkingTime()


# Analyse Tour Length

# In[1445]:


rec.histTourLength()


# In[1446]:


rec.plotTourLength()


# Analyse Number of Left-Over Parcels

# In[1447]:


rec.histLeftOver()


# In[1448]:


rec.plotLeftOver()


# Analyse Delayed Parcel Delivery

# In[1449]:


rec.tableParcelDelay()


# In[1450]:


rec.barplotParcelDelay()


# In[1451]:


rec.plotParcelDelay()


# Run Multiple Larger Simulation

# In[1452]:


def multiSimulation(seeds,        # seed values for simulation runs   
                    M, W, C,      # Geometry data
                    days=50,      # run simulation for number of days
                    p=0.2,        # parcels per day and customer, 
                    limit=30000,  # bike range limit
                    q=0.1,        # probability that the customer is not at home 
                    patience=60): # max wait time for customer answering door
    for seed in seeds:
        random.seed(seed)
        filename = f"rec days={days:d} p={p:2.2f} limit={limit:d} " +                    f"q={q:2.2f} patience={patience:d} seed={seed:d}"
        print(filename)
        rec = simulation(M, W, C, 
                         days=days, p=p, limit=limit, q=q, patience=patience, 
                         heuristic=True, ticks=True)
        with open('C:/Users/Sean/Documents'+filename, 'wb') as f:
            pickle.dump(rec, f)
        print()


# In[1453]:


def loadSimulations(seeds,        # seed values for simulation runs   
                    days=50,      # run simulation for number of days
                    p=0.2,        # parcels per day and customer, 
                    limit=30000,  # bike range limit
                    q=0.1,        # probability that the customer is not at home 
                    patience=60): # max wait time for customer answering door
    recs = []
    for seed in seeds:
        random.seed(seed)
        filename = f"rec days={days:d} p={p:2.2f} limit={limit:d} " +                    f"q={q:2.2f} patience={patience:d} seed={seed:d}"
        print(filename)
        with open('C:/Users/Sean/Documents'+filename, 'rb') as f:
            rec = pickle.load(f)
            recs.append(rec)
    return recs


# Run Simulations

# In[1454]:


import pickle
with open('data.pickled', 'rb') as f:
    M, W, C = pickle.load(f)


# In[1455]:


multiSimulation(range(10), M, W, C, limit=45000, days=50)


# Combine Simulation

# In[1456]:


recs = loadSimulations(range(10), limit=45000, days=50)


# In[1457]:


rec = combineRecorders(recs)


# Analysis

# In[1458]:


rec.histWorkingTime()


# In[1459]:


rec.histTourLength()


# In[1460]:


rec.histLeftOver()


# In[1461]:


rec.tableParcelDelay()


# In[1462]:


rec.barplotParcelDelay()


# In[ ]:




