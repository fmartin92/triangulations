from math import pi, cos, sin
from pyx import *
from random import randint

def generateTrees(nodes):
    if nodes == 0:
        return [None]
    else:
        prev = [generateTrees(i) for i in range(nodes)]
        trees = []
        for i in range(nodes):
            for leftBranch in prev[i]:
                for rightBranch in prev[nodes-i-1]:
                    trees.append(Tree(0, leftBranch, rightBranch))
        return trees

def anTree(nodes):
    if nodes == 0:
        return None
    else:
        return Tree(0, anTree(nodes-1), None)

def randomTree(nodes):
    if nodes == 0:
        return None
    left = randint(0, nodes-1)
    right = nodes-1-left
    return Tree(0, randomTree(left), randomTree(right))

def midpoint(x, y):
    return [(x[i] + y[i])/2 for i in range(2)]

def line(p1, p2):
    return path.line(p1[0], p1[1], p2[0], p2[1])

def iterate(f, n, x):
    if n == 1:
        return f(x)
    else:
        return f(iterate(f, n-1, x))

class Tree:
    def __init__(self, root, left=None, right=None):
        self.left = left
        self.right = right
        self.root = root
        self.nodes = 1 + (self.left.nodes if self.left else 0) + \
        (self.right.nodes if self.right else 0)

    def __len__(self):
        return self.nodes

class Triangle:
    def __init__(self, e1, e2,e3):
        self.edges = (e1,e2,e3)
        self.vertices = (e1[0], e1[1], e2[1])

    def __iter__(self):
        return self.edges.__iter__()

    def __contains__(self, item):
        return item in self.edges

    def otherVertex(self, edge):
        for v in self.vertices:
            if v not in edge:
                return v

class Quiver:
    def __init__(self, vertices, arrows):
        self.vertices = vertices
        self.arrows = arrows

    def inNeighbors(self, v):
        return [a[0] for a in self.arrows if a[1]==v]
    def outNeighbors(self, v):
        return [a[1] for a in self.arrows if a[0]==v]
    def neighbors(self, v):
        return self.outNeighbors(v) + self.inNeighbors(v)

class Triangulation:
    def findEdges(self, tree, start, end):
        if tree.left:
            cutoff = 1 + len(tree.left) + start
            if tree.right:
                return [(start, cutoff), (cutoff, end)] + \
                self.findEdges(tree.left, start, cutoff) + \
                self.findEdges(tree.right, cutoff, end)
            else:
                return [(start, cutoff)] + \
                self.findEdges(tree.left, start, cutoff)
        else:
            cutoff = 1 + start
            if tree.right:
                return [(cutoff, end)] + \
                self.findEdges(tree.right, cutoff, end)
            else:
                return []

    def findQuiver(self):
        quiver = []
        n = self.vertices
        for edge in self.edges:
            d = edge[1] - edge[0]
            for i in range(1, d):
                curEdge = (edge[0], edge[1]-i)
                if curEdge in self.edges:
                    quiver.append((edge, curEdge))
                    break
            for i in range(1, n-d):
                if edge[0]-i < 0:
                    curEdge = (edge[1], edge[0]-i+n)
                else:
                    curEdge = (edge[0]-i, edge[1])
                if curEdge in self.edges:
                    quiver.append((edge, curEdge))
                    break
        return Quiver(self.edges, quiver)

    def findTriangles(self):
        triangles = []
        n = self.vertices
        aux = [(i,i+1) for i in range(n-1)] + [(0,n-1)] + self.edges
        for edge1 in aux:
            for edge2 in [edge for edge in aux if edge[0] == edge1[1]]:
                edge3 = (edge1[0], edge2[1])
                if edge3 in aux:
                    triangles.append(Triangle(edge1, edge2, edge3))
        return triangles

    def findInnTriangles(self):
        isInner = lambda edge: edge in self.edges
        return [triangle for triangle in self.triangles if \
        reduce(lambda x, y: x and y, [isInner(edge) for edge in triangle])]

    def inCycle(self, v0, v1, v2):
        for i in self.innTriangles:
            if v0 in i and v1 in i and v2 in i:
                return True
        return False

    def generateProjectives(self):
        projs = ()
        for v in self.quiver.vertices:
            projs += (self.projective(v),)
        return projs

    def __init__(self, treeRep):
        self.treeRep = treeRep
        self.vertices = len(treeRep) + 2
        self.edges = self.findEdges(treeRep, 0, self.vertices-1)
        self.quiver = self.findQuiver()
        self.triangles = self.findTriangles()
        self.innTriangles = self.findInnTriangles()
        self.projectives = self.generateProjectives()

    def vCoords(self, n):
        return cos(2*n*pi/self.vertices), sin(2*n*pi/self.vertices)

    def draw(self):
        c = canvas.canvas()
        #poligono
        for i in range(self.vertices):
            c.stroke(line(self.vCoords(i), self.vCoords(i+1)))
        #aristas internas
        for edge in self.edges:
            c.stroke(line(self.vCoords(edge[0]), self.vCoords(edge[1])))
        #flechas del quiver
        for arrow in self.quiver.arrows:
            fro = midpoint(self.vCoords(arrow[0][0]), \
            self.vCoords(arrow[0][1]))
            to = midpoint(self.vCoords(arrow[1][0]), \
            self.vCoords(arrow[1][1]))
            c.stroke(line(fro, to), [deco.earrow(size=0.1), color.rgb.blue])
        return c

    def stringExtensions(self, string):
        """Finds all possible one-letter extensions for a string,
        including both direct and inverse letters"""
        candidates = self.quiver.neighbors(string[-1])
        candidates.remove(string[-2])
        return filter(lambda x: not self.inCycle(string[-1], \
        string[-2], x), candidates)

    def stringOutExtensions(self, string):
        """Finds all possible one-direct-letter extensions for a string"""
        candidates = self.quiver.outNeighbors(string[-1])
        if string[-2] in candidates:
            candidates.remove(string[-2])
        return filter(lambda x: not self.inCycle(string[-1], \
        string[-2], x), candidates)

    def simple(self, v):
        return Module(self, v, v)

    def projExtend(self, string):
        """Returns the longest direct-letter extension of the string"""
        ext = self.stringOutExtensions(string)
        if ext:
            string.append(ext[0])
            return self.projExtend(string)
        else:
            return string

    def projective(self, v):
        out = self.quiver.outNeighbors(v)
        if len(out) == 0:
            return self.simple(v)
        elif len(out) == 1:
            return Module(self, v, self.projExtend([v, out[0]])[-1])
        else:
            return Module(self, self.projExtend([v, out[1]])[-1], \
            self.projExtend([v, out[0]])[-1])

    def randomModule(self):
        start = self.quiver.vertices[randint(1, len(v)-1)]
        end =  self.quiver.vertices[randint(1, len(v)-1)]
        return Module(self, start, end)
            
class Module:
    def generateString(self, start, end):
        if start == end:
            return (start,)
        paths = [[start, v] for v in self.triang.quiver.neighbors(start)]
        while paths:
            curPath = paths.pop(0)
            if curPath[-1] == end:
                return curPath
            else:
                paths.extend([curPath+[i] for i in \
                self.triang.stringExtensions(curPath)])

    def getTop(self):
        if len(self.string) == 1:
            return (self.string[0],)
        else:
            top = ()
            arrows = self.triang.quiver.arrows
            if (self.string[0], self.string[1]) in arrows:
                top += (self.string[0],)
            for i in range(1, len(self.string)-1):
                if (self.string[i], self.string[i+1]) in arrows \
                and (self.string[i], self.string[i-1]) in arrows:
                    top += (self.string[i],)
            if (self.string[-1], self.string[-2]) in arrows:
                top += (self.string[-1],)
        return top

    def getSocle(self):
        if len(self.string) == 1:
            return (self.string[0],)
        else:
            socle = ()
            arrows = self.triang.quiver.arrows
            if (self.string[1], self.string[0]) in arrows:
                socle += (self.string[0],)
            for i in range(1, len(self.string)-1):
                if (self.string[i+1], self.string[i]) in arrows \
                and (self.string[i-1], self.string[i]) in arrows:
                    socle += (self.string[i],)
            if (self.string[-2], self.string[-1]) in arrows:
                socle += (self.string[-1],)
        return socle

    def __init__(self, triang, start, end):
        self.triang = triang
        self.string = self.generateString(start, end)
        self.top = self.getTop()
        self.socle = self.getSocle()

    def __eq__(self, other):
        return (self.string == other.string) or \
               (self.string[::-1] == other.string)

    def __repr__(self):
        rep = ''
        for i in range(len(self.string)-1):
            rep += str(self.string[i])
            if (self.string[i], self.string[i+1]) \
            in self.triang.quiver.arrows:
                rep += ' -> '
            else:
                rep += ' <- '
        rep += str(self.string[-1])
        return rep

    def draw(self, c=None, col=color.rgb.red):
        if c is None:
            c = self.triang.draw()
        #el caso simple se hace aparte
        if len(self.string) > 1:
            [t1] = [t for t in self.triang.triangles \
            if self.string[0] in t and self.string[1] not in t]
            [t2] = [t for t in self.triang.triangles \
            if self.string[-1] in t and self.string[-2] not in t]
        else:
            [t1, t2] = [t for t in self.triang.triangles \
                if self.string[0] in t]

        fro = self.triang.vCoords(t1.otherVertex(self.string[0]))
        to = self.triang.vCoords(t2.otherVertex(self.string[-1]))
        c.stroke(line(fro, to), [col])
        return c

    def isProjective(self):
        for proj in self.triang.projectives:
            if self == proj:
                return True
        return False

    def auxOmega(self, branch, stable):
        for v in branch:
            if v not in self.string:
                mod = Module(self.triang, v, branch[-1])
                return () if (mod.isProjective() and stable) else (mod,)
        return ()

    def omega(self, stable=False):
        mods = ()
        valleys = [v for v in self.socle \
        if v not in [self.string[0], self.string[-1]]]

        if len(self.top) > 1:
            leftmost = self.top[0]
            rightmost = self.top[-1]
            outLeft = self.triang.quiver.outNeighbors(leftmost)
            outRight = self.triang.quiver.outNeighbors(rightmost)

            if len(outLeft) == 2:
                outLeft.remove(self.string[self.string.index(leftmost)+1])
                leftBranch = self.triang.projExtend([leftmost, outLeft[0]])
                mods += self.auxOmega(leftBranch, stable)

            for valley in valleys:
                simple = self.triang.simple(valley)
                if not (simple.isProjective() and stable):
                    mods += (simple,)

            if len(outRight) == 2:
                outRight.remove(self.string[self.string.index(rightmost)-1])
                rightBranch = self.triang.projExtend([rightmost, outRight[0]])
                mods += self.auxOmega(rightBranch, stable)

        else:
            for v in self.triang.quiver.outNeighbors(self.top[0]):
                branch = self.triang.projExtend([self.top[0], v])
                mods += self.auxOmega(branch, stable)

        return mods

    def projectiveCover(self):
        return tuple(self.triang.projective(v) for v in self.top)

    def nthProjective(self, n):
        if n == 1:
            return self.projectiveCover()
        else:
            mods = ()
            for syz in self.omega(False):
                mods += syz.nthProjective(n-1)
            return mods

    def projResolution(self, n):
        return tuple(self.nthProjective(i) for i in range(1, n))

    def Hom(self, target):
        for v in self.top:
            if v not in target.string:
                return 0
        return 1

    def Ext1(self, target):
        omega = self.omega()
        dim = 0
        for m in omega:
            dim += m.Hom(target)
        return dim


def dr(obj):
    obj.draw().writePDFfile("random")

def ddr(mods):
    mods[1].draw(mods[0].draw(col=color.rgb.green)).writePDFfile("random")

if __name__ == '__main__':
    t = Triangulation(randomTree(15))
    t.draw().writePDFfile("random")
    v = t.quiver.vertices
    m1 = t.randomModule()
    m2 = t.randomModule()
    mm = [m1, m2]
    print m1
    print m2
    ddr(mm)
    #for n, mod in enumerate(t.quiver.vertices):
    #   t.projective(mod).draw().writePDFfile("proj" + str(n))