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

def flatten(lst):
	return [val for sublst in lst for val in sublst]

class Tree:
	def __init__(self, root, left=None, right=None):
		self.left = left
		self.right = right
		self.root = root

	def __len__(self):
		if self.left is None and self.right is None:
			return 1
		elif self.left is None:
			return 1 + len(self.right)
		elif self.right is None:
			return 1 + len(self.left)
		else:
			return 1 + len(self.left) + len(self.right)

class Triangle:
	def __init__(self, e1, e2,e3):
		self.edges = (e1,e2,e3)
		self.vertices = (e1[0], e1[1], e2[1])

	def __iter__(self):
		return self.edges.__iter__()

	def __contains__(self, item):
		if item in self.edges:
			return True
		else:
			return False

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
					curEdge	= (edge[1], edge[0]-i+n)
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

	def __init__(self, treeRep):
		self.treeRep = treeRep
		self.vertices = len(treeRep) + 2
		self.edges = self.findEdges(treeRep, 0, self.vertices-1)
		self.quiver = self.findQuiver()
		self.triangles = self.findTriangles()
		self.innTriangles = self.findInnTriangles()

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
		"""Returns the last letter of the longest direct-letter
		extension of the string"""
		ext = self.stringOutExtensions(string)
		if ext:
			return self.projExtend(string + ext)
		else:
			return string[-1]

	def projective(self, v):
		out = self.quiver.outNeighbors(v)
		if len(out) == 0:
			return self.simple(v)
		elif len(out) == 1:
			return Module(self, v, self.projExtend([v, out[0]]))
		else:
			return Module(self, self.projExtend([v, out[1]]), \
			self.projExtend([v, out[0]]))
			
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

	def top(self):
		if len(self.string) == 1:
			return (self.string[0],)
		else:
			n = len(self.string)
			flag = False
			top = ()
			for i in range(n-1):
				if (self.string[i], self.string[i+1]) \
				in self.triang.quiver.arrows:
					if not flag:
						top += (self.string[i],)
						flag = True
				else:
					if flag:
						flag = False
			if (self.string[-1], self.string[-2]) \
			in self.triang.quiver.arrows:
				top += (self.string[-1],)
		return top

	def __init__(self, triang, start, end):
		self.triang = triang
		self.string = self.generateString(start, end)
		self.top = self.top()

	def draw(self, c=None):
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
		c.stroke(line(fro, to), [color.rgb.red])
		return c

if __name__ == '__main__':
	t = Triangulation(randomTree(16))
	t.draw().writePDFfile("random")
	for n, mod in enumerate(t.quiver.vertices):
		t.projective(mod).draw().writePDFfile("proj" + str(n))
