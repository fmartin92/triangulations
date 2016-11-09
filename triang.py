from math import pi, cos, sin
from pyx import *

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
		return quiver

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
		for arrow in self.quiver:
			fro = midpoint(self.vCoords(arrow[0][0]), \
			self.vCoords(arrow[0][1]))
			to = midpoint(self.vCoords(arrow[1][0]), \
			self.vCoords(arrow[1][1]))
			c.stroke(line(fro, to), [deco.earrow()])
		return c

class Module:
	def __init__(self, triang, string):
		self.triang = triang
		self.string = string

	def draw(self):
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
	#test: triangulaciones del heptagono
	#n = 1
	#for i in generateTrees(6):
	#	Triangulation(i).draw().writePDFfile("triang" + str(n))
	#	n += 1

	string = [(5,7), (5,9)]
	a = Triangulation(generateTrees(8)[16])
	a.draw().writePDFfile("triang")
	mod = Module(a, string)
	mod.draw().writePDFfile("mod")