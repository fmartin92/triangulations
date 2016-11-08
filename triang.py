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

def getCoords(vertices, n):
	return cos(2*n*pi/vertices), sin(2*n*pi/vertices)

class Triang:
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


	def __init__(self, treeRep):
		self.treeRep = treeRep
		self.vertices = len(treeRep) + 2
		self.edges = self.findEdges(treeRep, 0, self.vertices-1)

	def draw(self):
		n = self.vertices
		c = canvas.canvas()
		for i in range(n):
			fro0, fro1 = getCoords(n, i)
			to0, to1 = getCoords(n, i+1)
			c.stroke(path.line(fro0, fro1, to0, to1))
		for edge in self.edges:
			fro0, fro1 = getCoords(n, edge[0])
			to0, to1 = getCoords(n, edge[1])
			c.stroke(path.line(fro0, fro1, to0, to1))
		return c

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


if __name__ == '__main__':
	#test: triangulaciones del heptagono
	n = 1
	for i in generateTrees(5):
		Triang(i).draw().writePDFfile("triang" + str(n))
		n += 1