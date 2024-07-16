import numpy as np
import scipy.linalg as la
from numpy.linalg import norm
from pdb import set_trace as brk


def find_centre(vertices):
	return sum(vertices)/len(vertices)


def recenter(vertices):
	"""Recentre and normalize the vertices"""
	centre = find_centre(vertices)
	for i, _ in enumerate(vertices):
		vertices[i] -= centre
		vertices[i] /= norm(vertices[i])


class Simplex:
	def __init__(self, dimensions):
		self.vertices = self._find_vertices(dimensions, 1.9)

	def _find_vertices(self, dimensions, edge_len):
		n = dimensions
		assert n >= 1

		# Our 0D simplex can just be a point at the origin
		vertices = [np.zeros((n, 1))]

		# To add another dimension start with the centre of the points you have
		# and extend it "upwards" into the new dimension. Keep the whole thing
		# equilateral as you go by re-centering and re-normalizing
		for i in range(n):
			new_vertex = np.zeros((n, 1))
			new_vertex[i] = 1.0
			vertices.append(new_vertex)
			recenter(vertices)

		return vertices


def main():
	s = Simplex(4)
	for i, v in enumerate(s.vertices):
		j = (i+1) % len(s.vertices)
		print(v.transpose(), norm(v), norm(s.vertices[j]-v))


if __name__ == "__main__":
	main()
