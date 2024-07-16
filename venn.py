import numpy as np
import itertools
from numpy.linalg import norm
from pdb import set_trace as brk


def solve_quadratic(a, b, c):
	"""Solve a^2x + bx + c = 0"""
	t = np.sqrt(b*b - 4*a*c)
	return (-b + t) / (2*a), (-b - t) / (2*a)


def find_centre(vertices):
	return sum(vertices)/len(vertices)


def recentre(vertices):
	centre = find_centre(vertices)
	for i, _ in enumerate(vertices):
		vertices[i] -= centre


class Simplex:
	def __init__(self, dimensions, edge_len):
		self.vertices = self._find_vertices(dimensions, edge_len)

	def _find_vertices(self, n, edge_len):
		"""Construct a nice equilateral simplex centred on the origin, with
		n dimensions and n+1 vertices"""
		assert n >= 1

		# Start with a 0D simplex, which is just a point at the origin
		vertices = [np.zeros(n)]

		# To add another dimension start at the centre of the points we have
		# and extend "upwards" into the new dimension, keeping things centered
		# and equilateral as we go.
		for i in range(n):
			# The centre is already at the origin, since we recentre at every
			# step
			new_vertex = np.zeros(n)

			# If this is the first edge, it's easy to make it right the length
			if len(vertices) == 1:
				new_vertex[i] = edge_len
			else:
				# Otherwise find the right length by solving the quadratic
				# equation to make the edge length to the previous vertex
				# correct. Because it's all symmetrical the others should all
				# be the right lengths too.
				new_vertex[i] = 1
				b = -2 * np.dot(new_vertex, vertices[-1])
				c = np.dot(vertices[-1], vertices[-1]) - edge_len**2
				new_vertex[i] = solve_quadratic(1, b, c)[0]

			vertices.append(new_vertex)
			recentre(vertices)

		return vertices

	def print(self):
		print("Simplex with {} vertices".format(len(self.vertices)))
		for i, v in enumerate(self.vertices):
			j = (i+1) % len(self.vertices)
			print(v)

	def find_shared_point(self, indices):
		"""Our "Venn diagram" consists of a unit sphere at each vertex of the
		simplex. Indices identifies a list of vertices and therefore spheres.
		So, e.g. [0,1] means the spheres at vertices 0 and 1. Find a point that
		is inside the specified spheres and not inside any others, and verify
		it"""
		ret = find_centre([self.vertices[i] for i in indices])

		# Which spheres is the point inside?
		inside = []
		for i, v in enumerate(self.vertices):
			if norm(ret - v) < 1.0:
				inside.append(i)

		if inside != sorted(indices):
			brk()

		print("{} is inside and only inside spheres {} ({} vertices)".format(
			ret, inside, len(self.vertices)))
		return ret


def test(dimensions, edge_len=np.sqrt(2)):
	s = Simplex(dimensions, edge_len)
	s.print()
	print()

	for length in range(1, dimensions+1):
		for c in itertools.combinations(range(dimensions+1), length):
			s.find_shared_point(c)
	print()


def main():
	for dimensions in range(2, 100):
		test(dimensions, np.sqrt(2))


if __name__ == "__main__":
	main()
