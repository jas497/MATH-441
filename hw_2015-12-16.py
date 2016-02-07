# MATH 441 final project
# James Starkman, jas497
# Paper: Agent-based model (the shorter one)
# Persons with whom the paper was discussed (excluding the professor): no one

import optparse
import time
from random import random as U # allows us to call U() ~ Uniform([0,1])
import math
import matplotlib.pyplot as plt

class Particle:
	def __init__(self, x, y, v=0.03):
		"""This represents a single particle in the system.  It is created with a
velocity of ``v`` in a random direction, at point (x,y).

		"""
		self.xpos = x
		self.ypos = y
		theta = 2 * math.pi * U() # pick a random angle
		self.xvel = v * math.cos(theta)
		self.yvel = v * math.sin(theta)

	def setVelocityAngleMagnitude(self, theta, speed):
		"""Sets ``self.xvel`` and ``self.yvel`` using the given angle and magnitude,
which is counter-clockwise if positive and with units of radians (standard
mathematical convention).

		"""
		self.xvel = speed * math.cos(theta)
		self.yvel = speed * math.sin(theta)

# end class Particle

class Field:
	def __init__(self, N, L, rho, eta):
		"""Creates and manages the field, which is not a lattice, but a continuous square
plane.  Only two parameters are needed; set the other to zero.

N
   Number of particles to create.

L
   Side length of the square field.

rho
   Density, equivalent to N/L^2.

		"""
		# Handle the choose-any-two nature of the parameters
		self.N = (N if N != 0 and N is not None else rho*L*L)
		self.width = (L if L != 0 and L is not None else int(math.sqrt(N/rho)))
		self.height = self.width
		self.eta = eta
		# generate the N agents
		self.particles = [Particle(self.width * U(),
		                           self.height * U()) for _ in range(self.N)]

	def getNeighbors(self, p, r):
		"""Get all other particles within ``r`` of this one.  This model uses periodic
boundary conditions (on a torus), so distance formula is nontrivial.

p
   The main Particle of which to get the neighbors.

r
   Radius around the particle to use.

		"""
		return (other for other in self.particles if self.normDeltaSq(p, other) < r*r)

	def normDeltaSq(self, p1, p2):
		"""Compute the square of the distance between Particles ``p1`` and ``p2``

		"""
		xdist = min(abs(p1.xpos-p2.xpos), self.width - abs(p1.xpos-p2.xpos))
		ydist = min(abs(p1.ypos-p2.ypos), self.height - abs(p1.ypos-p2.ypos))
		return xdist*xdist + ydist*ydist;

	def update(self, delta_t):
		"""Generates a new list of particles, computes their velocities, and then does
an integration step to determine the new particle's x and y position.

		"""
		new_particles = [Particle(0, 0) for _ in range(self.N)]
		# loop over both lists at once:
		for (p, P) in zip(self.particles, new_particles):
			# p is each current particle
			# P is each new particle
			neighbors = self.getNeighbors(p, 1)
			# updating velocity of new particle
			new_angle = math.atan2(sum(n.yvel for n in neighbors),
								   sum(n.xvel for n in neighbors)) + self.eta*(U()-0.5)
			P.setVelocityAngleMagnitude(new_angle, 0.03)
			# numerical integration
			P.xpos = (p.xpos + P.xvel*delta_t) % self.width
			P.ypos = (p.ypos + P.yvel*delta_t) % self.height

		# when finished updating everything, save the new agents and let the GC
		# handle the old
		self.particles = new_particles

# end class Field

def run(N, L, rho, eta):
	dt = 1
	f = Field(N, L, rho, eta)

	plt.ion()           # make interactive
	fig = plt.figure()  # create new figure
	plt.axis([0,L,0,L]) # measuring L by L

	# for each time step
	while True:
		# clear off last frame
		plt.clf()
		# do a time step
		plt.scatter([p.xpos for p in f.particles],
					[p.ypos for p in f.particles])
		f.update(dt)
		plt.pause(0.001)
		print("Drawn at:", time.clock())

def main():
	parser = optparse.OptionParser(usage="""
Usage: %prog [options] eta
  ETA is the width of the interval from which to draw the random perturbation
  Of the following options, two out of three are required (do not use the third).
""")
	parser.add_option("-L", "--side-length",
					  dest="L",
					  type="float",
					  help="specifies the side length of the square field",
					  default=0.0)
	parser.add_option("-N", "--number-particles",
					  dest="N",
					  type="int",
					  help="specifies the number of particles to use",
					  default=0)
	parser.add_option("-r", "--density",
					  dest="rho",
					  type="float",
					  help="specifies the density to use (rho=N/L^2)",
					  default=0.0)
	(options, args) = parser.parse_args()
	
	run(options.N,
		options.L,
		options.rho,
		float(args[0]) if len(args) == 1 else 0.0)

if __name__ == "__main__":
	main()
