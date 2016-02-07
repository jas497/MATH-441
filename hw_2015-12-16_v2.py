# MATH 441 final project
# James Starkman, jas497
# Paper: Agent-based model (the shorter one)
# Persons with whom the paper was discussed (excluding the professor): no one

import optparse
import time
from random import random as U # allows us to call U() ~ Uniform([0,1])
import math
import matplotlib.pyplot as plt
import collections # want the deque, a double-ended queue, for a circular buffer
import statistics # want the standard deviation

class Particle:
	def __init__(self, x, y, v=0.03):
		"""This represents a single particle in the system.  It is created with a
velocity of ``v`` (defaults to 0.03 if unspecified) in a random direction,
at point (x,y).

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
		"""Creates and manages the field, which is not a lattice, but a continuous
square plane.  Only two parameters are needed; set the other to zero.  Assumes a
radius around the particle of 1.

N
   Number of particles to create.

L
   Side length of the square field.  Must be an integer.

rho
   Density, approximately equivalent to N/L^2.

		"""
		# Handle the choose-any-two nature of the parameters
		self.N = (N if (N != 0 and N is not None) else int(rho*L*L))
		self.width = (L if (L != 0 and L is not None) else int(math.sqrt(N/rho)))
		self.height = self.width
		self.eta = eta
		# generate the N agents, stored in a triple-nested list.  Use this to
		# iterate over all particles in a given cell:
		# do_stuff(p) for p in self.particles[y][x]
		self.particles = [[[] for __ in range(self.width)] for _ in range(self.height)]
		x = y = 0
		for i in range(self.N):
			x = self.width * U()
			y = self.height * U()
			self.particles[int(y)][int(x)].append(Particle(x, y))

	def normDeltaSq(self, p1, p2):
		"""Computes the square of the distance between Particles ``p1`` and ``p2``, on
the torus that is the field.

		"""
		dx = abs(p1.xpos-p2.xpos)
		xdist = min(dx, self.width - dx)
		dy = abs(p1.ypos-p2.ypos)
		ydist = min(dy, self.height - dy)
		return xdist*xdist + ydist*ydist

	def getLocalAverageAngle(self, p, x, y):
		"""Computes the average angle of the surrounding Particles within one distance
unit of here.  See paper for full explanation.

		"""
		sum_xvel = 0; sum_yvel = 0;
		local_y = [(y-1)%self.height, y, (y+1)%self.height]
		local_x = [(x-1)%self.width , x, (x+1)%self.width ]
		for Y in local_y:
			for X in local_x:
				for P in self.particles[Y][X]:
					# P is each particle on each of the nine tiles
					if self.normDeltaSq(p, P) < 1: # if within radius of influence
						sum_xvel += P.xvel
						sum_yvel += P.yvel

		return math.atan2(sum_yvel, sum_xvel)

	def update(self):
		"""Generates a new list of particles, computes their velocities, and then does
an integration step (delta t = 1) to determine the new particle's x and y position.

		"""
		new_particles = [[[] for __ in range(self.width)] for _ in range(self.height)]
		x = 0; y = 0;
		while y < self.height:
			x = 0
			while x < self.width:
				for p in self.particles[y][x]: # p is each current particle
					# updating velocity of new particle
					new_angle = self.getLocalAverageAngle(p, x, y) + self.eta*(U()-0.5)
					P = Particle(0,0)
					P.setVelocityAngleMagnitude(new_angle, 0.03)
					# numerical integration with dt=1, and we are on a torus
					P.xpos = (p.xpos + P.xvel) % self.width
					P.ypos = (p.ypos + P.yvel) % self.height
					# store in correct tile
					new_particles[int(P.ypos)][int(P.xpos)].append(P)
				x += 1
			y += 1
		# when finished updating everything, save the new agents and let the GC
		# handle the old
		self.particles = new_particles

# end class Field

def tick(f, draw_plots=False):
	# compute statistics and prepare for plotting
	xs = []; ys = [];
	sum_xvel = 0; sum_yvel = 0;
	for horizontal_strip in f.particles:
		for tile in horizontal_strip:
			for p in tile:
				xs.append(p.xpos)
				ys.append(p.ypos)
				sum_xvel += p.xvel
				sum_yvel += p.yvel

	# compute next step
	f.update()

	if draw_plots:
		plt.clf()                        # clear last frame
		plt.axis([0,f.width,0,f.height]) # create main plot
		plt.scatter(xs, ys)              # fill main plot
		plt.pause(0.001)                 # force it to draw and wait at least this long

	# return v_a
	return math.sqrt(sum_xvel*sum_xvel + sum_yvel*sum_yvel) / (f.N * 0.03)


def run(N, L, rho, eta, draw_plots=False, eternal=False, log_level=None):
	f = Field(N, L, rho, eta)
	if draw_plots:
		plt.ion()           # make interactive
		fig = plt.figure()  # create new figure

	i = 0
	v_a_hist = collections.deque(maxlen=10) # to hold the moving average
	convergence = 0.01            # will be increased later
	while True:                   # use Control-C to escape if set to run forever
		v_a = tick(f, draw_plots) # compute next update
		v_a_hist.append(v_a)      # old value is automatically forced out
		i += 1
		if log_level: # show only if verbose
			print("Tick", i, "was drawn at", time.clock(), "with v_a = " + str(v_a))
		if not eternal and i > 10 and statistics.stdev(v_a_hist) < convergence:
			# end early if finite lifespan and strddev is low enough
			break
		if i % 500 == 0: # every 500 ticks
			# lower our standards, and accept a higher strddev
			convergence += 0.01

	return v_a_hist[-1] # return most recent value for v_a

def frange(x, y, jump):
	"""Generates a linear spread of values (like range, but allows floating-point
arguments)."""
	while x < y:
		yield x
		x += jump

def main():
	parser = optparse.OptionParser(usage="""%prog eta [options]
  Use Control-C to terminate the program.  ETA is the width of the interval from
  which to draw the random perturbation (default zero, i.e., deterministic).  Of
  the following options that are not flags, two out of three are required (if a
  third is used, it will be ignored).  Order of the parameters does not matter,
  but take care not to pass ETA to a flag.""")
	parser.add_option("-p", "--draw-plots",
					  dest="draw_plots",
					  action="store_true",
					  help="draw live scatter plot of the field",
					  default=False)
	parser.add_option("-f", "--run-forever",
					  dest="eternal",
					  action="store_true",
					  help="run simulation forever instead of killing when v_a converges",
					  default=False)
	parser.add_option("-q", "--quiet",
					  dest="log_level",
					  action="store_false",
					  help="log nothing")
	parser.add_option("-v", "--verbose",
					  dest="log_level",
					  action="store_true",
					  help="log every tick and frame")
	parser.add_option("-s", "--summarize",
					  dest="summarize",
					  type="int",
					  help="for replicating figure 2.  Argument is 1 for upper, 2 for lower",
					  default=0)
	parser.add_option("-L", "--side-length",
					  dest="side_length",
					  type="int",
					  help="integer specifying the side length of the square field",
					  default=0.0)
	parser.add_option("-N", "--number-particles",
					  dest="qty_particles",
					  type="int",
					  help="integer specifying the number of particles to use",
					  default=0)
	parser.add_option("-r", "--density",
					  dest="rho",
					  type="float",
					  help="floating-point number specifying the density to use (rho=N/L^2)",
					  default=0.0)
	(options, args) = parser.parse_args()

	if options.summarize == 1:
		xs = []; ys = [];
		for eta in frange(0.0, 5.0, 0.1):
			if options.log_level in {None, True}: # if unspecified or verbose
				print("\nDoing a triplet of runs where eta =", eta)
			# their figure 2a always used rho=4, so we will too
			v_a  = run(options.qty_particles, 0, 4, eta, log_level=options.log_level)
			v_a += run(options.qty_particles, 0, 4, eta, log_level=options.log_level)
			v_a += run(options.qty_particles, 0, 4, eta, log_level=options.log_level)
			xs.append(eta)
			ys.append(v_a / 3)

		plt.plot(xs, ys)
		plt.show()

	elif options.summarize == 2:
		# this figure runs fast enough (for small rho, at least) to plot the output live
		plt.ion()
		fig = plt.figure()
		for rho in frange(0.1, 5.0, 0.1):
			if options.log_level in {None, True}: # if unspecified or verbose
				print("\nDoing a run where rho =", rho)
			v_a  = run(0, options.side_length, rho, float(args[0]), log_level=options.log_level)
			v_a += run(0, options.side_length, rho, float(args[0]), log_level=options.log_level)
			v_a += run(0, options.side_length, rho, float(args[0]), log_level=options.log_level)
			plt.scatter(rho, v_a/3) # add new point (does not erase existing points)
			plt.pause(0.001)

		plt.show()

	else: # no -s option, so proceed as normal
		run(options.qty_particles,
			options.side_length,
			options.rho,
			float(args[0]) if len(args) == 1 else 0.0, # eta
			draw_plots = options.draw_plots,
			eternal    = options.eternal,
			log_level  = options.log_level)


if __name__ == "__main__":
	main()
