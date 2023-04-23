import sys
import numpy as np
import matplotlib.pyplot as plt
import random
from numpy import interp
import keyboard

"""
Tom Rice (~2023) has adapted Philip Mocz's code for a project with 
Sasha Solovyeva. Much of what follows was borrowed from
https://github.com/pmocz/nbody-python .

==========

Create Your Own N-body Simulation (With Python)
Philip Mocz (2020) Princeton Univeristy, @PMocz

Simulate orbits of stars interacting due to gravity
Code calculates pairwise forces according to Newton's Law of Gravity
"""

try:
	plt.rcParams['keymap.save'].remove('s')
	plt.rcParams['keymap.xscale'].remove('k')
	plt.rcParams['keymap.yscale'].remove('l')
	plt.rcParams['keymap.fullscreen'].remove('f')
except ValueError:
	pass


notes = [ 60, 62, 64, 65, 67, 69, 71, 72];
keys = [ 65, 83, 68, 70, 71, 72, 74, 75, 76 ];

pianoDict = {
	"a": 261.63,
	"s": 293.66,
	"d": 329.63,
	"f": 349.23,
	"g": 392.00,
	"h": 440.00,
	"j": 493.88,
	"k": 523.25
}

def getAcc( pos, mass, G, softening ):
	"""
	Calculate the acceleration on each particle due to Newton's Law 
	pos  is an N x 3 matrix of positions
	mass is an N x 1 vector of masses
	G is Newton's Gravitational constant
	softening is the softening length
	a is N x 3 matrix of accelerations
	"""
	# positions r = [x,y,z] for all particles
	x = pos[:,0:1]
	y = pos[:,1:2]
	z = pos[:,2:3]

	# matrix that stores all pairwise particle separations: r_j - r_i
	dx = x.T - x
	dy = y.T - y
	dz = z.T - z

	# matrix that stores 1/r^3 for all particle pairwise particle separations 
	inv_r3 = (dx**2 + dy**2 + dz**2 + softening**2)
	inv_r3[inv_r3>0] = inv_r3[inv_r3>0]**(-1.5)

	ax = G * (dx * inv_r3) @ mass
	ay = G * (dy * inv_r3) @ mass
	az = G * (dz * inv_r3) @ mass
	
	# pack together the acceleration components
	a = np.hstack((ax,ay,az))

	return a
	
def getEnergy( pos, vel, mass, G ):
	"""
	Get kinetic energy (KE) and potential energy (PE) of simulation
	pos is N x 3 matrix of positions
	vel is N x 3 matrix of velocities
	mass is an N x 1 vector of masses
	G is Newton's Gravitational constant
	KE is the kinetic energy of the system
	PE is the potential energy of the system
	"""
	# Kinetic Energy:
	KE = 0.5 * np.sum(np.sum( mass * vel**2 ))


	# Potential Energy:

	# positions r = [x,y,z] for all particles
	x = pos[:,0:1]
	y = pos[:,1:2]
	z = pos[:,2:3]

	# matrix that stores all pairwise particle separations: r_j - r_i
	dx = x.T - x
	dy = y.T - y
	dz = z.T - z

	# matrix that stores 1/r for all particle pairwise particle separations 
	inv_r = np.sqrt(dx**2 + dy**2 + dz**2)
	inv_r[inv_r>0] = 1.0/inv_r[inv_r>0]

	# sum over upper triangle, to count each interaction only once
	PE = G * np.sum(np.sum(np.triu(-(mass*mass.T)*inv_r,1)))
	
	return KE, PE;

# Tom writes a function which takes the current parameters and adds an object.
def add_object():
	pass

global globvar, globkey
globvar = 0
globkey = 0


class StarData:
	def __init__ (self, freq, dur, note, numNotes):
		# the data from the most recently played note(s)
		self.freq = freq
		self.mass = interp(self.freq, [27.5, 4186], [0.1, 10])
		self.age = 0
		self.yPos = interp(dur, [.5, 64], [-2, 2])
		self.vel = interp(numNotes, [1, 3], [0, 1])
		# tempo; determined at a later stage via the list, defaults to 0
		self.xPos = 0;
		self.zPos = 0;

global pianoNotes_list, currentSound
pianoNotes_list = []
currentSound = StarData(random.uniform(27.5, 4186), random.uniform(.5, 64), 7, random.randint(1, 3));

def on_press(event):
	print('press', event.key)
	sys.stdout.flush()
	global globvar, globkey 
	global currentSound
	globvar = 1
	globkey = event.key

	frequency = random.uniform(27.5, 4186)

	for key in pianoDict.keys():
		if key == event.key:
			frequency = pianoDict[key]

	# if event.key == 'x':
	#     visible = xl.get_visible()
	#     xl.set_visible(not visible)
	#     fig.canvas.draw()

	currentSound = StarData(frequency, random.uniform(.5, 64), 7, random.randint(1, 3));
	pianoNotes_list.append(currentSound)

	# taking the interval between the current and the previous note
	try:
		interval = abs(currentSound.freq - pianoNotes_list[-2].freq)
	except IndexError:
		interval = 0

	currentSound.xPos = interp(interval, [0, 4158.5], [-2, 2])
	print("xPos:", currentSound.xPos)

	# calculating the tempo and using it for the z position
	# TO BE COMPLETED

	print("Mass ", currentSound.mass)


def main():
	""" N-body simulation """
	
	# Simulation parameters
	N         = 5    # Number of particles
	t         = 0      # current time of the simulation
	tEnd      = 10.0   # time at which simulation ends
	dt        = 0.01   # timestep
	softening = 0.1    # softening length
	G         = 1.0    # Newton's Gravitational Constant
	plotRealTime = True # switch on for plotting as the simulation goes along
	
	# Generate Initial Conditions
	np.random.seed(17)            # set the random number generator seed
	
# 	mass = 20.0*np.ones((N,1))/N  # total mass of particles is 20
	mass = (np.arange(N).reshape(N,1)+1) * 1/100
	pos  = np.random.randn(N,3)   # randomly selected positions and velocities
	vel  = np.random.randn(N,3)
	
	# Convert to Center-of-Mass frame
	vel -= np.mean(mass * vel,0) / np.mean(mass)
	
	# calculate initial gravitational accelerations
	acc = getAcc( pos, mass, G, softening )
	
	# calculate initial energy of system
	KE, PE  = getEnergy( pos, vel, mass, G )
	
	# number of timesteps
	Nt = int(np.ceil(tEnd/dt))
	
	# # save energies, particle orbits for plotting trails
	# pos_save = np.zeros((N,3,Nt+1))
	# pos_save[:,:,0] = pos
	# KE_save = np.zeros(Nt+1)
	# KE_save[0] = KE
	# PE_save = np.zeros(Nt+1)
	# PE_save[0] = PE
	t_all = np.arange(Nt+1)*dt
	
	# prep figure
	fig = plt.figure(figsize=(4*2,5*2), dpi=80)
	grid = plt.GridSpec(3, 1, wspace=0.0, hspace=0.3)
	ax1 = plt.subplot(grid[0:2,0])
	# ax2 = plt.subplot(grid[2,0])

	fig.canvas.mpl_connect('key_press_event', on_press)	
	
	# Simulation Main Loop
	for i in range(Nt):
		# (1/2) kick
		vel += acc * dt/2.0
		
		# drift
		pos += vel * dt
		
		# update accelerations
		acc = getAcc( pos, mass, G, softening )
		
		# (1/2) kick
		vel += acc * dt/2.0
		
		# update time
		t += dt
		
		# get energy of system
		KE, PE  = getEnergy( pos, vel, mass, G )
		
		# # save energies, positions for plotting trail
		# pos_save[:,:,i+1] = pos
		# KE_save[i+1] = KE
		# PE_save[i+1] = PE
		
		# plot in real time
		if plotRealTime or (i == Nt-1):
			plt.sca(ax1)
			plt.cla()
			# xx = pos_save[:,0,max(i-50,0):i+1]
			# yy = pos_save[:,1,max(i-50,0):i+1]
			# plt.scatter(xx,yy,s=1,color=[.7,.7,1])
			plt.scatter(pos[:,0],pos[:,1],s=mass*10,color='blue')
			ax1.set(xlim=(-2, 2), ylim=(-2, 2))
			ax1.set_aspect('equal', 'box')
			ax1.set_xticks([-2,-1,0,1,2])
			ax1.set_yticks([-2,-1,0,1,2])
			
			# plt.sca(ax2)
			# plt.cla()
			# plt.scatter(t_all,KE_save,color='red',s=1,label='KE' if i == Nt-1 else "")
			# plt.scatter(t_all,PE_save,color='blue',s=1,label='PE' if i == Nt-1 else "")
			# plt.scatter(t_all,KE_save+PE_save,color='black',s=1,label='Etot' if i == Nt-1 else "")
			# ax2.set(xlim=(0, tEnd), ylim=(-300, 300))
			# ax2.set_aspect(0.007)
			
			plt.pause(0.001)

			# # alt version: Tom adds new objects on new keypress
			# def on_press(event):
			#     print('press', event.key)
			#     sys.stdout.flush()
			#     if event.key == 'x':
			#         visible = xl.get_visible()
			#         xl.set_visible(not visible)
			#         fig.canvas.draw()


			# Tom is adding new objects on fixed time
			# print(i)
			global globvar, globkey

			if globvar == 1:
				# global globvar
				globvar = 0


			# if i % 50 == 0:
				# print(" i  % 10  == 0 ")
				print("globvar == 1")
				print("N= ", N+1)

				# add an object
				N += 1
				try:
					mass_new = currentSound.mass
				except:
					mass_new = i/100 + np.mean(mass)
				pos_new = [currentSound.xPos, currentSound.yPos, currentSound.zPos]
				vel_new = [0, currentSound.vel, 0]
				acc_new = [0,0,np.sin(i/1000)]

				mass = np.append(mass, mass_new).reshape(N,1)
				pos = np.append(pos, pos_new).reshape(N,3)
				vel = np.append(vel, vel_new).reshape(N,3)
				acc = np.append(acc, acc_new).reshape(N,3)
   
	
	
	# add labels/legend
	# plt.sca(ax2)
	# plt.xlabel('time')
	# plt.ylabel('energy')
	# ax2.legend(loc='upper right')
	
	# Save figure
	plt.show()
		
	return 0
	


  
if __name__== "__main__":
  main()
