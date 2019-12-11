from mpl_toolkits.mplot3d import Axes3D
from matplotlib import patches
from matplotlib.collections import PatchCollection
import numpy as np
import RevoluteJoint
import PrismaticJoint
import World
import matplotlib.pyplot as plt
import pdb
import lemkelcp as lcp


def interpolate_frames(base_frames, dt, method='linear'):
	new_frames = []
	num_frames = len(base_frames)
	num_int_frames = np.int(1/dt)
	for i in range(num_frames-1):
		int_frames = []
		pose_1 = base_frames[i]
		pose_2 = base_frames[i+1]
		for j in range(num_int_frames):
			int_frame = ( pose_2*j + pose_1*(num_int_frames-j)) / num_int_frames
			int_frames.append(int_frame)
		new_frames = new_frames + int_frames
	return new_frames

def interpolate_velocities(frames, dt):
	n = len(frames)
	dof = frames[0].shape[0]
	vels = [np.zeros((dof,))]
	for i in range(n-1):
		vel = (frames[i+1] - frames[i])/dt
		vels.append(vel)
	vels[0] = vels[1]
	return vels

g = np.array([0,0,-10])
# g = np.array([0,0,0])
# same axis of rotation
origin=np.array([6,4,10,4], dtype=np.float32)
# Always ensure that stating axis is always inclined with cartesian axis

showText = False

pobj_0 = PrismaticJoint.PrismaticJoint('pobj_0', 6,4,10, 6,4,10,q_angle=np.pi/2, x_alpha=np.pi/2,r_length=0, color="red")
pobj_0.setMass(0)
pobj_0.showText = showText
pobj_1 = PrismaticJoint.PrismaticJoint('pobj_1', 6,4,10, 6,4,10,q_angle=np.pi/2, x_alpha=np.pi/2,r_length=0, color="blue")
pobj_1.setMass(0)
pobj_0.addChild(pobj_1)

l1obj_1 = RevoluteJoint.RevoluteJoint('l1obj_1', 6,4,10, 10,4,10, x_length=4, x_alpha=0,z_length=2, color="pink")
l1obj_1.showText = showText
pobj_1.addChild(l1obj_1)

l1obj_2 = RevoluteJoint.RevoluteJoint('l1obj_2', 10,4,10, 14,4,10, x_length=4, x_alpha=0,z_length=0, color="cyan")
l1obj_2.showText = showText
l1obj_1.addChild(l1obj_2)

l1obj_3 = RevoluteJoint.RevoluteJoint('l1obj_3', 14,4,10, 15,4,10, x_length=1, x_alpha=0,z_length=0, color="orange")
l1obj_3.showText = showText
l1obj_2.addChild(l1obj_3)

r1obj_1 = RevoluteJoint.RevoluteJoint('r1obj_1', 6,4,10, 10,4,10, x_length=4, x_alpha=0,z_length=2, color="pink")
r1obj_1.showText = showText
pobj_1.addChild(r1obj_1)

r1obj_2 = RevoluteJoint.RevoluteJoint('r1obj_2', 10,4,10, 14,4,10, x_length=4, x_alpha=0,z_length=0, color="cyan")
r1obj_2.showText = showText
r1obj_1.addChild(r1obj_2)

r1obj_3 = RevoluteJoint.RevoluteJoint('r1obj_3', 14,4,10, 15,4,10, x_length=1, x_alpha=0,z_length=0, color="orange")
r1obj_3.showText = showText
r1obj_2.addChild(r1obj_3)


l2obj_1 = RevoluteJoint.RevoluteJoint('l2obj_1', 6,4,10, 10,4,10, x_length=4, x_alpha=0,z_length=-2, color="pink")
l2obj_1.showText = showText
pobj_1.addChild(l2obj_1)

l2obj_2 = RevoluteJoint.RevoluteJoint('l2obj_2', 10,4,10, 14,4,10, x_length=4, x_alpha=0,z_length=0, color="cyan")
l2obj_2.showText = showText
l2obj_1.addChild(l2obj_2)

l2obj_3 = RevoluteJoint.RevoluteJoint('l2obj_3', 14,4,10, 15,4,10, x_length=1, x_alpha=0,z_length=0, color="orange")
l2obj_3.showText = showText
l2obj_2.addChild(l2obj_3)

r2obj_1 = RevoluteJoint.RevoluteJoint('r2obj_1', 6,4,10, 10,4,10, x_length=4, x_alpha=0,z_length=-2, color="pink")
r2obj_1.showText = showText
pobj_1.addChild(r2obj_1)

r2obj_2 = RevoluteJoint.RevoluteJoint('r2obj_2', 10,4,10, 14,4,10, x_length=4, x_alpha=0,z_length=0, color="cyan")
r2obj_2.showText = showText
r2obj_1.addChild(r2obj_2)

r2obj_3 = RevoluteJoint.RevoluteJoint('r2obj_3', 14,4,10, 15,4,10, x_length=1, x_alpha=0,z_length=0, color="orange")
r2obj_3.showText = showText
r2obj_2.addChild(r2obj_3)


# ground is the x-y plane
# obj_list = [pobj_0,pobj_1, lobj_1, lobj_2, lobj_3]
# obj_list = [pobj_0,pobj_1, lobj_1, lobj_2, robj_1, robj_2]
# obj_list = [pobj_0,pobj_1, l1obj_1, l1obj_2, l1obj_3, r1obj_1, r1obj_2, r1obj_3]
obj_list = [pobj_0,pobj_1, l1obj_1, l1obj_2, l1obj_3, r1obj_1, r1obj_2, r1obj_3,
			l2obj_1, l2obj_2, l2obj_3, r2obj_1, r2obj_2, r2obj_3]
n = len(obj_list)


# torque = np.array([0], dtype=np.float32)
# torque = np.array([0, 0], dtype=np.float32)
# torque = np.array([1, 1], dtype=np.float32)
torque = np.zeros((n,))
# torque[-1] = 1
# torque = np.array([0,0, 0, 0,0], dtype=np.float32)

if __name__ == '__main__':
	dt = 0.01
	world = World.World(obj_list)

	# some different pose to work with
	# world.setQ(np.array([0,0,np.pi-np.pi/8,0, np.pi/8-np.pi/2, np.pi+np.pi/8,0, -np.pi/8+np.pi/2, 
		# np.pi-np.pi/8,0, np.pi/8-np.pi/2, np.pi+np.pi/8,0, -np.pi/8+np.pi/2], dtype=np.float32))

	# only valid standing pose for dog with point contacts
	world.setQ(np.array([0,0,np.pi,0, -np.pi/2, np.pi,0, np.pi/2, 
		np.pi,0, -np.pi/2, np.pi,0, np.pi/2], dtype=np.float32))
	# world.setQ(np.array([0,0,np.pi,0, 0, np.pi,0, 0,np.pi,0, 0, np.pi,0, 0],dtype=np.float32))

	# world.setQ(np.array([0,0,np.pi-np.pi/8,np.pi/8,0, np.pi,0,0], dtype=np.float32))

	world.setdqdt(np.array([0,0,1,-1,1,-1,1,-1, 1,-1,1,-1,1,-1], dtype=np.float32))

	world.set_link_origin(origin)
	world.set_gravity(g)
	world.set_friction(2)


	on_ground = False
	while True:

		if on_ground:
			on_ground = world.advect(torque, dt/4)
		else:
			on_ground = world.advect(torque, dt)
		world.update()
		world.render()

		# if world.t>0.1:
		# 	break


