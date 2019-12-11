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
# set options for now I am not considering the
# world else this will be handled by the world

# t_list=[]
# # e_list=[]
# q1_list=[]
# q2_list=[]

g = np.array([0,0,-10])
# g = np.array([0,0,0])

# same axis of rotation
origin=np.array([4,4,2,4], dtype=np.float32)
# Always ensure that stating axis is always inclined with cartesian axis

# pobj_0.addChild(pobj_1)


pobj_0 = PrismaticJoint.PrismaticJoint('pobj_0', 4,4,1, 4,4,1, 	q_angle=np.pi/2, x_alpha=np.pi/2,r_length=0, color="red")
pobj_0.setMass(0)
pobj_0.showText = True
pobj_1 = PrismaticJoint.PrismaticJoint('pobj_1', 4,4,1, 4,4,1, 	q_angle=np.pi/2, x_alpha=np.pi/2,r_length=0, color="blue")
pobj_1.setMass(0)
# pobj_2 = PrismaticJoint.PrismaticJoint('pobj_2', 4,4,1, 4,4,1, 	q_angle=np.pi/2, x_alpha=np.pi/2,r_length=0, color="green")
# pobj_2.setMass(0)
pobj_0.addChild(pobj_1)
# pobj_1.addChild(pobj_2)

obj_1 = RevoluteJoint.RevoluteJoint('obj_1', 4,4,1, 6,4,1, x_length=2, x_alpha=0,z_length=0, color="pink")
# obj_1 = RevoluteJoint.RevoluteJoint('obj_1', 4,4,1, 6,4,1, x_length=2, x_alpha=np.pi/2,z_length=0, color="violet")
obj_1.showText = True
pobj_1.addChild(obj_1)

obj_2 = RevoluteJoint.RevoluteJoint('obj_2', 6,4,1, 8,4,1, x_length=2, x_alpha=0,z_length=0, color="cyan")
obj_2.showText = True
obj_1.addChild(obj_2)


# ground is the x-y plane
# obj_list = [pobj_0]
obj_list = [pobj_0,pobj_1, obj_1, obj_2]
# obj_list = [pobj_0,pobj_1, obj_1]
# obj_list = [pobj_0,pobj_1,pobj_2, obj_1, obj_2]
# obj_list = [obj_1, obj_2]
n = len(obj_list)
torque = np.zeros((n,))


if __name__ == '__main__':
	dt = 0.01
	world = World.World(obj_list)
	world.set_link_origin(origin)
	world.set_gravity(g)
	world.setQ(np.array([0,0,np.pi/2 - np.pi/8,2*np.pi/8], dtype=np.float32))
	world.friction = 1.4
	# world.setdqdt(np.array([10,0,1,1], dtype=np.float32))

	# t_list=[]
	# e_list=[]
	# q_list=[]

	collision = True
	# collision = False
	while True:
		# world.advect(torque, dt)
		world.advect(torque, dt)
		world.update()

		# debug statemetns
		# t_list.append(t)
		# q1_list.append(q[0])

		world.render()

		# if t>10:
		# 	break


# plt.title('energy graph')
# plt.xlabel('time ---->')
# plt.ylabel('energy ---->')
# # plt.plot(t_list, e_list)
# plt.plot(t_list, q1_list, color='r')
# plt.plot(t_list, q2_list, color='b')
# plt.show()

# Notes --------------------------------

# this is for revolute joint change at this point for prismatic joint
# there is some problem with this one
# since this value is dependent on theta which is quadratically integerated
# thus its value depend on how small dt is
# thus not giving an exact answer since area intergerated is not correct.
# since get transformation matrix depends on the value of current theta
