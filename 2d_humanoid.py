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

# g = np.array([0,0,-10])
g = np.array([0,0,0])
# same axis of rotation
origin=np.array([6,4,9,4], dtype=np.float32)
# Always ensure that stating axis is always inclined with cartesian axis

# pobj_0.addChild(pobj_1)


pobj_0 = PrismaticJoint.PrismaticJoint('pobj_0', 6,4,9, 6,4,9, 	q_angle=np.pi/2, x_alpha=np.pi/2,r_length=0, color="red")
pobj_0.setMass(0)
pobj_0.showText = True
pobj_1 = PrismaticJoint.PrismaticJoint('pobj_1', 6,4,9, 6,4,9, 	q_angle=np.pi/2, x_alpha=np.pi/2,r_length=0, color="blue")
pobj_1.setMass(0)
pobj_0.addChild(pobj_1)

lobj_1 = RevoluteJoint.RevoluteJoint('lobj_1', 6,4,9, 10,4,9, x_length=4, x_alpha=0,z_length=0, color="pink")
lobj_1.showText = True
pobj_1.addChild(lobj_1)

lobj_2 = RevoluteJoint.RevoluteJoint('lobj_2', 10,4,9, 14,4,9, x_length=4, x_alpha=0,z_length=0, color="cyan")
lobj_2.showText = True
lobj_1.addChild(lobj_2)

lobj_3 = RevoluteJoint.RevoluteJoint('lobj_3', 14,4,9, 16,4,9, x_length=2, x_alpha=0,z_length=0, color="orange")
lobj_3.showText = True
lobj_2.addChild(lobj_3)

robj_1 = RevoluteJoint.RevoluteJoint('robj_1', 6,4,9, 10,4,9, x_length=4, x_alpha=0,z_length=0, color="pink")
robj_1.showText = True
pobj_1.addChild(robj_1)

robj_2 = RevoluteJoint.RevoluteJoint('robj_2', 10,4,9, 14,4,9, x_length=4, x_alpha=0,z_length=0, color="cyan")
robj_2.showText = True
robj_1.addChild(robj_2)

robj_3 = RevoluteJoint.RevoluteJoint('robj_3', 14,4,9, 16,4,9, x_length=2, x_alpha=0,z_length=0, color="orange")
robj_3.showText = True
robj_2.addChild(robj_3)

# ground is the x-y plane
# obj_list = [pobj_0,pobj_1, lobj_1, lobj_2, lobj_3]
obj_list = [pobj_0,pobj_1, lobj_1, lobj_2, robj_1, robj_2]
# obj_list = [pobj_0,pobj_1, lobj_1, lobj_2, lobj_3, robj_1, robj_2, robj_3]
n = len(obj_list)

# set Axis for all the joints ----
# I have kept these axes as I want my initial axes to behave in this manner
# Remember super important how you intialize your axes
# axis_x, axis_y, axis_z = np.array([[1,0,0],[0,1,0],[0,0,1]])
axis_x, axis_y, axis_z = np.array([[1,0,0],[0,1,0],[0,0,1]])
for i in range(0,n):
	matrix = obj_list[i].get_transformation_matrix()[0:3,0:3]
	axis_x = matrix @ axis_x
	axis_y = matrix @ axis_y
	axis_z = matrix @ axis_z
	obj_list[i].set_xaxis(axis_x)
	obj_list[i].set_yaxis(axis_y)
	obj_list[i].set_zaxis(axis_z)

# torque = np.array([0], dtype=np.float32)
# torque = np.array([0, 0], dtype=np.float32)
# torque = np.array([1, 1], dtype=np.float32)
torque = np.zeros((n,))
# torque[-1] = 1
# torque = np.array([0,0, 0, 0,0], dtype=np.float32)

# set indices 
for i in range(n):
	obj_list[i].setIndex(i)


if __name__ == '__main__':
	dt = 0.01
	world = World.World(obj_list)
	# world.setQ(np.array([0,0,np.pi,0, 0, ], dtype=np.float32))
	world.setQ(np.array([0,0,np.pi-np.pi/8,0, np.pi+np.pi/8,0], dtype=np.float32))
	# world.setQ(np.array([0,0,np.pi-np.pi/8,0, 0, np.pi+np.pi/8,0, 0], dtype=np.float32))
	# world.setQ(np.array([0,0,np.pi/2 +3* np.pi/8 ,np.pi/16, -np.pi/8], dtype=np.float32))
	# world.setQ(np.array([0,0,np.pi/2 + np.pi/8,np.pi/8, -np.pi/8, -np.pi/8], dtype=np.float32))
	# world.setdqdt(np.array([1,5,1,1,1], dtype=np.float32))
	world.setdqdt(np.array([0,0,1,-1,1,-1], dtype=np.float32))
	# world.setdqdt(np.array([0,0,0,0,0,0,0,0], dtype=np.float32))

	world.set_link_origin(origin)
	world.set_gravity(g)

	# t_list=[]
	# e_list=[]
	# q_list=[]

	on_ground = False
	while True:

		# DO some inverse kinetics here.............
		# world.advect(torque, dt)
		# if on_ground:
		# 	torque[-1] = 1
		on_ground = world.advect(torque, dt)
		# print(on_ground)
		world.update()

		# debug statemetns
		# t_list.append(world.t)
		# e_list.append(energy)

		world.render()

		# if world.t>10:
		# 	break


# plt.title('energy graph')
# plt.xlabel('time ---->')
# plt.ylabel('energy ---->')
# # plt.plot(t_list, e_list)
# plt.plot(t_list, e_list, color='r')
# # plt.plot(t_list, q2_list, color='b')
# plt.show()

