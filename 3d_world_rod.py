from mpl_toolkits.mplot3d import Axes3D
from matplotlib import patches
from matplotlib.collections import PatchCollection
import numpy as np
import RevoluteJoint
import PrismaticJoint
import World
import matplotlib.pyplot as plt
import pdb
from special_matrices import *
from algorithms import *
import lemkelcp as lcp
# set options for now I am not considering the
# world else this will be handled by the world

# t_list=[]
# # e_list=[]
# q1_list=[]
# q2_list=[]

g = np.array([0,0,-10])
# g = np.array([0,-10,0])
# g = np.array([-10,0,0])
# g = np.array([-1,-1,-1])
# g = np.array([0,0,0])
# same axis of rotation
origin=np.array([4,4,3,4], dtype=np.float32)
# Always ensure that stating axis is always inclined with cartesian axis

# pobj_0.addChild(pobj_1)


pobj_0 = PrismaticJoint.PrismaticJoint('pobj_0', 4,4,3, 4,4,3, 	q_angle=np.pi/2, x_alpha=np.pi/2,r_length=0, color="red")
pobj_0.setMass(0)
pobj_0.showText = True
pobj_1 = PrismaticJoint.PrismaticJoint('pobj_1', 4,4,3, 4,4,3, 	q_angle=np.pi/2, x_alpha=np.pi/2,r_length=0, color="blue")
pobj_1.setMass(0)
# pobj_2 = PrismaticJoint.PrismaticJoint('pobj_2', 4,4,3, 4,4,3, 	q_angle=np.pi/2, x_alpha=np.pi/2,r_length=0, color="green")
# pobj_2.setMass(0)
pobj_0.addChild(pobj_1)
# pobj_1.addChild(pobj_2)

obj_1 = RevoluteJoint.RevoluteJoint('obj_1', 4,4,3, 6,4,3, x_length=2, x_alpha=0,z_length=0, color="pink")
# obj_1 = RevoluteJoint.RevoluteJoint('obj_1', 4,4,3, 6,4,1, x_length=2, x_alpha=np.pi/2,z_length=0, color="violet")
obj_1.showText = True
pobj_1.addChild(obj_1)

obj_2 = RevoluteJoint.RevoluteJoint('obj_2', 6,4,3, 8,4,3, x_length=2, x_alpha=0,z_length=0, color="cyan")
obj_2.showText = True
obj_1.addChild(obj_2)


# ground is the x-y plane
# obj_list = [pobj_0]
# obj_list = [pobj_0,pobj_1, obj_1, obj_2]
obj_list = [pobj_0,pobj_1, obj_1]
# obj_list = [pobj_0,pobj_1,pobj_2, obj_1, obj_2]
# obj_list = [obj_1, obj_2]
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
	dt = 0.01/2
	world = World.World(obj_list)
	world.set_link_origin(origin)
	world.set_gravity(g)

	# t_list=[]
	# e_list=[]
	# q_list=[]

	collision = True
	# collision = False
	while True:
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

						# vel = Jvi @ dqdt
						# final_vel = -coef_rest * vel
						# init_vel = vel
						# final_acc = (final_vel - init_vel)
						# pdb.set_trace()

						# calculate Impulse
						# end_effector_mass_inv = ( Jvi @ np.linalg.inv(Dq) @ np.transpose(Jvi) )
						# taustar = Dq @ dqdt + dt*(rhs)
						# A = dt* end_effector_mass_inv
						# qlcp = Jvi @ np.linalg.inv(Dq) @ taustar
						# fn = lcp.lemkelcp(A, qlcp)
						# pdb.set_trace()