from mpl_toolkits.mplot3d import Axes3D
from matplotlib import patches
from matplotlib.collections import PatchCollection
import numpy as np
import RevoluteJoint
import PrismaticJoint
import matplotlib.pyplot as plt
import pdb
from special_matrices import *
from algorithms import *

# set options for now I am not considering the
# world else this will be handled by the world

# rendering options
np.set_printoptions(precision=4)
plt.ion()
fig = plt.figure()
ax = plt.gca(projection="3d")
ax.set_xlim(0,10)
ax.set_ylim(0,10)
ax.set_zlim(0,10)
fig.canvas.mpl_connect("close_event", lambda event: exit())

def evaluate_potential_energy(obj_list):
	n = len(obj_list)
	P = 0
	O = np.zeros((n+1,3))
	matrix = np.eye(4)
	for i in range(1,n+1):
		matrix = matrix @ obj_list[i-1].get_transformation_matrix() 
		O[i] = (matrix @ np.append([O[0]],1))[0:3]


	for i in range(n):
		mi = obj_list[i].getMass()
		if obj_list[i].type == 'PrismaticJoint':
			rci = O[i+1]
		elif obj_list[i].type == 'RevoluteJoint':
			rci = (O[i] + O[i+1]) / 2.0
		P-= (np.transpose(g) @ rci ) * mi
	return P

def evaluate_potential_energy_derivative(obj_list):
	n = len(obj_list)
	phi = np.zeros((n,))

	O = np.zeros((n+1,3))
	matrix = np.eye(4)
	for i in range(1,n+1):
		matrix = matrix @ obj_list[i-1].get_transformation_matrix() 
		O[i] = (matrix @ np.append([O[0]],1))[0:3]

	for i in range(n):
		mi = obj_list[i].getMass()
		rci = (O[i] + O[i+1]) / 2.0
		for k in range(i+1):
			if k==0:
				Zkm1 = np.array([0,0,1])
			else:	
				Zkm1 = obj_list[k-1].z_axis
			if obj_list[k].type == 'PrismaticJoint':
				phi[k] -= mi * (np.transpose(g) @ Zkm1)
			elif obj_list[k].type == 'RevoluteJoint':				
				phi[k] -= mi * (np.transpose(g) @ np.cross( Zkm1 , (rci - O[k]) ))
	return phi

def evaluate_Dq(obj_list):
	# this function is coorect as per unit test for 2 rods
	n = len(obj_list)
	Dq = np.zeros((n,n))

	# store Oi values for all the link
	O = np.zeros((n+1,3))
	matrix = np.eye(4)
	for i in range(1,n+1):
		matrix = matrix @ obj_list[i-1].get_transformation_matrix() 
		O[i] = (matrix @ np.append([O[0]],1))[0:3]


	for i in range(n):

		# KE part ---------------
		mi = obj_list[i].getMass()
		Jvi = np.zeros((3,n))
		Oci = (O[i+1] + O[i]) / 2 #assuming uniform centre of mass
		for j in range(i+1):
			if j==0:
				Zjm1 = np.array([0,0,1])
			else:	
				Zjm1 = obj_list[j-1].z_axis
			if obj_list[j].type == 'PrismaticJoint':
				Jvi[:,j] = Zjm1
			elif obj_list[j].type == 'RevoluteJoint':
				Jvi[:,j] = np.cross(Zjm1, (Oci - O[j]))
		Dq += mi* (np.transpose(Jvi) @ Jvi)

		# PE part -----------------
		Ii = obj_list[i].Ic
		Ri = obj_list[i].get_rotation_matrix()[0:3,0:3]	
		Jwi = np.zeros((3,n))
		for j in range(i+1):
			if obj_list[j].type == 'PrismaticJoint':
				Jvi[:,j] = np.array([0,0,0], dtype=np.float32)
			elif obj_list[j].type == 'RevoluteJoint':
				if j==0:
					Jwi[:,j] = np.array([0,0,1])
				else:
					Jwi[:,j] = obj_list[j-1].z_axis
		Dq += np.transpose(Jwi) @ Ri @ Ii @ np.transpose(Ri) @ Jwi
	return Dq	


def evaluate_Dq_derivative(obj_list):
	# This derivative is good for analytic case for n=2 joints
	# I want to calculate the derivative of D[i,j] elem with resepect to Qk 
	n = len(obj_list)
	C = np.zeros((n,n,n))

	n = len(obj_list)

	# store Oi values for all the link
	O = np.zeros((n+1,3))
	matrix = np.eye(4)
	for i in range(1,n+1):
		matrix = matrix @ obj_list[i-1].get_transformation_matrix() 
		O[i] = (matrix @ np.append([O[0]],1))[0:3]

	for k in range(n):
		# here k represents qk element with which we are planning to differentiate
		if k==0:
			Zkm1 = np.array([0,0,1])
		else:	
			Zkm1 = obj_list[k-1].z_axis
		for obj_index in range(n):
			# this loops for different object Jvi that contribute to the Dq nXn matrix

			mi = obj_list[obj_index].getMass()
			Jvi = np.zeros((3,n))
			Jvidqk = np.zeros((3,n))
			Oci = (O[obj_index+1] + O[obj_index]) / 2 #assuming uniform centre of mass
			for j in range(obj_index+1):
				if j==0:
					Zjm1 = np.array([0,0,1])
				else:	
					Zjm1 = obj_list[j-1].z_axis
				if obj_list[j].type == 'PrismaticJoint':
					Jvi[:,j] = Zjm1
					# zero value for Jvidqk
				elif obj_list[j].type == 'RevoluteJoint':
					Jvi[:,j] = np.cross(Zjm1, (Oci - O[j]))
					if k<j:
						Jvidqk[:,j] = np.cross(Zjm1, np.cross(Zkm1, (Oci - O[j])) )
					elif(k<obj_index+1):
						Jvidqk[:,j] = np.cross(Zjm1, np.cross(Zkm1, (Oci - O[k])) )
					else:
						Jvidqk[:,j] = np.zeros((3))

			C[:,:,k] += mi * (np.transpose(Jvi) @ Jvidqk + np.transpose(Jvidqk) @ Jvi)
			# Dq += mi* (np.transpose(Jvi) @ Jvi)

			# PE part -----------------
			Ii = obj_list[obj_index].Ic
			Ri = obj_list[obj_index].get_rotation_matrix()[0:3,0:3]	
			Rdi = obj_list[obj_index].get_rotation_matrix_derivative()[0:3,0:3]	
			Jwi = np.zeros((3,n))
			for j in range(obj_index+1):
				if j==0:
					Jwi[:,j] = np.array([0,0,1])
				else:
					Jwi[:,j] = obj_list[j-1].z_axis
			# when k==i then only Ri will contain term related to qk
			if k==obj_index and obj_list[j].type == 'RevoluteJoint':
				C[:,:,k] += (np.transpose(Jwi) @ Rdi @ Ii @ np.transpose(Ri) @ Jwi) + (np.transpose(Jwi) @ Ri @ Ii @ np.transpose(Rdi) @ Jwi)

	return C

def createCArray(Dq_derivative, dqdt):
	n = dqdt.shape[0]
	C = np.zeros((n,))
	for k in range(n):
		for i in range(n):
			for j in range(n):
				C[k] += 0.5 * (Dq_derivative[k,j,i] + Dq_derivative[k,i,j] - Dq_derivative[i,j,k]) * dqdt[i] * dqdt[j]
	return C


def fixed_point_method(torque, dt, obj_list):
	n = len(obj_list)

	q = np.zeros((n))
	dqdt = np.zeros((n))
	for i in range(n):
		q[i] = obj_list[i].getQ()
		dqdt[i] = obj_list[i].getDq()

	Dq = evaluate_Dq(obj_list)
	phi = evaluate_potential_energy_derivative(obj_list)
	Dq_derivative = evaluate_Dq_derivative(obj_list)
	C = createCArray(Dq_derivative, dqdt)	
	rhs = torque.copy()
	for k in range(n):
		rhs[k] = rhs[k] - phi[k] - C[k]
	d2qkp1 = np.linalg.solve(Dq, rhs) # intial guess

	def fixed_point(d2qkp1): 
		dqkp1 = dqdt + d2qkp1*dt
		qkp1 = q + dt*(dqkp1) 
		for obj in obj_list:
			obj.setQ(qkp1[i])
			obj.setDq(dqkp1[i])

		# now these values are for 
		Dq1 = evaluate_Dq(obj_list)
		phi1 = evaluate_potential_energy_derivative(obj_list)
		Dq_derivative1 = evaluate_Dq_derivative(obj_list)
		C1 = createCArray(Dq_derivative1, dqkp1)
		rhs1 = torque.copy()
		for k in range(n):
			rhs1[k] = rhs1[k] - phi1[k] - C1[k]
		d2qkp1_new = np.linalg.solve(Dq1, rhs1)
		return d2qkp1_new

	d2qkp1_temp = fixed_point(d2qkp1)
	diff = (d2qkp1_temp - d2qkp1)
	d =  np.transpose(diff) @ diff
	eps = 1e-5
	itr=0
	while d>eps:
		d2qkp1_temp = d2qkp1.copy()
		d2qkp1 = fixed_point(d2qkp1)
		diff = (d2qkp1_temp - d2qkp1)
		d=np.transpose(diff) @ diff
		itr+=1

	print(itr, d)
	return d2qkp1

def detectCollision(obj):
	# for now lets assume the ground is at x-y plane at z=0
	x1,x2 = obj.lineX
	y1,y2 = obj.lineY 
	z1,z2 = obj.lineZ
	if z1*z2<0:
		return True
	else:
		return False

g = np.array([0,0,-10])
# g = np.array([0,-10,0])
# g = np.array([-10,0,0])
# same axis of rotation
origin=np.array([4,4,1,1], dtype=np.float32)


# pobj_0.addChild(pobj_1)


pobj_0 = PrismaticJoint.PrismaticJoint('pobj_0', 4,4,1, 4,4,1, 	q_angle=0, x_alpha=np.pi/2,r_length=0, color="red")
# pobj_0 = PrismaticJoint.PrismaticJoint('pobj_0', 4,4,1, 4,4,1, 	q_angle=np.pi/2, x_alpha=np.pi/2,r_length=0, color="red")
# pobj_1 = PrismaticJoint.PrismaticJoint('pobj_0', 4,4,1, 4,4,1, 	q_angle=np.pi/2, x_alpha=np.pi/2,r_length=0, color="blue")
# pobj_2 = PrismaticJoint.PrismaticJoint('pobj_0', 4,4,1, 4,4,1, 	q_angle=np.pi/2, x_alpha=np.pi/2,r_length=0, color="green")

# obj_0 = RevoluteJoint.RevoluteJoint('obj_0', 4,4,1, 4,4,1, x_length=0, x_alpha=np.pi/2,z_length=0, color="red")

obj_1 = RevoluteJoint.RevoluteJoint('obj_1', 4,4,1, 6,4,1, x_length=2, x_alpha=0,z_length=0, color="blue")
pobj_0.addChild(obj_1)

obj_2 = RevoluteJoint.RevoluteJoint('obj_2', 6,4,1, 8,4,1, x_length=2, x_alpha=0,z_length=0, color="green")
obj_1.addChild(obj_2)


# ground is the x-y plane
# obj_list = [pobj_0]
# obj_list = [pobj_0, obj_1]
obj_list = [pobj_0, obj_1, obj_2]
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
# torque = np.array([0,0, 0, 0,0], dtype=np.float32)

# set indices 
for i in range(n):
	obj_list[i].setIndex(i)


if __name__ == '__main__':
	dt = 0.01
	q = np.zeros((n))
	dqdt = np.zeros((n))
	d2qdt2 = np.zeros((n))

	# intiatlising the matrices
	# for linear chain
	# TODO generalize it for general tree like structure
	
	dx = 0.0
	dy = 0.0
	dtheta = 0.0
	t = 0

	for i in range(n):
		q[i] = obj_list[i].getQ()
		dqdt[i] = obj_list[i].getDq()
		d2qdt2[i] = obj_list[i].getD2q()

	# dqdt = np.array([5], dtype=np.float32)
	# dqdt = np.array([5,5], dtype=np.float32)
	# dqdt = np.array([1,1,1,1], dtype=np.float32)

	# t_list=[]
	# e_list=[]
	# q_list=[]
	while True:
		rhs = torque.copy()
		ke = 0
		Dq =  evaluate_Dq(obj_list)
		ke = (1/2) * np.transpose(dqdt) @ Dq @ dqdt
		pe = evaluate_potential_energy(obj_list)
		phi = evaluate_potential_energy_derivative(obj_list)
		# create phi array as derivative of pe with qi's

		Dq_derivative = evaluate_Dq_derivative(obj_list) # its is matrix of shape - nXnXn
		C = createCArray(Dq_derivative, dqdt)

		# create Mc matrix
		for k in range(n):
			rhs[k] = rhs[k] - phi[k] - C[k]

		d2qdt2 = np.linalg.solve(Dq, rhs)

		# detect collisions - 
		for obj in obj_list:
			detectCollision(obj)
			# pdb.set_trace()
			# resolve collisions


		# print(d2qdt2)
		dqdt += d2qdt2*dt
		q += dqdt*dt
		t += dt

		# pdb.set_trace()

		obj = obj_list[0]
		links = [obj]
		# update all the links in top down fashion
		# and update the q values for local 
		while len(links) > 0:
			obj = links[0]
			links.remove(links[0])
			i= obj.getIndex()
			obj.setQ(q[i])
			obj.setDq(dqdt[i])
			obj.setD2q(d2qdt2[i])

			obj.update(origin)
			if len(obj.child ) == 0:
				continue
			else:
				for i in range(len(obj.child)):
					links.append(obj.child[i])
		# pdb.set_trace()
		# print("t: ", np.array([t]), " q : ", q, "dqdt : " , dqdt, "d2qdt2 : ", d2qdt2, "rhs : ", rhs, "ke : ", np.array([ke]), "pe : ", np.array([pe]))
		# print("t: ", np.array([t]), "d2qdt2 : ", d2qdt2)
		print("t: ", np.array([t]), "energy : ", np.array([ke+pe]))
		# print("t: ", np.array([t]), "Dq : ", Dq)
		# print("t: ", np.array([t]), "C : ", C)
		# print("t: ", np.array([t]), "Phi : ", phi)
		# print()
		# t_list.append(t)
		# e_list.append(ke+pe)
		# q_list.append(q[0])

		# render object
		for artist in plt.gca().lines + plt.gca().collections + plt.gca().texts:
			artist.remove()

		for i in range(n):
			obj_list[i].drawObject(ax)


		# debug code ----------------------------------
		# O = np.zeros((n+1,3))
		# matrix = np.eye(4)
		# for i in range(1,n+1):
		# 	matrix = matrix @obj_list[i-1].get_transformation_matrix() 
		# 	O[i] = (matrix @ np.append([O[0]],1))[0:3]
		# 	ax.plot( [O[i-1][0], O[i][0]] ,[O[i-1][1], O[i][1]], [O[i-1][2], O[i][2]], marker = '*', color='red')
		# debug code ends ------------------------------

		# draw ground at xy plane
		xx, yy = np.meshgrid(range(0,11,10), range(0,11,10))
		zz = xx * 0
		ax.plot_surface(xx, yy, zz, alpha=0.2, color='r')



		# draw axis ->
		axis_x, axis_y, axis_z = np.zeros((3,3))
		axis_length = 10
		axis_u, axis_v, axis_w = np.array([[axis_length*1,0,0],[0,axis_length*1,0],[0,0,axis_length*1]])
		ax.quiver(axis_x,axis_y,axis_z,axis_u,axis_v,axis_w,arrow_length_ratio=0.1)


		fig.canvas.draw()
		fig.canvas.flush_events()

		# if t>20:
		# 	break

# plt.title('energy graph')
# plt.xlabel('time ---->')
# plt.ylabel('energy ---->')
# # plt.plot(t_list, e_list)
# plt.plot(t_list, e_list)
# plt.show()

# Notes --------------------------------

# this is for revolute joint change at this point for prismatic joint
# there is some problem with this one
# since this value is dependent on theta which is quadratically integerated
# thus its value depend on how small dt is
# thus not giving an exact answer since area intergerated is not correct.
# since get transformation matrix depends on the value of current theta
