from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import Rod1_3D
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
		rci = (O[i] + O[i-1]) / 2.0
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

	for i in range(1,n+1):
		mi = obj_list[i-1].getMass()
		rci = (O[i] + O[i-1]) / 2.0
		for k in range(1,i+1):
			Zkm1 = obj_list[k-1].z_axis
			phi[k-1] -= mi * (np.transpose(g) @np.cross( Zkm1 , (rci - O[k-1]) ))
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
			Zjm1 = obj_list[j].z_axis
			Jvi[:,j] = np.cross(Zjm1, (Oci - O[j]))
		Dq += mi* (np.transpose(Jvi) @ Jvi)

		# PE part -----------------
		Ii = obj_list[i].Ic
		Ri = obj_list[i].get_rotation_matrix()[0:3,0:3]	
		Jwi = np.zeros((3,n))
		for j in range(i+1):
			Jwi[:,j] = obj_list[j].z_axis
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
		Zkm1 = obj_list[k].z_axis
		for obj_index in range(n):
			# this loops for different object Jvi that contribute to the Dq nXn matrix

			mi = obj_list[obj_index].getMass()
			Jvi = np.zeros((3,n))
			Jvidqk = np.zeros((3,n))
			Oci = (O[obj_index+1] + O[obj_index]) / 2 #assuming uniform centre of mass
			for j in range(obj_index+1):
				Zjm1 = obj_list[j].z_axis
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
				Jwi[:,j] = obj_list[j].z_axis
			
			# when k==i then only Ri will contain term related to qk
			if k==obj_index:
				C[:,:,k] += (np.transpose(Jwi) @ Rdi @ Ii @ np.transpose(Ri) @ Jwi) + (np.transpose(Jwi) @ Ri @ Ii @ np.transpose(Rdi) @ Jwi)
			# Dq += np.transpose(Jwi) @ Ri @ Ii @ np.transpose(Ri) @ Jwi

	return C

def createCArray(Dq_derivative, dqdt):
	n = dqdt.shape[0]
	C = np.zeros((n,1))
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

	# print(itr, d)
	return d2qkp1


g = np.array([0,-10,0])
w = 0
obj_0_dq = np.array([w])
obj_1_dq = np.array([w])
obj_2_dq = np.array([w])

# x_alpha=0, z_length=0 means both the joint lies on the same plane and have
# same axis of rotation

# obj is made as per the theta = 45'
obj_0 = Rod1_3D.Rod1_3D('obj_0', 4,4,0, 6,4,0, x_alpha=0,z_length=0, color="red")
obj_0.showText = False
obj_0.setDq(obj_0_dq)

obj_1 = Rod1_3D.Rod1_3D('obj_1', 6,4,0, 6,4,2, x_alpha=np.pi/4,z_length=0, color="blue")
obj_1.showText = False
obj_1.setDq(obj_1_dq)
obj_0.addChild(obj_1)

# obj_2 = Rod1_3D.Rod1_3D('obj_2', 8,4,0, 10,4,0, x_alpha=0,z_length=0, color="green")
# obj_2.showText = False
# obj_2.setDq(obj_2_dq)
# obj_1.addChild(obj_2)

# obj_list = [obj_0]
obj_list = [obj_0, obj_1]
# obj_list = [obj_0, obj_1, obj_2]
n = len(obj_list)

# torque = np.array([0], dtype=np.float32)
torque = np.array([0, 0], dtype=np.float32)
# torque = np.array([0,0, 0], dtype=np.float32)

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

	# dqdt = np.array([0,1,1], dtype=np.float32)
	# dqdt = np.array([-1,0], dtype=np.float32)
	# dqdt = np.array([np.sqrt(20 - 5)], dtype=np.float32)

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
		# d2qdt2 = fixed_point_method(torque.copy(), dt, obj_list.copy())
		# more of a damping procedure-----------------------

		dqdt += d2qdt2*dt
		q += dqdt*dt
		t += dt

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

			obj.update()
			if len(obj.child ) == 0:
				continue
			else:
				for i in range(len(obj.child)):
					links.append(obj.child[i])
		# print("t: ", np.array([t]), " q : ", q * 180 / np.pi % 360, "dqdt : " , dqdt, "d2qdt2 : ", d2qdt2, "rhs : ", rhs, "ke : ", np.array([ke]), "pe : ", np.array([pe]))
		# print("t: ", np.array([t]), "d2qdt2 : ", d2qdt2)
		print("t: ", np.array([t]), "energy : ", np.array([ke+pe]))

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


		# draw axis ->
		axis_x, axis_y, axis_z = np.zeros((3,3))
		axis_length = 10
		axis_u, axis_v, axis_w = np.array([[axis_length*1,0,0],[0,axis_length*1,0],[0,0,axis_length*1]])
		ax.quiver(axis_x,axis_y,axis_z,axis_u,axis_v,axis_w,arrow_length_ratio=0.1)


		fig.canvas.draw()
		fig.canvas.flush_events()




# Notes --------------------------------

# this is for revolute joint change at this point for prismatic joint
# there is some problem with this one
# since this value is dependent on theta which is quadratically integerated
# thus its value depend on how small dt is
# thus not giving an exact answer since area intergerated is not correct.
# since get transformation matrix depends on the value of current theta
