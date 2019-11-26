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
import lemkelcp as lcp

def axisAngleRotationMatrix(angle, axis):
	mat = np.zeros((3,3))
	c = np.cos(angle)
	s = np.sin(angle)
	ux = axis[0]
	uy = axis[1]
	uz = axis[2]

	mat[0,0] = c + (1-c)*ux**2
	mat[0,1] = ux*uy*(1-c) - uz*s
	mat[0,2] = ux*uz*(1-c) + uy*s

	mat[1,0] = uy*ux*(1-c) + uz*s
	mat[1,1] = c + (1-c)*uy**2
	mat[1,2] = uy*uz*(1-c) - ux*s
	
	mat[2,0] = uz*ux*(1-c) - uy*s
	mat[2,1] = uz*uy*(1-c) + ux*s
	mat[2,2] = c + (1-c)*uz**2

	return mat


def evaluate_potential_energy(obj_list, g):
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

def evaluate_potential_energy_derivative(obj_list, g):
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

def evaluateJacobian(obj_list, index, point_affect):
	# index is the end point of obj_list[index] at which we are computing the jacobian

	# store Oi values for all the link
	n=len(obj_list)
	O = np.zeros((n+1,3))
	matrix = np.eye(4)
	for i in range(1,n+1):
		matrix = matrix @ obj_list[i-1].get_transformation_matrix() 
		O[i] = (matrix @ np.append([O[0]],1))[0:3]

	i = index
	mi = obj_list[i].getMass()
	Jvi = np.zeros((3,n))
	if point_affect==1:
		Oci = O[i+1] #it has to be the end effector and not the mid-point
	elif point_affect==0:
		Oci = O[i] #it has to be the end effector and not the mid-point
	else:
		Oci = (O[i] + O[i+1]) / 2 # it has to be the mid-point as every point is touching

	for j in range(i+1):
		if j==0:
			Zjm1 = np.array([0,0,1])
		else:	
			Zjm1 = obj_list[j-1].z_axis
		if obj_list[j].type == 'PrismaticJoint':
			Jvi[:,j] = Zjm1
		elif obj_list[j].type == 'RevoluteJoint':
			Jvi[:,j] = np.cross(Zjm1, (Oci - O[j]))
	return Jvi

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

	Ri = np.eye(3) # it grows with time
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

		# RE part -----------------
		Ii = obj_list[i].Ic
		Ri = Ri @ obj_list[i].get_rotation_matrix()[0:3,0:3]	
		Jwi = np.zeros((3,n))
		for j in range(i+1):
			if obj_list[j].type == 'PrismaticJoint':
				Jwi[:,j] = np.array([0,0,0], dtype=np.float32)
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
		
		# compute gradient of z_axis with respect to the derivative of k
		dz_axis_list = [np.array([0,0,0])]
		z_axis_list = [np.array([0,0,1])]
		matrix = np.eye(3)
		dmatrix = np.eye(3)
		R_list=[]
		dR_list=[]		
		for i in range(0,n):
			if i==k:
				dmatrix = dmatrix @ obj_list[i].get_rotation_matrix_derivative()[0:3,0:3]
			else:
				dmatrix = dmatrix @ obj_list[i].get_rotation_matrix()[0:3,0:3]
			matrix = matrix @ obj_list[i].get_rotation_matrix()[0:3,0:3]
			R_list.append(matrix)
			dR_list.append(dmatrix)
			if i<k:
				dz_axis_list.append(np.array([0,0,0]))
			else:
				dz_axis_list.append(dmatrix @ np.array([0,0,1]))
			z_axis_list.append(matrix @ np.array([0,0,1]))

		for obj_index in range(n):
			# this loops for different object Jvi that contribute to the Dq nXn matrix

			mi = obj_list[obj_index].getMass()
			Jvi = np.zeros((3,n))
			dJvi = np.zeros((3,n))
			Oci = (O[obj_index+1] + O[obj_index]) / 2 #assuming uniform centre of mass
			for j in range(obj_index+1):
				Zjm1 = z_axis_list[j]
				dZjm1 = dz_axis_list[j]
				if obj_list[j].type == 'PrismaticJoint':
					Jvi[:,j] = Zjm1
					if k<j:
						if obj_list[k].type == 'RevoluteJoint':
							dJvi[:,j] = dZjm1
						# otherwise no effect on this axis if is prismatic joint
				elif obj_list[j].type == 'RevoluteJoint':
					Jvi[:,j] = np.cross(Zjm1, (Oci - O[j]))
					if k<j:
						dZjm1 = dZjm1
						if obj_list[k].type == 'RevoluteJoint':
							dJvi[:,j] = np.cross(Zjm1, np.cross(Zkm1, (Oci - O[j])) ) + np.cross(dZjm1, (Oci - O[j]))
						# otherwise no effect on this axis if is prismatic joint
					elif(k<obj_index+1):
						# no derivative of Zjm1 as k>j
						dZjm1 = dZjm1
						if obj_list[k].type == 'RevoluteJoint':
							dJvi[:,j] = np.cross(Zjm1, np.cross(Zkm1, (Oci - O[k])) )
						# otherwise no effect on this axis if is prismatic joint
					else:
						dJvi[:,j] = np.zeros((3))

			C[:,:,k] += mi * (np.transpose(Jvi) @ dJvi + np.transpose(dJvi) @ Jvi)
			# Dq += mi* (np.transpose(Jvi) @ Jvi)

			# RE part -----------------
			Ii = obj_list[obj_index].Ic
			# TODO
			Ri = R_list[obj_index]
			Rdi = dR_list[obj_index]


			Jwi = np.zeros((3,n))
			dJwi = np.zeros((3,n))
			for j in range(obj_index+1):
				if obj_list[j].type == 'PrismaticJoint':
					Jwi[:,j] = np.array([0,0,0], dtype=np.float32)
					dJwi[:,j] = np.array([0,0,0], dtype=np.float32)

				elif obj_list[j].type == 'RevoluteJoint':
					Jwi[:,j] = z_axis_list[j]
					dJwi[:,j] = dz_axis_list[j]

			if k<=obj_index:
				C[:,:,k] += \
				(np.transpose(Jwi) @ Rdi @ Ii @ np.transpose(Ri) @ Jwi) + \
				(np.transpose(Jwi) @ Ri @ Ii @ np.transpose(Rdi) @ Jwi) 
				+ (np.transpose(dJwi) @ Ri @ Ii @ np.transpose(Ri) @ Jwi) + \
				(np.transpose(Jwi) @ Ri @ Ii @ np.transpose(Ri) @ dJwi)
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

def calculate_joint_effector_velocities(obj_list):
	v_list = []
	n = len(obj_list)
	for i in range(n):
		obj = obj_list[i]
		if i==0:
			z_axis = np.array([0,0,1])
		else:
			z_axis = obj_list[i-1].z_axis

		if obj.type == 'PrismaticJoint':
			vel = z_axis * obj.getDq()
			if not i==0:
				vel = v_list[-1] + vel
			v_list.append(vel)
		else:
			w = obj.getDq()
			v_dir = np.cross(z_axis, obj.getRadialVector()) * w
			if i==0:
				vel = v_dir
			else:
				vel = v_list[-1] + v_dir
			v_list.append(vel)
	return v_list

def detectCollision(obj, plane):
	# for now lets assume the ground is at x-y plane at z=0
	x1,x2 = obj.lineX
	y1,y2 = obj.lineY 
	z1,z2 = obj.lineZ

	nx,ny,nz = plane[0]
	x0,y0,z0 = plane[1]

	eps = 1e-1

	if ( abs(nx*(x2-x1) + ny*(y2-y1) + nz*(z2-z1) ) ) <= eps:
		# line very close to being parralel to the plane
		if ( abs(nx*(x0-x1) + ny*(y0-y1) + nz*(z0-z1)) ) <= eps:
			# line lies on the plane therefore inf points of contact. 
			# I will return just the end point (Remember)
			point_of_contact = np.array([(x1+x2)/2,(y1+y2)/2,(z1+z2)/2])
			return True, point_of_contact, plane[0], -1
		# elif ( abs(nx*(x0-x2) + ny*(y0-y2) + nz*(z0-z2)) ) <= eps:
			# lam=-1 means every point is touching
		else:
			# line is parallel to plane and does not lie on the plane
			return False, None, None, None

	lam = (nx*(x0-x1) + ny*(y0-y1) + nz*(z0-z1) ) / ( (nx*(x2-x1) + ny*(y2-y1) + nz*(z2-z1) ) )

	if 0<=lam<=1:
		# point_affect -> 0 means first 1 means second
		val = nx*(x1-x0) + ny*(y1-y0) + nz*(z1-z0)
		if val < 0:
			point_affect = 0
		else:
			point_affect = 1
		point_of_contact = np.array([
				x1 + lam*(x2-x1), 
				y1 + lam*(y2-y1), 
				z1 + lam*(z2-z1), 
			])
		return True, point_of_contact, plane[0], point_affect
	else:
		return False, None, None, None

def detectCollisions(obj_list):
	# there is a plane at z=0
	# with which we have to compute the collision dynamics
	plane = [np.array([0,0,1]), np.array([0,0,0])]
	point_of_contacts = []
	contact_normals = []
	obj_indices = []
	point_affect_list = []
	# plane is given by its normal vector and some point on the plane
	n=len(obj_list)
	for obj in obj_list:
		if obj.type == 'RevoluteJoint':
			isContact, point_of_contact, normal_contact, point_affect = detectCollision(obj, plane)
			if isContact:
				point_of_contacts.append(point_of_contact)
				contact_normals.append(normal_contact)
				obj_indices.append(obj.index)
				point_affect_list.append(point_affect)

	return obj_indices, point_of_contacts, contact_normals, point_affect_list

class World(object):
	"""docstring for World"""
	def __init__(self, obj_list):
		super(World, self).__init__()
		self.obj_list = obj_list
		# rendering options
		np.set_printoptions(precision=4)
		plt.ion()
		self.fig = plt.figure()
		self.ax = plt.gca(projection="3d")
		self.ax.set_xlim(0,10)
		self.ax.set_ylim(0,10)
		self.ax.set_zlim(0,10)
		self.fig.canvas.mpl_connect("close_event", lambda event: exit())
		self.origin = np.array([0,0,0,1])
		n = len(obj_list)
		self.q = np.zeros((n))
		self.dqdt = np.zeros((n))
		self.d2qdt2 = np.zeros((n))
		self.gravity = np.array([0,0,-10], dtype=np.float32)

		# for linear chain
		# TODO generalize it for general tree like structure
		for i in range(n):
			self.q[i] = obj_list[i].getQ()
			self.dqdt[i] = obj_list[i].getDq()
			self.d2qdt2[i] = obj_list[i].getD2q()

		# q = np.array([0,0,np.pi/2], dtype=np.float32)
		# self.q = np.array([0,0,3*np.pi/4], dtype=np.float32)
		self.q = np.array([0,0,np.pi/2 - np.pi/8,2*np.pi/8], dtype=np.float32)
		# self.q = np.array([0,0,0 ,2*np.pi/8], dtype=np.float32)
		# self.dqdt = np.array([0,0,1], dtype=np.float32)
		# self.dqdt = np.array([0,0,0,1], dtype=np.float32)
		# dqdt = np.array([2,0,1,0], dtype=np.float32)

		for i in range(n):
			self.obj_list[i].setQ(self.q[i] )
			self.obj_list[i].setDq(self.dqdt[i] )
			self.obj_list[i].setD2q(self.d2qdt2[i] )

		self.t = 0


	def set_gravity(self, gravity):
		self.gravity = gravity

	def set_link_origin(self, origin):
		self.origin = origin

	def advect1(self, torque, dt):
		obj_list = self.obj_list
		n = len(obj_list)

		# assigning it like this makses a shallow copy
		q = self.q.copy()
		dqdt = self.dqdt.copy()
		collision = True
		ke = 0
		Dq =  evaluate_Dq(obj_list)
		ke = (1/2) * np.transpose(dqdt) @ Dq @ dqdt
		pe = evaluate_potential_energy(obj_list, self.gravity)
		phi = evaluate_potential_energy_derivative(obj_list, self.gravity)
		# create phi array as derivative of pe with qi's
		Dq_derivative = evaluate_Dq_derivative(obj_list) # its is matrix of shape - nXnXn
		C = createCArray(Dq_derivative, dqdt)
		# create Mc matrix
		rhs = torque - phi - C
		rhs = rhs * dt

		def collision_response():
			vel_diff = np.zeros((n,))
			if collision:
				obj_indices, contact_points, contact_normals, point_affect_list = detectCollisions(obj_list)
				if len(contact_points) > 0:
					# collision is detected therefore apply collision constraints
					num_contact_point = len(contact_points)
					N_lcp = np.zeros((n,num_contact_point))
					coef_rest = 0.9
					for i in obj_indices: 
						obj = obj_list[i]
						list_index = obj_indices.index(i)
						contact_point = contact_points[list_index]
						contact_normal = contact_normals[list_index]
						point_affect = point_affect_list[list_index]
						Jvi = evaluateJacobian(obj_list, i, point_affect)
						# solve the constraint equation ---> using LCP
						# N_lcp = np.transpose(Jvi) @ contact_normal
						N_lcp_i = np.transpose(Jvi) @ contact_normal
						N_lcp[:,list_index] = N_lcp_i

					tau_star = rhs + Dq @ dqdt
					lcp_A1 = dt* (np.transpose(N_lcp) @ np.linalg.inv(Dq) @ N_lcp)
					lcp_q1 = np.transpose(N_lcp) @ np.linalg.inv(Dq) @ tau_star
					# other friction terms will also come for now ignoring
					lcp_A = np.zeros((2+num_contact_point,2+num_contact_point))
					lcp_A[0:num_contact_point,0:num_contact_point] = lcp_A1
					lcp_A[0+num_contact_point,0+num_contact_point] = 1
					lcp_A[1+num_contact_point,1+num_contact_point] = 1

					lcp_q = np.append(lcp_q1, [0, 0])
					sol = lcp.lemkelcp(lcp_A,lcp_q)
					fn = sol[0][0:num_contact_point] * (1+coef_rest)
					val = dt* N_lcp @ fn
					return np.reshape(val, (n,))

				else:
					return np.zeros((n,))
			else:
				return np.zeros((n,))

		cr = collision_response()
		rhs += cr
		new_momentum = rhs + Dq @ dqdt
		dqdt = np.linalg.solve(Dq, new_momentum )
		q += dqdt*dt

		self.q = q.copy()
		self.dqdt = dqdt.copy()

		t = self.t
		self.t += dt

		print("t: ", np.array([t]), "dqdt : ", dqdt)
		self.update()
		# Jci = evaluateJacobian(obj_list, 2, -1)
		# Jvi = evaluateJacobian(obj_list, 2, 1)
		# print("t: ", np.array([t]), "com vel : ", Jci@dqdt)
		# print("t: ", np.array([t]), "end vel : ", Jvi@dqdt)
		# print("t: ", np.array([t]), "q : ", q)
		print("t: ", np.array([t]), "energy : ", np.array([ke+pe]))
		# # print("t: ", np.array([t]), " Dq : ", '\n',  Dq)
		# # print("t: ", np.array([t]), "C : ", C)	
		print("t: ", np.array([t]), "Phi : ", phi)
		print("t: ", np.array([t]), "rhs : ", rhs)
		print()



	def advect(self, torque, dt):
		obj_list = self.obj_list
		n = len(obj_list)

		# assigning it like this makses a shallow copy
		q = self.q.copy()
		dqdt = self.dqdt.copy()
		d2qdt2 = self.d2qdt2.copy()
		collision = True
		ke = 0
		Dq =  evaluate_Dq(obj_list)
		ke = (1/2) * np.transpose(dqdt) @ Dq @ dqdt
		pe = evaluate_potential_energy(obj_list, self.gravity)
		phi = evaluate_potential_energy_derivative(obj_list, self.gravity)
		# create phi array as derivative of pe with qi's
		Dq_derivative = evaluate_Dq_derivative(obj_list) # its is matrix of shape - nXnXn
		C = createCArray(Dq_derivative, dqdt)
		# create Mc matrix
		rhs = torque - phi - C
		d2qdt2 = np.linalg.solve(Dq, rhs)
		
		# best integerating system as of now now just add more 
		# damping terms to kill of the motion rather than improving this further
		q += dqdt*dt + 0.5 * dt*dt*d2qdt2
		dqdt += d2qdt2*dt
		# q += dqdt*dt
	
		self.q = q.copy()
		self.dqdt = dqdt.copy()
		self.d2qdt2 = d2qdt2.copy()
		self.update()

		vel_diff = np.zeros((n,))
		if collision:
			obj_indices, contact_points, contact_normals, point_affect_list = detectCollisions(obj_list)
			if len(contact_points) > 0:
				# collision is detected therefore apply collision constraints
				m = n	
				K = np.zeros((m,m))
				coef_rest = 0.4
				acc_diff = np.zeros((m,))
				d2qdt2_end_effectors = calculate_joint_effector_accelerations(obj_list, d2qdt2)
				for i in range(m): # for second dof
					g_i = 1  # force/torque at ith joint
					test_force = rhs.copy()
					test_force[i] += g_i
					d2qdt2_test = np.linalg.solve(Dq, test_force)
					# apply this g_i force to calculate kij's
					for j in range(m):
						q_t = d2qdt2_test[j]
						q_0 = d2qdt2[j]
						k_ij = (q_t - q_0) / g_i
						K[i,j] = k_ij

					obj = obj_list[i]
					# since I know for now point of contact is only 1
					# TODO generalie this block
					if i in obj_indices:
						list_index = obj_indices.index(i)
						contact_point = contact_points[list_index]
						contact_normal = contact_normals[list_index]
						point_affect = point_affect_list[list_index]
						Jvi = evaluateJacobian(obj_list, i, point_affect)
						Jci = evaluateJacobian(obj_list, i, -1)
						# solve the constraint equation ---> using LCP
						print("Collision of object index ", i , "with ground at this point ", contact_point)
						print(point_affect)
						print(Jvi @ dqdt)

						impulse_lhs = np.transpose(contact_normal) @ Jvi @ dqdt
						if impulse_lhs<0:
							impulse_lhs = -(1+coef_rest)*impulse_lhs
						else:
							impulse_lhs = 0
						impulse_rhs = np.transpose(contact_normal) @ Jvi @ np.linalg.inv(Dq) @ np.transpose(Jvi)
						impulse = np.array([0,0,impulse_lhs/impulse_rhs[2]], dtype=np.float32)
						# for inelastic collision with e
						# pdb.set_trace()
						vel_diff = np.linalg.solve(Dq, np.transpose(Jvi)@ impulse )
						dqdt += vel_diff	
						print(impulse)
						print(vel_diff)
						print('vel of com', Jci @ vel_diff)
						# pdb.set_trace()

						# acc_diff = vel_diff / dt
						# f_generalized = np.linalg.solve(K,acc_diff)
						# rhs_new = rhs.copy()
						# rhs_new += f_generalized
						# d2qdt2 = np.linalg.solve(Dq, rhs_new)
					# new acceleration which will give this new velocity value


		# q += dqdt*dt
		# dqdt += d2qdt2*dt		
		# self.q = q.copy()
		self.dqdt = dqdt.copy()
		# self.d2qdt2 = d2qdt2.copy()
		# self.update()


		t = self.t
		self.t += dt

		# print("t: ", np.array([t]), " q : ", q, "dqdt : " , dqdt, "d2qdt2 : ", d2qdt2, "rhs : ", rhs, "ke : ", np.array([ke]), "pe : ", np.array([pe]))
		print("t: ", np.array([t]), "d2qdt2 : ", d2qdt2)
		print("t: ", np.array([t]), "dqdt : ", dqdt)
		self.update()
		Jci = evaluateJacobian(obj_list, 2, -1)
		Jvi = evaluateJacobian(obj_list, 2, 1)
		print("t: ", np.array([t]), "com vel : ", Jci@dqdt)
		print("t: ", np.array([t]), "end vel : ", Jvi@dqdt)
		# print("t: ", np.array([t]), "q : ", q)
		print("t: ", np.array([t]), "energy : ", np.array([ke+pe]))
		# # print("t: ", np.array([t]), " Dq : ", '\n',  Dq)
		# # print("t: ", np.array([t]), "C : ", C)	
		print("t: ", np.array([t]), "Phi : ", phi)
		print("t: ", np.array([t]), "rhs : ", rhs)
		print()
		# pdb.set_trace()

	def update(self):
		obj_list = self.obj_list
		n = len(obj_list)
		obj = obj_list[0]
		links = [obj]
		origin = self.origin
		# update all the links in top down fashion
		# and update the q values for local 
		while len(links) > 0:
			obj = links[0]
			links.remove(links[0])
			i= obj.getIndex()
			obj.setQ(self.q[i])
			if obj.Dq_flag:
				obj.setDq(self.dqdt[i])
				obj.setD2q(self.d2qdt2[i])	
			else:
				obj.Dq_flag = True

			obj.update(origin)

			if len(obj.child ) == 0:
				continue
			else:
				for i in range(len(obj.child)):
					links.append(obj.child[i])

		# update axes
		axis_x, axis_y, axis_z = np.array([[1,0,0],[0,1,0],[0,0,1]])
		for i in range(0,n):
			# matrix = obj_list[i].get_transformation_matrix()[0:3,0:3]
			# axis_x_new = rotate_axis_x with axis_z with theta_i
			# axis_z_new = rotate_axis_z with axis_x_new with alpha_i

			alpha = obj_list[i].getAlpha()
			theta = obj_list[i].getTheta()
			axis_x_new = axisAngleRotationMatrix(theta, axis_z) @ axis_x
			axis_z_new = axisAngleRotationMatrix(alpha, axis_x_new) @ axis_z
			axis_y_new = np.cross(axis_z_new, axis_x_new)


			axis_x = axis_x_new
			axis_y = axis_y_new
			axis_z = axis_z_new
			obj_list[i].set_xaxis(axis_x)
			obj_list[i].set_yaxis(axis_y)
			obj_list[i].set_zaxis(axis_z)
			# print("x_axis ", axis_x)
			# # print(axis_y)
			# print("z_axis ", axis_z)
			# print()
		# pdb.set_trace()

	def render(self):
		n = len(self.obj_list)
		# render objects
		for artist in plt.gca().lines + plt.gca().collections + plt.gca().texts:
			artist.remove()

		for i in range(n):
			self.obj_list[i].drawObject(self.ax)

		# debug code ----------------------------------
		O = np.zeros((n+1,3))
		Oc = np.zeros((3,))
		matrix = np.eye(4)
		mass = 0
		for i in range(1,n+1):
			matrix = matrix @ self.obj_list[i-1].get_transformation_matrix() 
			O[i] = self.origin[0:3] + (matrix @ np.append([O[0]],1))[0:3]
			if i==1:
				Oc += self.obj_list[i-1].getMass() * (O[i] + self.origin[0:3]) / 2.0
			else:
				Oc += self.obj_list[i-1].getMass() * (O[i] + O[i-1]) / 2.0
			mass += self.obj_list[i-1].getMass()
		# Oc = Oc/ mass
		Oc = (O[2] + O[3]) / 2 
		# Oc = ((O[2] + O[3]) / 2 + (O[3] + O[4]) / 2) / 2
			# ax.plot( [O[i-1][0], O[i][0]] ,[O[i-1][1], O[i][1]], [O[i-1][2], O[i][2]], marker = '*', color='red')
		self.ax.plot( [Oc[0]] ,[Oc[1]], [Oc[2]], marker = '*', color='red')

		# draw ground at xy plane
		xx, yy = np.meshgrid(range(0,11,10), range(0,11,10))
		zz = xx * 0
		self.ax.plot_surface(xx, yy, zz, alpha=0.2, color='black')

		# draw axis ->
		axis_x, axis_y, axis_z = np.zeros((3,3))
		axis_length = 10
		axis_u, axis_v, axis_w = np.array([[axis_length*1,0,0],[0,axis_length*1,0],[0,0,axis_length*1]])
		self.ax.quiver(axis_x,axis_y,axis_z,axis_u,axis_v,axis_w,arrow_length_ratio=0.1)

		self.fig.canvas.draw()
		self.fig.canvas.flush_events()