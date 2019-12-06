from mpl_toolkits.mplot3d import Axes3D
from matplotlib import patches
from matplotlib.collections import PatchCollection
import numpy as np
import RevoluteJoint
import PrismaticJoint
import matplotlib.pyplot as plt
import pdb
import lemkelcp as lcp

np.set_printoptions(suppress=True)


def createBasisVector(contact_normal, contact_point, num_vectors):
	t1 = np.cross(contact_normal, np.array([1,0,0])) # some direction 
	# since it is a cross product it will always be perpendicular to normal
	t2 = np.cross(t1, contact_normal) # some direction 
	t1 = t1 / np.linalg.norm(t1)
	t2 = t2 / np.linalg.norm(t2)

	vectors_list = []
	for i in range(num_vectors):
		theta_i = 2*np.pi*i/(num_vectors)
		t_i = np.cos(theta_i)*t2+ np.sin(theta_i)*t1
		vectors_list.append(t_i)
	vectors_list = np.transpose(np.array(vectors_list))
	return vectors_list


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

def get_origin_coordinates(obj_list):
	n = len(obj_list)
	O = np.zeros((n+1,3))
	for i in range(n):
		matrix = np.eye(4)
		node = obj_list[i]
		while not node is None:
			# pdb.set_trace()
			matrix = node.get_transformation_matrix() @ matrix
			node = node.parent
		O[i+1] = (matrix @ np.append([O[0]],1))[0:3]

	return O

def get_com_coordinates(obj_list):

	n = len(obj_list)
	Oc = np.zeros((3,))

	rci_list = []
	for i in range(n):
		O_parent = np.zeros((3,))
		matrix = np.eye(4)
		node = obj_list[i]
		nodes_matrix = node.get_transformation_matrix()
		while not node is None:
			matrix = node.get_transformation_matrix() @ matrix
			node = node.parent

		O_current = (matrix @ np.array([0,0,0,1]))[0:3]
		O_parent = (matrix @ np.linalg.inv(nodes_matrix) @ np.array([0,0,0,1]))[0:3]

		rci = (O_parent + O_current) / 2.0
		rci_list.append(rci)
	return rci_list

def get_index_list(obj_list, i):
	node = obj_list[i]
	i_list = []
	while not node is None:
		i_list.append(node.getIndex())
		node = node.parent
	i_list.reverse()
	return i_list

def evaluate_potential_energy(obj_list, g):
	n = len(obj_list)
	P = 0
	rci_list = get_com_coordinates(obj_list)
	for i in range(n):
		mi = obj_list[i].getMass()
		rci = rci_list[i]
		P-= (np.transpose(g) @ rci ) * mi
	return P

def evaluate_potential_energy_derivative(obj_list, g):
	n = len(obj_list)
	phi = np.zeros((n,))

	O = get_origin_coordinates(obj_list)
	# check if parent will come here nor not.......
	rci_list = get_com_coordinates(obj_list)
	# remember derivative will come only when part of the same chain
	for i in range(n):
		mi = obj_list[i].getMass()
		rci = rci_list[i]
		i_list = get_index_list(obj_list, i)
		for k in i_list:
			if k==0:
				Zkm1 = np.array([0,0,1])
			else:	
				Zkm1 = obj_list[k].parent.z_axis

			# k's parent index is different
			if k>0:
				kp_index = obj_list[k].parent.getIndex() + 1
			else:
				kp_index = 0

			if obj_list[k].type == 'PrismaticJoint':
				phi[k] -= mi * (np.transpose(g) @ Zkm1)
			elif obj_list[k].type == 'RevoluteJoint':				
				phi[k] -= mi * (np.transpose(g) @ np.cross( Zkm1 , (rci - O[kp_index]) ))
	return phi

def evaluateJacobian(obj_list, index, point_affect):
	# index is the end point of obj_list[index] at which we are computing the jacobian
	# store Oi values for all the link
	n=len(obj_list)
	O = get_origin_coordinates(obj_list)
	# 0th coordinate is origin which is zero

	i = index
	p_index = obj_list[index].parent.getIndex()
	mi = obj_list[i].getMass()
	Jvi = np.zeros((3,n))

	i_list = get_index_list(obj_list, i)

	if point_affect==1:
		Oci = O[i+1] #it has to be the end effector and not the mid-point
	elif point_affect==0:
		Oci = O[p_index+1] #it has to be the end effector and not the mid-point
	else:
		Oci = (O[p_index+1] + O[i+1]) / 2 # it has to be the mid-point as every point is touching

	for j in i_list:
		if j==0:
			Zjm1 = np.array([0,0,1])
		else:	
			Zjm1 = obj_list[j].parent.z_axis
		
		# j's parent index is different
		if j>0:
			jp_index = obj_list[j].parent.getIndex() + 1
		else:
			jp_index = 0

		if obj_list[j].type == 'PrismaticJoint':
			Jvi[:,j] = Zjm1
		elif obj_list[j].type == 'RevoluteJoint':
			Jvi[:,j] = np.cross(Zjm1, (Oci - O[jp_index]))
	return Jvi

def evaluate_Dq(obj_list):
	# this function is coorect as per unit test for 2 rods
	n = len(obj_list)
	Dq = np.zeros((n,n))

	# store Oi values for all the link
	O = get_origin_coordinates(obj_list)
	Oci_list = get_com_coordinates(obj_list)

	Ri = np.eye(3) # it grows with time
	for i in range(n):
		# KE part ---------------
		mi = obj_list[i].getMass()
		Jvi = np.zeros((3,n))
		Oci = Oci_list[i]
		i_list = get_index_list(obj_list, i)
		for j in i_list:
			if j==0:
				Zjm1 = np.array([0,0,1])
			else:	
				Zjm1 = obj_list[j].parent.z_axis

			# j's parent index is different
			if j>0:
				jp_index = obj_list[j].parent.getIndex() + 1
			else:
				jp_index = 0

			if obj_list[j].type == 'PrismaticJoint':
				Jvi[:,j] = Zjm1
			elif obj_list[j].type == 'RevoluteJoint':
				Jvi[:,j] = np.cross(Zjm1, (Oci - O[jp_index]))
		Dq += mi* (np.transpose(Jvi) @ Jvi)

		# RE part -----------------
		Ii = obj_list[i].Ic
		Ri = Ri @ obj_list[i].get_rotation_matrix()[0:3,0:3]	
		Jwi = np.zeros((3,n))
		for j in i_list:
			if obj_list[j].type == 'PrismaticJoint':
				Jwi[:,j] = np.array([0,0,0], dtype=np.float32)
			elif obj_list[j].type == 'RevoluteJoint':
				if j==0:
					Jwi[:,j] = np.array([0,0,1])
				else:
					Jwi[:,j] = obj_list[j].parent.z_axis
		Dq += np.transpose(Jwi) @ Ri @ Ii @ np.transpose(Ri) @ Jwi
	return Dq	


def evaluate_Dq_derivative(obj_list):
	# This derivative is good for analytic case for n=2 joints
	# I want to calculate the derivative of D[i,j] elem with resepect to Qk 
	n = len(obj_list)
	C = np.zeros((n,n,n))

	# store Oi values for all the link
	O = get_origin_coordinates(obj_list)
	Oci_list = get_com_coordinates(obj_list)

	for k in range(n):
		# here k represents qk element with which we are planning to differentiate
		if k==0:
			Zkm1 = np.array([0,0,1])
		else:	
			Zkm1 = obj_list[k].parent.z_axis
		
		# k's parent index is different
		if k>0:
			kp_index = obj_list[k].parent.getIndex() + 1
		else:
			kp_index = 0

		# compute gradient of z_axis with respect to the derivative of k
		dz_axis_list = []
		z_axis_list = []		

		R_list=[]
		dR_list=[]

		for i in range(n):
			R_list.append(np.eye(3))
			dR_list.append(np.eye(3))

		for i in range(n):
			obj = obj_list[i]
			par = obj.parent

			if i==k:
				dR_list[i] = obj.get_rotation_matrix_derivative()[0:3,0:3]
			else:
				dR_list[i] = obj.get_rotation_matrix()[0:3,0:3]
			R_list[i] = obj.get_rotation_matrix()[0:3,0:3]

			while not par is None:
				par_index = par.getIndex()
				if par_index==k:
					dR_list[i] = par.get_rotation_matrix_derivative()[0:3,0:3] @ dR_list[i]
				else:
					dR_list[i] = par.get_rotation_matrix()[0:3,0:3] @ dR_list[i]

				R_list[i] = par.get_rotation_matrix()[0:3,0:3] @ R_list[i]
				par = par.parent

			# defining the derivative of z_axis
			if i<k:
				dz_axis_list.append(np.array([0,0,0]))
			else:
				dz_axis_list.append(dR_list[i] @ np.array([0,0,1]))

			# defining the z_axis
			if i==0:
				z_axis_list.append(np.array([0,0,1]))
			else:
				z_axis_list.append(obj_list[i].parent.z_axis)

		for i in range(n):
			# this loops for different object Jvi that contribute to the Dq nXn matrix
			mi = obj_list[i].getMass()
			Jvi = np.zeros((3,n))
			dJvi = np.zeros((3,n))
			Oci = Oci_list[i]

			i_list = get_index_list(obj_list, i)
			for j in i_list:
				Zjm1 = z_axis_list[j]
				dZjm1 = dz_axis_list[j]
				# j's parent index is different
				if j>0:
					jp_index = obj_list[j].parent.getIndex() + 1
				else:
					jp_index = 0
				if obj_list[j].type == 'PrismaticJoint':
					Jvi[:,j] = Zjm1
					if k<j and k in i_list:
						if obj_list[k].type == 'RevoluteJoint':
							dJvi[:,j] = dZjm1
						# otherwise no effect on this axis if is prismatic joint
				elif obj_list[j].type == 'RevoluteJoint':
					Jvi[:,j] = np.cross(Zjm1, (Oci - O[jp_index]))
					if k<j and k in i_list:
						dZjm1 = dZjm1
						if obj_list[k].type == 'RevoluteJoint':
							dJvi[:,j] = np.cross(Zjm1, np.cross(Zkm1, (Oci - O[jp_index])) ) + np.cross(dZjm1, (Oci - O[jp_index]))
						# otherwise no effect on this axis if is prismatic joint
					elif k<i+1 and k in i_list:
						# no derivative of Zjm1 as k>j
						dZjm1 = dZjm1
						if obj_list[k].type == 'RevoluteJoint':
							dJvi[:,j] = np.cross(Zjm1, np.cross(Zkm1, (Oci - O[kp_index])) )
						# otherwise no effect on this axis if is prismatic joint
					else:
						dJvi[:,j] = np.zeros((3))

			# if k==2 or k==4:
			# 	pdb.set_trace()
			C[:,:,k] += mi * (np.transpose(Jvi) @ dJvi + np.transpose(dJvi) @ Jvi)
			# Dq += mi* (np.transpose(Jvi) @ Jvi)

			# RE part -----------------
			Ii = obj_list[i].Ic
			Ri = R_list[i]
			Rdi = dR_list[i]


			Jwi = np.zeros((3,n))
			dJwi = np.zeros((3,n))
			for j in range(i+1):
				if obj_list[j].type == 'PrismaticJoint':
					Jwi[:,j] = np.array([0,0,0], dtype=np.float32)
					dJwi[:,j] = np.array([0,0,0], dtype=np.float32)

				elif obj_list[j].type == 'RevoluteJoint':
					Jwi[:,j] = z_axis_list[j]
					dJwi[:,j] = dz_axis_list[j]

			if k<=i and k in i_list:
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
		self.ax.set_xlim(0,15)
		self.ax.set_ylim(0,15)
		self.ax.set_zlim(0,15)
		self.fig.canvas.mpl_connect("close_event", lambda event: exit())
		self.origin = np.array([0,0,0,1])
		n = len(obj_list)
		self.q = np.zeros((n))
		self.dqdt = np.zeros((n))
		self.d2qdt2 = np.zeros((n))
		self.gravity = np.array([0,0,-10], dtype=np.float32)
		self.friction = 1.0

		for i in range(n):
			self.q[i] = obj_list[i].getQ()
			self.dqdt[i] = obj_list[i].getDq()
			self.d2qdt2[i] = obj_list[i].getD2q()

		self.t = 0

		self.isCollision = True

	def set_gravity(self, gravity):
		self.gravity = gravity

	def set_friction(self, fr):
		self.friction = fr

	def set_link_origin(self, origin):
		self.origin = origin

	def setQ(self, q):
		self.q = q
		n = len(self.obj_list)
		for i in range(n):
			self.obj_list[i].setQ(self.q[i] )

	def setdqdt(self, dqdt):
		self.dqdt = dqdt
		n = len(self.obj_list)
		for i in range(n):
			self.obj_list[i].setDq(self.dqdt[i])

	def iadvect(self, qdot_des, dt, torque_list):
		on_ground = False
		obj_list = self.obj_list
		n = len(obj_list)

		# previous copies are used to compute the matrices ---> 
		q = self.q.copy()
		dqdt = self.dqdt.copy()
		collision = self.isCollision
		ke = 0
		Dq =  evaluate_Dq(obj_list)
		ke = (1/2) * np.transpose(dqdt) @ Dq @ dqdt
		pe = evaluate_potential_energy(obj_list, self.gravity)
		phi = evaluate_potential_energy_derivative(obj_list, self.gravity)
		# create phi array as derivative of pe with qi's
		Dq_derivative = evaluate_Dq_derivative(obj_list) # its is matrix of shape - nXnXn
		# C = createCArray(Dq_derivative, dqdt)
		C = createCArray(Dq_derivative, dqdt)

		tau_star = Dq@dqdt - dt*(C + phi)

		# selection matrix #P		
		num_dof = len(dqdt)
		num_actuated_dof = len(qdot_des)
		diff = num_dof - num_actuated_dof
		P = np.zeros((num_actuated_dof, num_dof))
		for i in range(num_actuated_dof):
			for j in range(num_actuated_dof):
				if i==j:
					P[i,j+diff] = 1
		torque = np.zeros((num_actuated_dof,))
		v_new = np.zeros((num_dof,))
		# with the constraint that P @ v_final = q_des

		# don't make it zero numerical unstability
		if collision:
			obj_indices, contact_points, contact_normals, point_affect_list = detectCollisions(obj_list)
			if len(contact_points) > 0:
				# collision is detected therefore apply collision constraints
				number_of_basis_vectors = 4
				coef_rest = 0.0
				coef_fric = self.friction
				num_contact_points = len(contact_points)
				d = number_of_basis_vectors
				p = num_contact_points
				N_lcp = np.zeros((n,p))
				B_lcp = np.zeros((n,d*p))
				E = np.zeros((p*d,p))

				for i in range(p):
					E[i*d:(i+1)*d,i] = np.ones((d,))

				pdb.set_trace()

				for i in obj_indices: 
					obj = obj_list[i]
					list_index = obj_indices.index(i)
					contact_point = contact_points[list_index]
					contact_normal = contact_normals[list_index]
					point_affect = point_affect_list[list_index]
					# chain_index = chain_index_list[list_index]
					Jvi = evaluateJacobian(obj_list, i, point_affect)
					# solve the constraint equation ---> using LCP
					N_lcp_i = np.transpose(Jvi) @ contact_normal
					N_lcp[:,list_index] = N_lcp_i

					D_i = createBasisVector(contact_normal, contact_point, d)
					B_lcp_i = np.transpose(Jvi) @ D_i
					d = number_of_basis_vectors
					B_lcp[:,d*list_index:d*(list_index+1)] = B_lcp_i

				lcp_A11 = dt* (np.transpose(N_lcp) @ np.linalg.inv(Dq) @ N_lcp)
				lcp_A12 = dt* (np.transpose(N_lcp) @ np.linalg.inv(Dq) @ B_lcp)
				lcp_A13 = np.zeros((p, p))

				lcp_A21 = dt* (np.transpose(B_lcp) @ np.linalg.inv(Dq) @ N_lcp)
				lcp_A22 = dt* (np.transpose(B_lcp) @ np.linalg.inv(Dq) @ B_lcp)
				lcp_A23 = E

				lcp_A31 = np.eye(p) * coef_fric
				lcp_A32 = -np.transpose(E)
				lcp_A33 = np.zeros((p, p))

				lcp_q1 = np.transpose(N_lcp) @ np.linalg.inv(Dq) @ tau_star
				lcp_q2 = np.transpose(B_lcp) @ np.linalg.inv(Dq) @ tau_star
				lcp_q3 = np.zeros((p,1))
				# other friction terms will also come for now ignoring
				dim = 2*p + p*d
				lcp_A = np.zeros((dim,dim))
				lcp_A[0:p,0:p] = lcp_A11
				lcp_A[0:p,p:p+p*d] = lcp_A12
				lcp_A[0:p,p+p*d:2*p+p*d] = lcp_A13

				lcp_A[p:p+p*d,0:p] = lcp_A21
				lcp_A[p:p+p*d,p:p+p*d] = lcp_A22
				lcp_A[p:p+p*d,p+p*d:2*p+p*d] = lcp_A23

				lcp_A[p+p*d:2*p+p*d,0:p] = lcp_A31
				lcp_A[p+p*d:2*p+p*d,p:p+p*d] = lcp_A32
				lcp_A[p+p*d:2*p+p*d,p+p*d:2*p+p*d] = lcp_A33

				lcp_q = np.zeros((2*p+p*d,1))
				lcp_q[0:p] = np.reshape(lcp_q1, (p,1))
				lcp_q[p:p+p*d] = np.reshape(lcp_q2, (p*d,1))
				lcp_q[p+p*d:2*p+p*d] = np.reshape(lcp_q3, (p,1))

				sol = lcp.lemkelcp(lcp_A,lcp_q)

				print(sol[2])
				if sol[2] == 'Solution Found':
					fn = sol[0][0:p] * (1+coef_rest)
					fd = sol[0][p:p+p*d]
					val = dt* (N_lcp @ fn + B_lcp@fd)
					return np.reshape(val, (n,)), True
				elif sol[2] == 'Secondary ray found':
					# it is the detaching case where there is contact but the collision is moving apart
					return np.zeros((n,)), True
				else:
					return np.zeros((n,)), True

			else:
				# since no collision takes place now here is how we will obtain the torque values
				# thus solution is similar to the case below
				v_final = np.zeros((n,))
				diff = num_dof- num_actuated_dof
				A_11 = Dq[0:diff,0:diff]
				A_12 = Dq[0:diff:,diff:]
				v_final_rest = qdot_des
				v_final_start = np.linalg.solve(A_11 , tau_star[0:diff] -A_12 @ v_final_rest)
				v_final[diff:,] = qdot_des
				v_final[0:diff,] = v_final_start
				torque = P @ Dq @ v_final - P @ tau_star
		else:
			# assuming we have switched off the collision interactoin
			# then we just need to solve two equations----
			# 1) P @ v_new = qdot_des
			# 2) Dq@v_new - Dq@v_old = dt*(np.tranpose(P)@tau - C - phi)
			v_final = np.zeros((n,))
			diff = num_dof- num_actuated_dof
			A_11 = Dq[0:diff,0:diff]
			A_12 = Dq[0:diff:,diff:]
			v_final_rest = qdot_des
			v_final_start = np.linalg.solve(A_11 , tau_star[0:diff] -A_12 @ v_final_rest)
			v_final[diff:,] = qdot_des
			v_final[0:diff,] = v_final_start
			torque = (P @ Dq @ v_final - P @ tau_star) / dt

		actual_torque = np.transpose(P) @ torque
		torque_list.append(actual_torque)
		dqdt = v_final
		q += dqdt*dt

		self.q = q.copy()
		self.dqdt = dqdt.copy()

		t = self.t
		self.t += dt

		print("t: ", np.array([t]), "torque : ",actual_torque)
		print("t: ", np.array([t]), "dqdt : ", dqdt)
		self.update()
		print("t: ", np.array([t]), "q : ", q)
		print("t: ", np.array([t]), "energy : ", np.array([ke+pe]))
		print("t: ", np.array([t]), "Phi : ", phi)
		print("t: ", np.array([t]), "C : ", C)
		print("t: ", np.array([t]), " Dq : ", '\n',  Dq)
		print()

		return torque_list

	def advect(self, torque, dt):
		on_ground = False
		obj_list = self.obj_list
		n = len(obj_list)

		# assigning it like this makses a shallow copy
		q = self.q.copy()
		dqdt = self.dqdt.copy()
		collision = self.isCollision
		ke = 0
		Dq =  evaluate_Dq(obj_list)
		ke = (1/2) * np.transpose(dqdt) @ Dq @ dqdt
		pe = evaluate_potential_energy(obj_list, self.gravity)
		phi = evaluate_potential_energy_derivative(obj_list, self.gravity)
		# create phi array as derivative of pe with qi's
		Dq_derivative = evaluate_Dq_derivative(obj_list) # its is matrix of shape - nXnXn
		C = createCArray(Dq_derivative, dqdt)

		rhs = torque - phi - C
		rhs = rhs * dt

		def collision_response():
			coef_rest = 0.0
			# coef_fric = 0.001
			# coef_fric = 1.19
			coef_fric = self.friction
			# because of numerics answer is slightly greater than 1
			# otherwise optimal value for mu is 1.0

			# don't make it zero numerical unstability
			number_of_basis_vectors = 4
			vel_diff = np.zeros((n,))
			tau_star = rhs + Dq @ dqdt
			if collision:
				obj_indices, contact_points, contact_normals, point_affect_list = detectCollisions(obj_list)
				if len(contact_points) > 0:
					# collision is detected therefore apply collision constraints
					num_contact_points = len(contact_points)
					d = number_of_basis_vectors
					p = num_contact_points
					N_lcp = np.zeros((n,p))
					B_lcp = np.zeros((n,d*p))
					E = np.zeros((p*d,p))
					for i in range(p):
						E[i*d:(i+1)*d,i] = np.ones((d,))

					for i in obj_indices: 
						obj = obj_list[i]
						list_index = obj_indices.index(i)
						contact_point = contact_points[list_index]
						contact_normal = contact_normals[list_index]
						point_affect = point_affect_list[list_index]
						# chain_index = chain_index_list[list_index]
						Jvi = evaluateJacobian(obj_list, i, point_affect)
						# solve the constraint equation ---> using LCP
						N_lcp_i = np.transpose(Jvi) @ contact_normal
						N_lcp[:,list_index] = N_lcp_i

						D_i = createBasisVector(contact_normal, contact_point, d)
						B_lcp_i = np.transpose(Jvi) @ D_i
						d = number_of_basis_vectors
						B_lcp[:,d*list_index:d*(list_index+1)] = B_lcp_i

					lcp_A11 = dt* (np.transpose(N_lcp) @ np.linalg.inv(Dq) @ N_lcp)
					lcp_A12 = dt* (np.transpose(N_lcp) @ np.linalg.inv(Dq) @ B_lcp)
					lcp_A13 = np.zeros((p, p))

					lcp_A21 = dt* (np.transpose(B_lcp) @ np.linalg.inv(Dq) @ N_lcp)
					lcp_A22 = dt* (np.transpose(B_lcp) @ np.linalg.inv(Dq) @ B_lcp)
					lcp_A23 = E

					lcp_A31 = np.eye(p) * coef_fric
					lcp_A32 = -np.transpose(E)
					lcp_A33 = np.zeros((p, p))

					lcp_q1 = np.transpose(N_lcp) @ np.linalg.inv(Dq) @ tau_star
					lcp_q2 = np.transpose(B_lcp) @ np.linalg.inv(Dq) @ tau_star
					lcp_q3 = np.zeros((p,1))
					# other friction terms will also come for now ignoring
					dim = 2*p + p*d
					lcp_A = np.zeros((dim,dim))
					lcp_A[0:p,0:p] = lcp_A11
					lcp_A[0:p,p:p+p*d] = lcp_A12
					lcp_A[0:p,p+p*d:2*p+p*d] = lcp_A13

					lcp_A[p:p+p*d,0:p] = lcp_A21
					lcp_A[p:p+p*d,p:p+p*d] = lcp_A22
					lcp_A[p:p+p*d,p+p*d:2*p+p*d] = lcp_A23

					lcp_A[p+p*d:2*p+p*d,0:p] = lcp_A31
					lcp_A[p+p*d:2*p+p*d,p:p+p*d] = lcp_A32
					lcp_A[p+p*d:2*p+p*d,p+p*d:2*p+p*d] = lcp_A33

					lcp_q = np.zeros((2*p+p*d,1))
					lcp_q[0:p] = np.reshape(lcp_q1, (p,1))
					lcp_q[p:p+p*d] = np.reshape(lcp_q2, (p*d,1))
					lcp_q[p+p*d:2*p+p*d] = np.reshape(lcp_q3, (p,1))

					sol = lcp.lemkelcp(lcp_A,lcp_q)

					print(sol[2])
					if sol[2] == 'Solution Found':
						fn = sol[0][0:p] * (1+coef_rest)
						fd = sol[0][p:p+p*d]
						val = dt* (N_lcp @ fn + B_lcp@fd)
						return np.reshape(val, (n,)), True
					elif sol[2] == 'Secondary ray found':
						# it is the detaching case where there is contact but the collision is moving apart
						return np.zeros((n,)), True
					else:
						return np.zeros((n,)), True

				else:
					return np.zeros((n,)), False
			else:
				return np.zeros((n,)), False

		cr, on_ground = collision_response()
		rhs += cr
		new_momentum = rhs + Dq @ dqdt
		dqdt = np.linalg.solve(Dq, new_momentum )
		q += dqdt*dt

		self.q = q.copy()
		self.dqdt = dqdt.copy()

		t = self.t
		self.t += dt

		print("t: ", np.array([t]), "torque : ", torque)
		print("t: ", np.array([t]), "dqdt : ", dqdt)
		self.update()
		# Jci = evaluateJacobian(obj_list, 2, -1)
		# Jvi = evaluateJacobian(obj_list, 2, 1)
		# print("t: ", np.array([t]), "com vel : ", Jci@dqdt)
		# print("t: ", np.array([t]), "end vel : ", Jvi@dqdt)
		# print("t: ", np.array([t]), "pe : ", np.array([pe]))
		# print("t: ", np.array([t]), "ke : ", np.array([ke]))
		print("t: ", np.array([t]), "q : ", q)
		print("t: ", np.array([t]), "C : ", C)	
		print("t: ", np.array([t]), "Phi : ", phi)
		print("t: ", np.array([t]), "energy : ", np.array([ke+pe]))
		print("t: ", np.array([t]), " Dq : ", '\n',  Dq)
		# print("t: ", np.array([t]), "Dq_2 : \n", Dq_derivative[:,:,2])	
		# print("t: ", np.array([t]), "Dq_4 : \n", Dq_derivative[:,:,4])	
		# print("t: ", np.array([t]), "rhs : ", rhs)
		print()

		return on_ground


	def update_node(self, obj):
		i = obj.getIndex()
		obj.setQ(self.q[i])
		obj.setDq(self.dqdt[i])
		obj.setD2q(self.d2qdt2[i])	
		obj.update(self.origin)

		for child in obj.child:
			self.update_node(child)

	def update_axis(self, obj, paxis_x, paxis_y, paxis_z):
		alpha = obj.getAlpha()
		theta = obj.getTheta()
		axis_x = axisAngleRotationMatrix(theta, paxis_z) @ paxis_x
		axis_z = axisAngleRotationMatrix(alpha, axis_x) @ paxis_z
		axis_y = np.cross(axis_z, axis_x)

		obj.set_xaxis(axis_x)
		obj.set_yaxis(axis_y)
		obj.set_zaxis(axis_z)

		for child in obj.child:
			self.update_axis(child, axis_x, axis_y, axis_z)

	def update(self):
		obj_list = self.obj_list
		n = len(obj_list)
		obj = obj_list[0]

		# dfs to update
		self.update_node(obj)

		# update axes for a tree like structure in dfs format
		axis_x, axis_y, axis_z = np.array([[1,0,0],[0,1,0],[0,0,1]])
		self.update_axis(obj, axis_x, axis_y, axis_z)

	def render(self):
		n = len(self.obj_list)
		# render objects
		for artist in plt.gca().lines + plt.gca().collections + plt.gca().texts:
			artist.remove()

		for i in range(n):
			self.obj_list[i].drawObject(self.ax)

		# debug code ----------------------------------

		# O = np.zeros((n+1,3))
		# Oc = np.zeros((3,))
		# mass = 0
		# for i in range(0,n):

		# 	matrix = np.eye(4)
		# 	node = self.obj_list[i]
		# 	while not node is None:
		# 		matrix = node.get_transformation_matrix() @ matrix
		# 		node = node.parent
		# 	# matrix = matrix @ self.obj_list[i].get_transformation_matrix() 
		# 	O[i+1] = self.origin[0:3] + (matrix @ np.append([O[0]],1))[0:3]
		# 	if i==0:
		# 		Oc += self.obj_list[i].getMass() * (O[i+1] + self.origin[0:3]) / 2.0
		# 	else:
		# 		Oc += self.obj_list[i].getMass() * (O[i] + O[i+1]) / 2.0
		# 	mass += self.obj_list[i].getMass()


		# O = self.origin[0:3] + get_origin_coordinates(self.obj_list)
		# Oc = Oc/ mass
		# Oc = (O[2] + O[3]) / 2 
		# Oc = ((O[2] + O[3]) / 2 + (O[3] + O[4]) / 2) / 2
		# self.ax.plot( O[:,0] ,O[:,1], O[:,2], marker = '*', color='red')
		# self.ax.plot( [O[i-1][0], O[i][0]] ,[O[i-1][1], O[i][1]], [O[i-1][2], O[i][2]], marker = '*', color='red')
		# self.ax.plot( [Oc[0]] ,[Oc[1]], [Oc[2]], marker = '*', color='red')

		# draw ground at xy plane
		xx, yy = np.meshgrid(range(0,16,15), range(0,16,15))
		zz = xx * 0
		self.ax.plot_surface(xx, yy, zz, alpha=0.2, color='black')

		# draw axis ->
		axis_x, axis_y, axis_z = np.zeros((3,3))
		axis_length = 15
		axis_u, axis_v, axis_w = np.array([[axis_length*1,0,0],[0,axis_length*1,0],[0,0,axis_length*1]])
		self.ax.quiver(axis_x,axis_y,axis_z,axis_u,axis_v,axis_w,arrow_length_ratio=0.1)

		self.fig.canvas.draw()
		self.fig.canvas.flush_events()