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
# t_list=[]
# # e_list=[]
# q1_list=[]
# q2_list=[]

def skew_matrix(axis):
	x,y,z = axis
	mat = np.zeros((3,3))
	mat[1,0]=z
	mat[0,1]=-z
	mat[2,0]=-y
	mat[0,2]=y
	mat[1,2]=-x
	mat[2,1]=x
	return mat

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

def evaluateJacobian(obj_list, index):
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
				# dmatrix = dmatrix @ skew_matrix(z_axis_list[-1]) @ obj_list[i].get_rotation_matrix()[0:3,0:3]
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

	if ( (nx*(x2-x1) + ny*(y2-y1) + nz*(z2-z1) ) ) == 0:
		# line parralel tp the plane
		if (nx*(x0-x1) + ny*(y0-y1) + nz*(z0-z1) ) == 0:
			# line lies on the plane therefore inf points of contact. 
			# I will return just the end point (Remember)
			point_of_contact = np.array([x2,y2,z2])
			return True, point_of_contact, plane[0]
		else:
			# line is parallel to plane and does not lie on the plane
			return False, None, None

	lam = (nx*(x0-x1) + ny*(y0-y1) + nz*(z0-z1) ) / ( (nx*(x2-x1) + ny*(y2-y1) + nz*(z2-z1) ) )

	if 0<=lam<=1:
		point_of_contact = np.array([
				x1 + lam*(x2-x1), 
				y1 + lam*(y2-y1), 
				z1 + lam*(z2-z1), 
			])
		return True, point_of_contact, plane[0]
	else:
		return False, None, None

def detectCollisions(obj_list):
	# there is a plane at z=0
	# with which we have to compute the collision dynamics
	plane = [np.array([0,0,1]), np.array([0,0,0])]
	point_of_contacts = []
	contact_normals = []
	obj_indices = []
	# plane is given by its normal vector and some point on the plane
	n=len(obj_list)
	for obj in obj_list:
		if obj.type == 'RevoluteJoint':
			isContact, point_of_contact, normal_contact = detectCollision(obj, plane)
			if isContact:
				point_of_contacts.append(point_of_contact)
				contact_normals.append(normal_contact)
				obj_indices.append(obj.index)

	return obj_indices, point_of_contacts, contact_normals

def calculate_joint_effector_accelerations(obj_list, d2qdt2):
	a_list = []
	n = len(obj_list)
	for i in range(n):
		obj = obj_list[i]
		if i==0:
			z_axis = np.array([0,0,1])
		else:
			z_axis = obj_list[i-1].z_axis

		if obj.type == 'PrismaticJoint':
			acc = z_axis * d2qdt2[i]
			if not i==0:
				acc = a_list[-1] + acc
			a_list.append(acc)
		else:
			alpha = d2qdt2[i]
			a_dir = np.cross(z_axis, obj.getRadialVector()) * alpha
			if i==0:
				acc = a_dir
			else:
				acc = a_list[-1] + a_dir
			a_list.append(acc)
	return a_list


g = np.array([0,0,-10])
# g = np.array([0,-10,0])
# g = np.array([-10,0,0])
# g = np.array([-1,-1,-1])
# g = np.array([0,0,0])
# same axis of rotation
origin=np.array([4,4,4,4], dtype=np.float32)
# Always ensure that stating axis is always inclined with cartesian axis

# pobj_0.addChild(pobj_1)


pobj_0 = PrismaticJoint.PrismaticJoint('pobj_0', 4,4,4, 4,4,4, 	q_angle=np.pi/2, x_alpha=np.pi/2,r_length=0, color="red")
pobj_0.setMass(0)
pobj_1 = PrismaticJoint.PrismaticJoint('pobj_1', 4,4,4, 4,4,4, 	q_angle=np.pi/2, x_alpha=np.pi/2,r_length=0, color="blue")
pobj_1.setMass(0)
# pobj_2 = PrismaticJoint.PrismaticJoint('pobj_2', 4,4,4, 4,4,4, 	q_angle=np.pi/2, x_alpha=np.pi/2,r_length=0, color="green")
# pobj_2.setMass(0)
pobj_0.addChild(pobj_1)
# pobj_1.addChild(pobj_2)

obj_1 = RevoluteJoint.RevoluteJoint('obj_1', 4,4,4, 6,4,4, x_length=2, x_alpha=0,z_length=0, color="pink")
# obj_1 = RevoluteJoint.RevoluteJoint('obj_1', 4,4,4, 6,4,1, x_length=2, x_alpha=np.pi/2,z_length=0, color="violet")
obj_1.showText = True
pobj_1.addChild(obj_1)

obj_2 = RevoluteJoint.RevoluteJoint('obj_2', 6,4,4, 8,4,4, x_length=2, x_alpha=0,z_length=0, color="cyan")
obj_2.showText = True
obj_1.addChild(obj_2)


# ground is the x-y plane
# obj_list = [pobj_0]
obj_list = [pobj_0,pobj_1, obj_1, obj_2]
# obj_list = [pobj_0,pobj_1, obj_1]
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
	dt = 0.01
	q = np.zeros((n))
	dqdt = np.zeros((n))
	d2qdt2 = np.zeros((n))

	# intiatlising the matrices
	# for linear chain
	# TODO generalize it for general tree like structure
	t = 0

	for i in range(n):
		q[i] = obj_list[i].getQ()
		dqdt[i] = obj_list[i].getDq()
		d2qdt2[i] = obj_list[i].getD2q()

	# q = np.array([0,0,np.pi/2], dtype=np.float32)
	# q = np.array([0,0,0,0], dtype=np.float32)
	q = np.array([0,0,np.pi/2,np.pi/8], dtype=np.float32)
	# dqdt = np.array([0,0,-1,0], dtype=np.float32)
	# dqdt = np.array([-2,0,1,0], dtype=np.float32)

	for i in range(n):
		obj_list[i].setQ(q[i] )
		obj_list[i].setDq(dqdt[i] )
		obj_list[i].setD2q(d2qdt2[i] )



	# t_list=[]
	# e_list=[]
	# q_list=[]
	collision = True
	# collision = False
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

		rhs = rhs - phi - C
		# for k in range(n):
		# 	rhs[k] = rhs[k] - phi[k] - C[k]

		d2qdt2 = np.linalg.solve(Dq, rhs)

		# impulse_vector = np.concatenate([impulse, np.array([0,0,0])] ,axis=0)
		# if collision:
		# 	obj_indices, contact_points, contact_normals = detectCollisions(obj_list)
		# 	if len(contact_points) > 0:
		# 		# collision is detected therefore apply collision constraints
		# 		m = 4		
		# 		K = np.zeros((m,m))
		# 		coef_rest = 0.5
		# 		acc_diff = np.zeros((m,))
		# 		d2qdt2_end_effectors = calculate_joint_effector_accelerations(obj_list, d2qdt2)
		# 		for i in range(m): # for second dof
		# 			g_i = 1  # force/torque at ith joint
		# 			test_force = rhs.copy()
		# 			test_force[i] += g_i
		# 			d2qdt2_test = np.linalg.solve(Dq, test_force)
		# 			# apply this g_i force to calculate kij's
		# 			for j in range(m):
		# 				q_t = d2qdt2_test[j]
		# 				q_0 = d2qdt2[j]
		# 				k_ij = (q_t - q_0) / g_i
		# 				K[i,j] = k_ij

		# 			d2qdt2_end_effectors_test = calculate_joint_effector_accelerations(obj_list, d2qdt2_test)
		# 			obj = obj_list[i]
		# 			if i not in obj_indices:
		# 				acc_diff[i] = 0 # no change in acc
		# 			else:
		# 				list_index = obj_indices.index(i)
		# 				contact_point = contact_points[list_index]
		# 				contact_normal = contact_normals[list_index]
		# 				# solve the constraint equation ---> using LCP
		# 				print("Collision of object index ", i , "with ground at this point ", contact_point)
		# 				vels = calculate_joint_effector_velocities(obj_list)
		# 				vel = vels[i]
		# 				final_vel = -coef_rest * vel
		# 				init_vel = vel
		# 				final_acc = (final_vel - init_vel)/dt
		# 				final_acc_z = np.transpose(final_acc) @ contact_normal
		# 				# pdb.set_trace()
		# 				acc_diff[0] = final_acc_z - d2qdt2[0]
		# 				acc_diff[i] = 0
		# 		print(acc_diff, vel)
		# 		pdb.set_trace()
		# 		print(K)

		# 		f = np.linalg.solve(K,acc_diff)
		# 		# pdb.set_trace()
		# 		rhs_new = rhs.copy()
		# 		rhs_new[0:m] += f[0:m]
		# 		d2qdt2 = np.linalg.solve(Dq, rhs_new)

		# to constraint the root joint to have zero acceleration
		# q += dqdt*dt 
		q += dqdt*dt + 0.5 * dt*dt*d2qdt2
		dqdt += d2qdt2*dt
		t += dt

		# collision_dqdt = np.zeros((n,))
		# if collision:
		# 	obj_indices, contact_points, contact_normals = detectCollisions(obj_list)
		# 	coef_rest = 0.6
		# 	if len(contact_points) > 0:
		# 		# collision is detected therefore apply collision constraints
		# 		# for now this is hard coded since there is only 1 collision point_of_contact
		# 		# TODO correct it
		# 		i = obj_indices[0]
		# 		list_index = obj_indices.index(i)
		# 		Jvi = evaluateJacobian(obj_list, list_index)
		# 		contact_point = contact_points[list_index]
		# 		contact_normal = contact_normals[list_index]
		# 		# solve the constraint equation ---> using LCP
		# 		print("Collision of object index ", i , "with ground at this point ", contact_point)
		# 		vels = calculate_joint_effector_velocities(obj_list)
		# 		vel = vels[i]
		# 		final_vel = -coef_rest * vel
		# 		init_vel = vel
		# 		impulse = (final_vel - init_vel)
		# 		collision_dqdt = np.linalg.solve(Dq, np.transpose(Jvi)@impulse )

		dqdt += collision_dqdt


		obj = obj_list[0]
		links = [obj]
		# update all the links in top down fashion
		# and update the q values for local 
		while len(links) > 0:
			obj = links[0]
			links.remove(links[0])
			i= obj.getIndex()
			obj.setQ(q[i])
			if obj.Dq_flag:
				obj.setDq(dqdt[i])
				obj.setD2q(d2qdt2[i])	
			else:
				obj.Dq_flag = True

			obj.update(origin)

			# update axes
			axis_x, axis_y, axis_z = np.array([[1,0,0],[0,1,0],[0,0,1]])
			for i in range(0,n):
				matrix = obj_list[i].get_transformation_matrix()[0:3,0:3]
				axis_x = matrix @ axis_x
				axis_y = matrix @ axis_y
				axis_z = matrix @ axis_z 
				obj_list[i].set_xaxis(axis_x)
				obj_list[i].set_yaxis(axis_y)
				obj_list[i].set_zaxis(axis_z)

			if len(obj.child ) == 0:
				continue
			else:
				for i in range(len(obj.child)):
					links.append(obj.child[i])
		# pdb.set_trace()
		# print("t: ", np.array([t]), " q : ", q, "dqdt : " , dqdt, "d2qdt2 : ", d2qdt2, "rhs : ", rhs, "ke : ", np.array([ke]), "pe : ", np.array([pe]))
		print("t: ", np.array([t]), "d2qdt2 : ", d2qdt2)
		# print("t: ", np.array([t]), "d2qdt2_new : ", d2qdt2_new)
		print("t: ", np.array([t]), "dqdt : ", dqdt)
		# print("t: ", np.array([t]), "q : ", q)
		print("t: ", np.array([t]), "energy : ", np.array([ke+pe]))
		# print("t: ", np.array([t]), "Angular momentum : ", '\n',  Dq@dqdt, np.sum(Dq@dqdt))
		# # print("t: ", np.array([t]), " Dq : ", '\n',  Dq)
		# # print("t: ", np.array([t]), "C : ", C)	
		print("t: ", np.array([t]), "Phi : ", phi)
		print("t: ", np.array([t]), "rhs : ", rhs)
		# print("t: ", np.array([t]), "rhs_new : ", rhs_new)
		print()
		# pdb.set_trace()

		# debug statemetns
		# t_list.append(t)
		# q1_list.append(q[0])
		# q2_list.append(q[1])

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
		ax.plot_surface(xx, yy, zz, alpha=0.2, color='black')



		# draw axis ->
		axis_x, axis_y, axis_z = np.zeros((3,3))
		axis_length = 10
		axis_u, axis_v, axis_w = np.array([[axis_length*1,0,0],[0,axis_length*1,0],[0,0,axis_length*1]])
		ax.quiver(axis_x,axis_y,axis_z,axis_u,axis_v,axis_w,arrow_length_ratio=0.1)


		fig.canvas.draw()
		fig.canvas.flush_events()

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
