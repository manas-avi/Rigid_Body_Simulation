import numpy as np
import pdb
from special_matrices import *
from scipy.spatial.transform import Rotation as R

class PrismaticJoint(object):
	"""docstring for PrismaticJoint"""
	def __init__(self, name, x1, y1,z1, x2, y2,z2, q_angle, x_alpha, r_length, color=None, typ=None):
		super(PrismaticJoint, self).__init__()
		self.index = 0
		self.type = 'PrismaticJoint'
		self.name = name
		self.lineX, self.lineY, self.lineZ = [x1,x2] , [y1,y2], [z1,z2]

		self.x1_orignal = x1
		self.y1_orignal = y1
		self.z1_orignal = z1

		self.x2_orignal = x2
		self.y2_orignal = y2
		self.z2_orignal = z2

		# DH parameters
		self.q_angle = q_angle #theta parameter
		self.x_alpha = x_alpha # alpha parameter
		self.r_length = r_length  #r_length
		self.q = np.array([0]) # parameter z_length that varies


		# default axis
		self.x_axis = np.array([1,0,0])
		self.y_axis = np.array([0,1,0])
		self.z_axis = np.array([0,0,1])
	
		self.isMovable = True
		self.mass = 1.0

		self.Ic = np.zeros((3,3))
		# no rotation component associated

		self.showText = False
		self.showCom = True

		self.dq = np.array([0])
		self.d2q = np.array([0])

		self.rotation = np.eye(4)

		if color == None:
			self.color = np.random.rand(3,)
		else: 
			self.color = color

		self.parent = None
		self.child = []

		# flags
		self.Dq_flag = True


	def get_rotation_matrix(self):
		rotation_matrix = np.eye(4)
		theta = self.q_angle
		alpha = self.x_alpha

		rotation_matrix[0,0] = np.cos(theta)
		rotation_matrix[1,0] = np.sin(theta)

		rotation_matrix[0,1] = -np.sin(theta) * np.cos(alpha)
		rotation_matrix[1,1] =  np.cos(theta) * np.cos(alpha)
		rotation_matrix[2,1] =  np.sin(alpha)

		rotation_matrix[0,2] =  np.sin(theta) * np.sin(alpha)
		rotation_matrix[1,2] = -np.cos(theta) * np.sin(alpha)
		rotation_matrix[2,2] =  np.cos(alpha)

		return rotation_matrix

	def get_rotation_matrix_derivative(self):
		rotation_matrix = np.zeros((4,4))
		return rotation_matrix

	def get_translation_matrix(self):
		theta = self.q_angle
		a = self.r_length
		d = self.q[0]

		translation_matrix = np.eye(4)
		translation_matrix[0,3] = a * np.cos(theta)
		translation_matrix[1,3] = a * np.sin(theta)
		translation_matrix[2,3] = d

		return translation_matrix

	def get_transformation_matrix(self):
		return  self.get_translation_matrix() @ self.get_rotation_matrix()

	def update(self, origin):
		x1 = self.x1_orignal
		y1 = self.y1_orignal
		z1 = self.z1_orignal

		if self.parent == None:
			parent_rotation = np.eye(4)
			x1_new,y1_new,z1_new = x1,y1,z1
		else:
			parent_rotation = self.parent.rotation
			x1_new,y1_new,z1_new = np.array([self.parent.lineX[1],self.parent.lineY[1],self.parent.lineZ[1]]) 
			
		x1_new,y1_new,z1_new,_ = origin + parent_rotation @ (np.array([0,0,0,1]))
		x2_new,y2_new,z2_new,_ = origin + parent_rotation @ self.get_transformation_matrix() @ (np.array([0,0,0,1]))

		self.rotation = parent_rotation @ self.get_transformation_matrix()
		self.lineX, self.lineY, self.lineZ = [x1_new, x2_new] , [y1_new, y2_new], [z1_new, z2_new]

	def set_xaxis(self, x_axis):
		self.x_axis = x_axis

	def set_yaxis(self, y_axis):
		self.y_axis = y_axis

	def set_zaxis(self, z_axis):
		self.z_axis = z_axis

	def addChild(self, node):
		self.child.append(node)
		node.parent = self
		self.rotation = np.eye(4)


	def getIndex(self):
		return self.index

	def setIndex(self, idx):
		self.index = idx

	def setDq(self, dq):
		self.dq = np.array([dq])

	def getDq(self):
		return self.dq

	def setD2q(self, d2q):
		self.d2q = np.array([d2q])

	def getD2q(self):
		return self.d2q

	def setQ(self, q):
		self.q = np.array([q])

	def getQ(self):
		return self.q

	def getIc(self):
		return self.Ic

	def getMassMatrix(self):
		return self.massMatrix

	def setMass(self, mass):
		self.mass = mass

	def getMass(self):
		return self.mass

	def drawObject(self, ax):
		if self.isMovable and self.showText:
			w = self.dq[0]
			ax.text(self.lineX[0], self.lineY[0], self.lineZ[0], '  w:' + str(round(w, 4)) ) 
		ax.plot( self.lineX,self.lineY, self.lineZ, marker = 'o', color=self.color)

		if self.showCom:
			x1,x2 = self.lineX
			y1,y2 = self.lineY 
			z1,z2 = self.lineZ
			ax.plot( [(x1+x2)/2],[(y1+y2)/2], [(z1+z2)/2], marker = '*', color=self.color )

	def getLength(self):
		return np.sqrt( (self.lineX[0]-self.lineX[1])**2 + (self.lineY[0]-self.lineY[1])**2 
						+ (self.lineZ[0]-self.lineZ[1])**2 )

