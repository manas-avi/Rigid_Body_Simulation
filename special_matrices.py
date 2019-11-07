import numpy as np

def jacobian_matrix(r, theta):
	mat = np.eye(2)
	mat[0,0] = np.cos(theta)
	mat[0,1] = -r * np.sin(theta)
	mat[1,0] = np.sin(theta)
	mat[1,1] = r * np.cos(theta)
	return mat

def djacobian_matrix(r, theta, rvel, w):
	mat1 = np.eye(2)
	mat1[0,0] = 0
	mat1[0,1] = -np.sin(theta)
	mat1[1,0] = 0
	mat1[1,1] = np.cos(theta)
	mat1 = mat1 * rvel

	mat2 = np.eye(2)
	mat2[0,0] = -np.sin(theta)
	mat2[0,1] = -r * np.cos(theta)
	mat2[1,0] = np.cos(theta)
	mat2[1,1] = -r * np.sin(theta)
	mat2 = mat2 * w

	mat = mat1 + mat2
	return mat


def jacobian_matrix3(x, y, theta):
	mat = np.eye(3)
	return mat

def djacobian_matrix3(x, y, theta, vx, vy, w):
	return np.zeros((3,3))

def getRotationMatrix2X2(theta):
	mat = np.eye(3)
	mat[0,0] = np.cos(theta)
	mat[0,1] = -np.sin(theta)
	mat[1,0] = np.sin(theta)
	mat[1,1] = np.cos(theta)
	return mat

def getTranslationMatrix2X2(x,y):
	mat = np.eye(3)
	mat[0,2] = x
	mat[1,2] = y
	return mat

def getTranslationMatrix4X4(x,y,z):
	mat = np.eye(4)
	mat[0,3] = x
	mat[1,3] = y
	mat[2,3] = z
	return mat

def getRotationMatrix4X4(q):
	theta = q
	mat = np.eye(4)
	mat[0,0] = np.cos(theta)
	mat[0,1] = -np.sin(theta)
	mat[1,0] = np.sin(theta)
	mat[1,1] = np.cos(theta)

	return mat

	# TODO write the out of plane rotation also 

	# new_mat = np.eye(4)
	# new_mat[0:3,0:3] = mat
	# return new_mat

def getVector3X(vector4X):
	return np.array([vector4X[0], vector4X[1], vector4X[2]]) / vector4X[3]

def runge_kutta_method(v, vdot, dt):
	pass