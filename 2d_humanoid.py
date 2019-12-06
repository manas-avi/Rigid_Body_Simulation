import numpy as np
import RevoluteJoint
import PrismaticJoint
import World
import matplotlib.pyplot as plt
import pdb
import os


def readFile(filename):

	f = open(filename)
	lines = f.readlines()
	new_lines = []
	for line in lines:
		line = line.split(',')
		new_elems = []
		for elem in line:
			new_elems.append(float(elem))
		new_lines.append(new_elems)
	new_lines = np.array(new_lines)
	return new_lines


def interpolate_frames(base_frames, dt, method='linear'):
	new_frames = []
	num_frames = len(base_frames)
	num_int_frames = np.int(1/dt)
	for i in range(num_frames-1):
		int_frames = []
		pose_1 = base_frames[i]
		pose_2 = base_frames[i+1]
		for j in range(num_int_frames):
			int_frame = ( pose_2*j + pose_1*(num_int_frames-j)) / num_int_frames
			int_frames.append(int_frame)
		new_frames = new_frames + int_frames
	return new_frames

def interpolate_velocities(frames, dt):
	n = len(frames)
	dof = frames[0].shape[0]
	vels = [np.zeros((dof,))]
	for i in range(n-1):
		vel = (frames[i+1] - frames[i])
		vels.append(vel)
	vels[0] = vels[1]
	return vels

def interpolate_acceleration(frames, dt):
	n = len(frames)
	dof = frames[0].shape[0]
	accs = [np.zeros((dof,))]
	for i in range(n-1):
		# acc = (frames[i+1] - frames[i]) / dt
		acc = frames[i+1] - frames[i]
		accs.append(acc)
	return accs

g = np.array([0,0,-10])
# g = np.array([0,0,0])
# same axis of rotation
origin=np.array([6,4,8,4], dtype=np.float32)
# Always ensure that stating axis is always inclined with cartesian axis

# pobj_0.addChild(pobj_1)


pobj_0 = PrismaticJoint.PrismaticJoint('pobj_0', 6,4,8, 6,4,8,q_angle=np.pi/2, x_alpha=np.pi/2,r_length=0, color="red")
pobj_0.setMass(0)
pobj_0.showText = True
pobj_1 = PrismaticJoint.PrismaticJoint('pobj_1', 6,4,8, 6,4,8,q_angle=np.pi/2, x_alpha=np.pi/2,r_length=0, color="blue")
pobj_1.setMass(0)
pobj_0.addChild(pobj_1)

lobj_1 = RevoluteJoint.RevoluteJoint('lobj_1', 6,4,8, 10,4,8, x_length=4, x_alpha=0,z_length=0, color="pink")
lobj_1.showText = True
pobj_1.addChild(lobj_1)

lobj_2 = RevoluteJoint.RevoluteJoint('lobj_2', 10,4,8, 14,4,8, x_length=4, x_alpha=0,z_length=0, color="cyan")
lobj_2.showText = True
lobj_1.addChild(lobj_2)

lobj_3 = RevoluteJoint.RevoluteJoint('lobj_3', 14,4,8, 15,4,8, x_length=1, x_alpha=0,z_length=0, color="orange")
lobj_3.showText = True
lobj_2.addChild(lobj_3)

robj_1 = RevoluteJoint.RevoluteJoint('robj_1', 6,4,8, 10,4,8, x_length=4, x_alpha=0,z_length=0, color="pink")
robj_1.showText = True
pobj_1.addChild(robj_1)

robj_2 = RevoluteJoint.RevoluteJoint('robj_2', 10,4,8, 14,4,8, x_length=4, x_alpha=0,z_length=0, color="cyan")
robj_2.showText = True
robj_1.addChild(robj_2)

robj_3 = RevoluteJoint.RevoluteJoint('robj_3', 14,4,8, 15,4,8, x_length=1, x_alpha=0,z_length=0, color="orange")
robj_3.showText = True
robj_2.addChild(robj_3)

# ground is the x-y plane
# obj_list = [pobj_0,pobj_1, lobj_1, lobj_2, lobj_3]
# obj_list = [pobj_0,pobj_1, lobj_1, lobj_2, robj_1, robj_2]
obj_list = [pobj_0,pobj_1, lobj_1, lobj_2, lobj_3, robj_1, robj_2, robj_3]
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
	dt = 0.01*5
	world = World.World(obj_list)

	# rest pose
	# world.setQ(np.array([0,0,np.pi,0,-np.pi/2,np.pi,0,-np.pi/2], dtype=np.float32))

	# contact pose
	# world.setQ(np.array([0,0,5*np.pi/6,np.pi/8,-3*np.pi/5,   21*np.pi/20, np.pi/10,-np.pi/2], dtype=np.float32))	

	# down pose
	# world.setQ(np.array([0,0,3*np.pi/4,np.pi/3,-np.pi/2,    19*np.pi/20,2*np.pi/5,-2*np.pi/5], dtype=np.float32))	

	# passing pose
	# world.setQ(np.array([0,0,5*np.pi/6,np.pi/3,-np.pi/3,    21*np.pi/20,np.pi/20,-3*np.pi/5], dtype=np.float32))	

	# up pose
	# world.setQ(np.array([0,0,7*np.pi/10,np.pi/2,-np.pi/3,    21*np.pi/20,np.pi/20,-2*np.pi/5], dtype=np.float32))	

	rest_pose = np.array([0,0,np.pi,0,-np.pi/2,np.pi,0,-np.pi/2], dtype=np.float32)
	contact_pose = np.array([0,0,5*np.pi/6,np.pi/8,-3*np.pi/5,   11*np.pi/10, np.pi/10,-np.pi/2], dtype=np.float32)
	down_pose = np.array([0,0,3*np.pi/4,np.pi/3,-np.pi/2,    19*np.pi/20,2*np.pi/5,-2*np.pi/5], dtype=np.float32)
	passing_pose = np.array([0,0,5*np.pi/6,np.pi/3,-np.pi/3,    21*np.pi/20,np.pi/20,-3*np.pi/5], dtype=np.float32)
	up_pose = np.array([0,0,21*np.pi/20,np.pi/20,-2*np.pi/5,    7*np.pi/10,np.pi/2,-np.pi/3], dtype=np.float32)

	rcontact_pose = np.array([0,0,11*np.pi/10, np.pi/10,-np.pi/2,  5*np.pi/6,np.pi/8,-3*np.pi/5], dtype=np.float32)
	rdown_pose = np.array([0,0,19*np.pi/20,2*np.pi/5,-2*np.pi/5  ,3*np.pi/4,np.pi/3,-np.pi/2], dtype=np.float32)
	rpassing_pose = np.array([0,0,21*np.pi/20,np.pi/20,-3*np.pi/5,  5*np.pi/6,np.pi/3,-np.pi/3], dtype=np.float32)
	rup_pose = np.array([0,0,7*np.pi/10,np.pi/2,-np.pi/3,  21*np.pi/20,np.pi/20,-2*np.pi/5], dtype=np.float32)

	base_frames = [contact_pose, down_pose, passing_pose, up_pose, rcontact_pose,
					 rdown_pose, rpassing_pose, rup_pose, contact_pose]


	new_frames = interpolate_frames(base_frames, dt)
	vel_frames = interpolate_velocities(new_frames, dt)
	new_vel_frames = []
	for f in vel_frames:
		f[1] = dt
		new_vel_frames.append(f)

	vel_frames = new_vel_frames
	# pdb.set_trace()
	# this trajectory for all purposes is good enough
	acc_frames = interpolate_acceleration(vel_frames, dt)
	world.setQ(contact_pose)

	# world.setQ(passing_pose)
	# world.setQ(np.array([0,0,np.pi-np.pi/8,0, np.pi/8, np.pi+np.pi/8,0, -np.pi/8], dtype=np.float32))
	# dof_values = readFile('frame_values.txt')
	# world.setQ(dof_values[1])
	# world.setQ(np.array([0,0,np.pi-np.pi/8,np.pi/8,0, np.pi,0,0], dtype=np.float32))

	# world.setdqdt(np.array([1,5,1,1,1], dtype=np.float32))
	# world.setdqdt(np.array([0,0,1,-1,-1,1], dtype=np.float32))
	# world.setdqdt(np.array([0,0,-1,1,-1,1,-1,1], dtype=np.float32))
	# world.setdqdt(np.array([0,0,1,-1,1,-1,1,-1], dtype=np.float32))

	world.set_link_origin(origin)
	world.set_gravity(g)
	world.isCollision = False


	# t_list=[]
	# e_list=[]
	# q_list=[]

	num_poses = 8

	on_ground = False
	switch = 0
	frame_count = 0
	frame_rate = 1
	fps = 20
	frame_index = 0
	torque_list = []
	save_torque_value = True
	while True:
		# animate interpolated frames
		world.setQ(new_frames[frame_index])
		if 'torque_list.npy' in os.listdir('./'):
			torque_list = np.load('torque_list.npy')
			torque = torque_list[frame_index]
			world.advect(torque, dt)

		else:
			torque_list = world.iadvect(vel_frames[frame_index][2:], dt, torque_list)
		frame_index = (frame_index+1) % len(new_frames)
		world.update()

		# debug statemetns
		# t_list.append(world.t)
		# e_list.append(energy)

		world.render()

		if world.t>10:
			if save_torque_value:
				np.save('torque_list.npy', np.array(torque_list))
			break

	



exit()
# plt.title('energy graph')
# plt.xlabel('time ---->')
# plt.ylabel('energy ---->')
# # plt.plot(t_list, e_list)
# plt.plot(t_list, e_list, color='r')
# # plt.plot(t_list, q2_list, color='b')
# plt.show()


# ####### TRASH CODE ###############

if on_ground:
	on_ground = world.advect(torque, dt/4)
else:
	on_ground = world.advect(torque, dt)

# # animate world with constant poses
# if switch == 0: 
# 	pose = contact_pose
# elif switch == 1: 
# 	pose = down_pose
# elif switch == 2: 
# 	pose = passing_pose
# elif switch == 3: 
# 	pose = up_pose
# elif switch == 4: 
# 	pose = rcontact_pose
# elif switch == 5: 
# 	pose = rdown_pose
# elif switch == 6: 
# 	pose = rpassing_pose
# elif switch == 7: 
# 	pose = rup_pose			


# frame_count+=frame_rate
# if frame_count>fps:
# 	world.setQ(pose)
# 	switch = (switch + 1) % num_poses
# 	frame_count = 0
# 	print('new switch is ', switch)