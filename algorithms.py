import numpy as np
import pdb

def runge_kutta_method(y, yhat, dt):
	return y+yhat*dt