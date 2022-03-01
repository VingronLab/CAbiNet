import numpy #as np
#import numpy.random
import torch
from scipy.sparse.linalg import eigsh


def svd_torch(x):
	x2 = torch.from_numpy(numpy.array(x))
	u, s, v = torch.svd(x2)
	u1 = u.numpy()
	s1 = s.numpy()
	v1 = v.numpy()
	return u1,s1,v1

def eig_torch(x):
	x1 = torch.from_numpy(numpy.array(x))
	e, u = torch.eig(x1,eigenvectors=True)
	e = e.numpy()
	u = u.numpy()
	return e, u

def eigsh_scipy(x, k=6,which = 'SM'):
	e,u = eigsh(x, k = int(k), which = which)
	return e, u
