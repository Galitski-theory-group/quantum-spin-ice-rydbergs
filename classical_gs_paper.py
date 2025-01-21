from math import *
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


e0 = 1/sqrt(2)*np.array([1,1,1]); 
e1 = 1/sqrt(2)*np.array([1,-1,-1]);
e2 = 1/sqrt(2)*np.array([-1,1,-1]);
e3 = 1/sqrt(2)*np.array([-1,-1,1]);
# with these vectors, the side of the tetrahedron is 1.
e = [e0,e1,e2,e3];

A1 = e0-e1;
A2 = e0-e2;
A3 = e0-e3;
A = [A1,A2,A3]

a = np.linalg.norm(e0/2-e1/2)
print("a = ",a)
vol = np.dot(A1,np.cross(A2,A3))

B1 = 2* pi* np.cross(A2,A3) / vol
B2 = 2* pi* np.cross(A3,A1) / vol
B3 = 2* pi* np.cross(A1,A2) / vol
B = [B1,B2,B3]

eps = np.array([1,1,-1,-1])

N = 30
Nuc = N**3

#alpha = np.array([0,0,0]) # alpha[i] is in {-0.5,0.5}
#alpha1 = alpha[0]
#alpha2 = alpha[1]
#alpha3 = alpha[2]

V = np.zeros((4,4))
for r1 in range(-N,N):
    for r2 in range(-N,N):
        for r3 in range(-N,N):
            r = r1* A1 + r2* A2 + r3* A3
            for i in range(4):
                for j in range(4):
                    rel_r = (e[i] - e[j])/2 + r
                    if rel_r[0] == 0 and rel_r[1] == 0 and rel_r[2] == 0:
                        continue
                    else :
                        V[i][j] += 1/( rel_r[0]**2 + rel_r[1]**2 + rel_r[2]**2 )**3 
                
V = V / 2

print("V = ")
print(V)
temp, temp2 =  np.linalg.eig(V)

print("evls = ",temp)
print("evecs = ", temp2)
print(np.transpose(temp2)[0])
print(sum(np.transpose(temp2)[1]))
print(sum(np.transpose(temp2)[2]))
print(sum(np.transpose(temp2)[3]))








