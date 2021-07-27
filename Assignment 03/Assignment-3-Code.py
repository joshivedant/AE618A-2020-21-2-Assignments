# IEN Matrix

IEN = [[ 0 for _ in range(336)] for _ in range(4)]

ien = 1
for e in range(120):
    if (ien%21)==0:
        ien += 1
    IEN[0][e] = ien
    IEN[1][e] = ien+1
    IEN[2][e] = ien+22
    IEN[3][e] = ien+21
    ien += 1

for e in range(120,216):
    if (e==126):
        ien += 7
    if ((e-120)%6)==0:
        ien += 1
    IEN[0][e] = ien
    IEN[1][e] = ien+1
    if e<126 or e>=210:
        IEN[2][e] = ien+22
        IEN[3][e] = ien+21
    else:
        IEN[2][e] = ien+15
        IEN[3][e] = ien+14
    ien += 1

for e in range(216, 336):
    if ((e-216)%20)==0:
        ien += 1
    IEN[0][e] = ien
    IEN[1][e] = ien+1
    IEN[2][e] = ien+22
    IEN[3][e] = ien+21
    ien += 1

# ID Array

ID = [-1 for _ in range(392)]

boundary = list(range(1, 23)) + [42, 43, 63, 64, 84, 85, 105, 106, 126, 127, 133, 134, 135, 136, 137, 138, 139, 140, 141, 147, 148, 154, 155, 161, 162, 168, 169, 175, 176, 182, 183, 189, 190, 196, 197, 203, 204, 210, 211, 217, 218, 224, 225, 231, 232, 238, 239, 245, 246, 252, 253, 254, 255, 256, 257, 258, 259, 260, 266, 267, 287, 288, 308, 309, 329, 330, 350, 351] + list(range(371, 393))

for b in boundary:
    ID[b-1] = 0

id = 1
for i in range(392):
    if ID[i] != 0:
        ID[i] = id
        id += 1

# LM Matrix

LM = [[ 0 for _ in range(336)] for _ in range(4)]

for a in range(4):
    for e in range(336):
        LM[a][e] = ID[IEN[a][e]-1]

# Element level stiffness matrix
import math

k = [[0 for _ in range(4)] for _ in range(4)]

N = [-1*math.sqrt(3/5) , -1*math.sqrt(3/5), -1*math.sqrt(3/5), 0 , 0 , 0 , math.sqrt(3/5) , math.sqrt(3/5) , math.sqrt(3/5)]
E = [-1*math.sqrt(3/5) , 0 , math.sqrt(3/5) , -1*math.sqrt(3/5) , 0 , math.sqrt(3/5) , -1*math.sqrt(3/5) , 0 , math.sqrt(3/5) ]
weights = [ 25/81 , 40/81 , 25/81 , 40/81 , 64/81 , 40/81 , 25/81 , 40/81 ,  25/81]

def Nderx( a , y ):
    if a == 0:
        return -0.5*(1-y)
    if a == 1:
        return 0.5*(1-y)
    if a == 2:
        return 0.5*(1+y)
    if a == 3:
        return -0.5*(1+y)

def Ndery( a , x ):
    if a == 0:
        return -0.5*(1-x)
    if a == 1:
        return -0.5*(1+x)
    if a == 2:
        return 0.5*(1+x)
    if a == 3:
        return 0.5*(1-x)

for a in range(4):
    for b in range(4):
        k_val = 0
        for i in range(9):
            k_val += (Nderx( a , N[i] ) * Nderx( b , N[i] ) + Ndery( a , E[i] ) * Ndery( b , E[i] ))*weights[i]
        k[a][b] = k_val

# printing the 4*4 element level stiffness matrix
for row in k:
    print(row)

# Global stiffness matrix

K = [[0 for _ in range(280)] for _ in range(280)]

for elem in range(336):
    for a in range(4):
        for b in range(4):
            p = LM[a][elem]
            q = LM[b][elem]
            if (p == 0) or (q == 0 ):
                continue
            else:
                K[p-1][q-1] = K[p-1][q-1] + k[a][b]

# Global Force Vector

BOUND = [ 106 , 107 , 108 , 109 , 110 , 111 , 112 , 113 , 114 , 115 , 126 , 127 , 138 , 139 , 150 , 151 , 162 , 163 , 174 , 175 , 186 , 187 , 198 , 199 , 210 , 211 , 222 , 223 , 224 , 225 , 226 , 227 , 228 , 229 , 230 , 231 ]

F = [0 for _ in range(280)]

for elem in BOUND:
    
    f = [0 for _ in range(4)]
    
    if elem == 106 :
        for j in range(4):
            f[j] = - k[j][2]
    if elem in [107 , 108 , 109 , 110 , 111 , 112 , 113 , 114]:
        for j in range(4):
            f[j] = -k[j][2]-k[j][3]
    if elem == 115 :
        for j in range(4):
            f[j] = -k[j][3]
    if elem in [126 , 138 , 150 , 162 , 174 , 186 , 198 , 210]:
        for j in range(4):
            f[j] = -k[j][1]-k[j][2]
    if elem in [127 , 139 , 151 , 163 , 175 , 187 , 199 , 211]:
        for j in range(4):
            f[j] = -k[j][0]-k[j][3]
    if elem == 222:
        for j in range(4):
            f[j] = -k[j][1]
    if elem in [223 , 224 , 225 , 226 , 227 , 228 , 229 , 230]:
        for j in range(4):
            f[j] = -k[j][0]-k[j][1]
    if elem == 231:
        for j in range(4):
            f[j] = -k[j][0]
    
    for a in range(4):
        p = LM[a][elem-1]
        if p != 0:
            F[p-1] = F[p-1] + f[a]

#Solving for the D vector

import numpy as np

K0 = np.matrix(K)
F0 = np.array(F)

kinv = np.linalg.inv(K0)

d0 = np.matmul( kinv , F0 )

X = [ 0 for _ in range(392)]
Y = [ 0 for _ in range(392)]

X0 = -10
for x in range(147):
    if (X0 == 11) :
        X0 = -10
    X[x] = X0
    X0 += 1

X0 = -10
for x in range(147, 245):
    if (X0 == -3):
        X0 = 4
    if (X0 == 11):
        X0 = -10
    X[x] = X0
    X0 += 1
  
X0 = -10
for x in range(245, 392):
    if (X0 == 11):
        X0 = -10
    X[x] = X0
    X0 += 1

Y0 = -11
for y in range(147):
    if (y%21)==0:
        Y0 += 1
    Y[y] = Y0

for y in range(147, 245):
    if (y-147)%14 == 0:
        Y0 += 1
    Y[y] = Y0
    
for y in range(245, 392):
    if (y-245)%21 == 0:
        Y0 += 1
    Y[y] = Y0

value = [ 0 for _ in range(392)]

for i in range(392):
    eq = ID[i]
    if (eq != 0):
        value[i] = d0.__getitem__((0,eq-1))

B1 = list(range(1, 23)) + [42, 43, 63, 64, 84, 85, 105, 106, 126, 127, 147, 148, 161, 162, 175, 176, 189, 190, 203, 204, 217, 218, 231, 232, 245, 246, 266, 267, 287, 288, 308, 309, 329, 330, 350, 351] + list(range(371, 393))
B2 = [ 133, 134, 135, 136, 137, 138, 139, 140, 141, 154, 155, 168, 169, 182, 183, 196, 197, 210, 211, 224, 225, 238, 239, 252, 253, 254, 255, 256, 257, 258, 259, 260 ]

for b in B2:
    value[b-1] = 1

from mpl_toolkits import mplot3d

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import seaborn as sns 

x = np.array(X)
y = np.array(Y)
z = np.array(value)

df = pd.DataFrame.from_dict(np.array([x,y,z]).T)
df.columns = ['X_value','Y_value','Z_value']
df['Z_value'] = pd.to_numeric(df['Z_value'])

pivotted= df.pivot('Y_value','X_value','Z_value')

sns.heatmap(pivotted,cmap='RdBu')
sns.set(rc={'figure.figsize':(20,20)})