# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 10:49:06 2023

@author: H. Morita
"""

import numpy as np
import matplotlib.pyplot as plt

#原子数
N = 10

#座標を生成するパラメータ
r = np.linspace(0,2*np.pi,N+1)[:-1]
s = 0.5

#リングの原子座標
xyz = np.c_[np.cos(r), 0.5*np.sin(r), (np.random.rand(N)-0.5)*s]

plt.figure(figsize=(8,4))
plt.subplot(1,2,1)
plt.scatter(xyz[:,0],xyz[:,1])
plt.xlabel('x')
plt.ylabel('y')
plt.grid()
plt.axis('equal')
plt.subplot(1,2,2)
plt.scatter(xyz[:,0],xyz[:,2])
plt.axis('equal')
plt.xlabel('x')
plt.ylabel('z')
plt.grid()
plt.tight_layout()
plt.show()

#ingの中心座標をゼロにする
ring_centers = xyz.mean(axis=0)
xyz0 = xyz - ring_centers

#特異値分解
svd_u, svd_s, svd_uh = np.linalg.svd(xyz0)

#特異値（=分散共分散行列の固有値の平方根）
ring_ellipsoid_lengths = svd_s

#リングの丸さ
roundness = ring_ellipsoid_lengths[1]/ring_ellipsoid_lengths[0]

#リングの厚み
roughness = ring_ellipsoid_lengths[2]/np.sqrt(ring_ellipsoid_lengths[0]*ring_ellipsoid_lengths[1])

print([roundness, roughness])
