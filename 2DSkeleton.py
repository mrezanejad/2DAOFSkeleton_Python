#!/usr/bin/env python
# coding: utf-8

# # 2D AOF Skeleton
#This is a jupyter notebook for 2D AOF Skeletonization code

# In[ ]:





# In[58]:


import matplotlib.pyplot as plt
from scipy.misc import imresize


# In[59]:


import matplotlib.image as mpimg
import math
import numpy as np

def rgb2gray(rgb):
    return np.dot(rgb[...,:3], [0.2989, 0.5870, 0.1140])


# In[60]:


fileName = "horse.png"


# In[61]:


I=mpimg.imread(fileName)


# In[62]:


I = rgb2gray(I)
I = I[100:500,100:500]
I = imresize(I,0.3)
print(I.shape)


# In[63]:


imgplot = plt.imshow(I)


# In[64]:


plt.show()


# In[65]:


number_of_samples = 60
epsilon = 1 
flux_threshold = 18


# In[66]:


import scipy.ndimage.morphology as morphOps


# In[67]:


distImage,IDX = morphOps.distance_transform_edt(I,return_indices=True);


# In[68]:


plt.imshow(distImage)


# In[69]:


def sample_sphere_2D(number_of_samples):
    sphere_points = np.zeros((number_of_samples,2))
    alpha = (2*math.pi)/(number_of_samples)
    for i in range(number_of_samples):
        sphere_points[i][0] = math.cos(alpha*(i-1))
        sphere_points[i][1] = math.sin(alpha*(i-1))
    return sphere_points


# In[70]:


print(number_of_samples)


# In[71]:


sphere_points = sample_sphere_2D(number_of_samples)


# In[72]:


def sub2ind(array_shape, rows, cols):
    ind = rows*array_shape[1] + cols
    ind[ind < 0] = -1
    ind[ind >= array_shape[0]*array_shape[1]] = -1
    return ind

def ind2sub(array_shape, ind):
    ind[ind < 0] = -1
    ind[ind >= array_shape[0]*array_shape[1]] = -1
    rows = (ind.astype('int') / array_shape[1])
    cols = ind % array_shape[1]
    return (rows, cols)

def compute_aof(distImage ,IDX,sphere_points,epsilon):

    m = distImage.shape[0]
    n = distImage.shape[1]
    normals = np.zeros(sphere_points.shape)
    fluxImage = np.zeros((m,n))
    for t in range(0,number_of_samples):
        normals[t] = sphere_points[t]
    sphere_points = sphere_points * epsilon
    
    XInds = IDX[0]
    YInds = IDX[1]
    
    for i in range(0,m):
        print(i)
        for j in range(0,n):       
            flux_value = 0
            if (distImage[i][j] > -1.5):
                if( i > epsilon and j > epsilon and i < m - epsilon and j < n - epsilon ):
#                   sum over dot product of normal and the gradient vector field (q-dot)
                    for ind in range (0,number_of_samples):
                                                
#                       a point on the sphere
                        px = i+sphere_points[ind][0]+0.5;
                        py = j+sphere_points[ind][1]+0.5;
                        
                        
                        
                        
#                       the indices of the grid cell that sphere points fall into 
                        cI = math.floor(i+sphere_points[ind][0]+0.5)
                        cJ = math.floor(j+sphere_points[ind][1]+0.5)
                                               

#                       closest point on the boundary to that sphere point

                        bx = XInds[cI][cJ]
                        by = YInds[cI][cJ]
#                       the vector connect them
                        qq = [bx-px,by-py]
                    
                        d = np.linalg.norm(qq)
                        if(d!=0):
                            qq = qq / d
                        else:
                            qq = [0,0]                        
                        flux_value = flux_value + np.dot(qq,normals[ind])
            fluxImage[i][j] = flux_value  
    return fluxImage


# In[73]:


fluxImage = compute_aof(distImage,IDX,sphere_points,epsilon)


# In[74]:


print(fluxImage.shape)


# In[75]:


plt.imshow(fluxImage)


# In[ ]:




