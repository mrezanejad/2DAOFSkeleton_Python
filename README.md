
## 2D AOF Skeleton
This is a step by step jupyter notebook for 2D AOF Skeletonization code



```python
import matplotlib.pyplot as plt
from scipy.misc import imresize
import matplotlib.image as mpimg
import math
import numpy as np
```

```python
def rgb2gray(rgb):
    return np.dot(rgb[...,:3], [0.2989, 0.5870, 0.1140])
```


```python
fileName = "horse.png"
```


```python
I=mpimg.imread(fileName)
```


```python
I = rgb2gray(I)
I = 1-I
I = imresize(I,0.3)
print(I.shape)
```

    (89, 90)


```python
imgplot = plt.imshow(I)
```


![png](output_8_0.png)



```python
plt.show()
```


```python
number_of_samples = 60
epsilon = 1 
flux_threshold = 18

```


```python
import scipy.ndimage.morphology as morphOps
```


```python
distImage,IDX = morphOps.distance_transform_edt(I,return_indices=True);
```


```python
plt.imshow(distImage)
```




    <matplotlib.image.AxesImage at 0xb24fc0d30>




![png](output_13_1.png)



```python
def sample_sphere_2D(number_of_samples):
    sphere_points = np.zeros((number_of_samples,2))
    alpha = (2*math.pi)/(number_of_samples)
    for i in range(number_of_samples):
        sphere_points[i][0] = math.cos(alpha*(i-1))
        sphere_points[i][1] = math.sin(alpha*(i-1))
    return sphere_points
```


```python
print(number_of_samples)
```

    60



```python
sphere_points = sample_sphere_2D(number_of_samples)
```


```python
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
```


```python
fluxImage = compute_aof(distImage,IDX,sphere_points,epsilon)
```

```python
print(fluxImage.shape)
```

    (89, 90)



```python
plt.imshow(fluxImage)
```




    <matplotlib.image.AxesImage at 0x11dda9438>




![png](output_20_1.png)



```python
skeletonImage = fluxImage
skeletonImage[fluxImage < flux_threshold] = 0
skeletonImage[fluxImage > flux_threshold] = 1

```


```python
plt.imshow(skeletonImage)
```




    <matplotlib.image.AxesImage at 0x11de47b00>




![png](output_22_1.png)



```python

```
