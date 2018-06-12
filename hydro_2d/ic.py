import numpy as np
from PIL import Image

from matplotlib import pyplot as plt

im = Image.open('alex.jpg')
a = np.array(list(im.getdata()))
a = a.reshape(im.height, im.width, 3)

brightness = np.sum(a, axis = 2)

print('nx = ', im.width)
print('ny = ', im.height)

plt.imshow(brightness)

u_prim = np.zeros((im.height, im.width, 7))

u_prim[:,:,0] = 1.0
u_prim[:,:,1] = a[:,:,0] / a[:,:,0].max() - 0.5
u_prim[:,:,2] = a[:,:,1] / a[:,:,1].max() - 0.5
u_prim[:,:,3] = 0.01 * brightness / brightness.max()
u_prim[:,:,4] = 1.0
u_prim[:,:,5] = 0.0
u_prim[:,:,6] = 0.0

np.savetxt('ic.dat', u_prim.flatten())
