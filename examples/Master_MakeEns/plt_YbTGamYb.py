#!/scratch/singadm/pkg_hera/Python/2.7.10/bin/python

import numpy as np
import matplotlib.pyplot as plt

f = open('YbTGamYb.txt',  'r')

z = []
for i in range(30):
  z.append(f.readline().strip().split())

z = np.array(z)
z = z.astype(np.float) 
f.close()

plt.imshow(z, cmap='jet',vmin=-0.05,vmax=1)
#plt.contourf(z, cmpa='jet') 
plt.savefig('YbTGamYb')
plt.show()
