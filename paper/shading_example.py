import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource

# example from http://matplotlib.org/1.4.2/examples/pylab_examples/shading_example.html

# test data
X,Y=np.mgrid[-5:5:0.05,-5:5:0.05]
Z=np.sqrt(X**2+Y**2)+np.sin(X**2+Y**2)
# create light source object.
ls = LightSource(azdeg=0,altdeg=65)
# shade data, creating an rgb array.
#rgb = ls.shade(Z,plt.cm.copper)
rgb = ls.shade(Z,plt.cm.gray)
# plot un-shaded and shaded images.
plt.figure(figsize=(12,5))
plt.subplot(121)
#plt.imshow(Z,cmap=plt.cm.copper)
plt.imshow(Z,cmap=plt.cm.gray)
plt.title('imshow')
plt.xticks([]); plt.yticks([])
plt.subplot(122)
plt.imshow(rgb)
plt.title('imshow with shading')
plt.xticks([]); plt.yticks([])
plt.show()
