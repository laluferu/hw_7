import matplotlib.pyplot as plt
import numpy as np

tray = np.genfromtxt("poblaciones.dat",delimiter=",")

a = tray[:,0]

b = tray[:,1]

c = tray[:,2]

d = tray[:,3]


fig = plt.figure(figsize = (20,20))
plt.subplot(2,3,1)
plt.scatter(a,b)
plt.xlabel(r'$\alpha$' )
plt.ylabel(r'$\beta$' )

plt.subplot(2,3,2)
plt.scatter(a,c)
plt.xlabel(r'$\alpha$' )
plt.ylabel(r'$\gamma$' )

plt.subplot(2,3,3)
plt.scatter(a,d)
plt.xlabel(r'$\alpha$' )
plt.ylabel(r'$\delta$' )

plt.subplot(2,3,4)
plt.scatter(b,c)
plt.xlabel(r'$\beta$' )
plt.ylabel(r'$\gamma$' )

plt.subplot(2,3,5)
plt.scatter(b,d)
plt.xlabel(r'$\beta$')
plt.ylabel(r'$\delta$')

plt.subplot(2,3,3)
plt.scatter(c,d)
plt.xlabel(r'$\gamma$' )
plt.ylabel(r'$\delta$' )

plt.savefig("poblaciones.pdf",dpi = 400)
