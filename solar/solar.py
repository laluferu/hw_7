import matplotlib.pyplot as plt
import numpy as np

tray = np.genfromtxt("solar.dat",delimiter=",")


a = tray[:,0]

b = tray[:,1]

c = tray[:,2]

d = tray[:,3]


fig = plt.figure(figsize = (20,20))
plt.subplot(2,3,1)
plt.scatter(a,b)
plt.title('Grafica a vs b ')
plt.xlabel('a' )
plt.ylabel('b' )

plt.subplot(2,3,2)
plt.scatter(a,c)
plt.title('Grafica a vs c ')
plt.xlabel('a' )
plt.ylabel('c' )

plt.subplot(2,3,3)
plt.scatter(a,d)
plt.title('Grafica a vs d ' )
plt.xlabel('a' )
plt.ylabel('d' )

plt.subplot(2,3,4)
plt.scatter(b,c)
plt.title('Grafica b vs c ' )
plt.xlabel('b' )
plt.ylabel('c' )

plt.subplot(2,3,5)
plt.scatter(b,d)
plt.title('Grafica b vs d ' )
plt.xlabel('b' )
plt.ylabel('d' )

plt.subplot(2,3,6)
plt.scatter(c,d)
plt.title('Grafica c vs d ' )
plt.xlabel('c' )
plt.ylabel('d' )

plt.savefig("solar.pdf",dpi = 400)
