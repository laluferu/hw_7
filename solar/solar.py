from pylab import *

tray = genfromtxt("solar.dat",delimiter=",")


a = tray[:,1]

b = tray[:,2]

c = tray[:,3]

d = tray[:,4]


fig = figure((figsize(7,6)))
subplot(231)
plot(a,b)
title('Grafica a vs b ')
xlabel('a' )
ylabel('b' )

subplot(232)
plot(a,c)
title('Grafica a vs c ')
xlabel('a' )
ylabel('c' )

subplot(233)
plot(a,d)
title('Grafica a vs d ' )
xlabel('a' )
ylabel('d' )

subplot(234)
plot(b,c)
title('Grafica b vs c ' )
xlabel('b' )
ylabel('c' )

subplot(235)
plot(b,a)
title('Grafica b vs a ' )
xlabel('b' )
ylabel('a' )

subplot(236)
plot(b,d)
title('Grafica b vs d ' )
xlabel('b' )
ylabel('d' )


savefig("solar.pdf")
