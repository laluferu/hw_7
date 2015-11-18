rom pylab import *

tray = genfromtxt("poblaciones.dat",delimiter=",")


a = tray[:,1]

b = tray[:,2]

c = tray[:,3]

d = tray[:,4]


fig = figure((figsize(10,10)))
subplot(231)
plot(a,b)

xlabel(r'$\alpha$' )
ylabel(r'$\beta$' )

subplot(232)
plot(a,c)

xlabel(r'$\alpha$' )
ylabel(r'$\gamma$' )

subplot(233)
plot(a,d)

xlabel(r'$\alpha$' )
ylabel(r'$\delta$' )

subplot(234)
plot(b,c)

xlabel('b' )
ylabel(r'$\gamma$' )

subplot(235)
plot(b,a)

xlabel(r'$\beta$' )
ylabel(r'$\alpha$')

subplot(236)
plot(b,d)

xlabel(r'$\beta$' )
ylabel(r'$\delta$' )


savefig("poblaciones.pdf")
