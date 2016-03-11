from __future__ import division
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib import pyplot
import numpy


###variable declarations


def euler(t):

###Assign initial conditions

    nr = 101
    nt = 1000
    dr = 2/(nr-1)
    rb=1
    dt = 0.1*dr
    g=1.4
    max_time=dt*nt
    print max_time, dr, dt
   

    r = numpy.linspace(0,10,nr)
    l = numpy.empty(nr)  ##rho  create a 1xn vector of 1's
    l.fill(1)
    j = numpy.empty(nr)
    j.fill(1) ##assumindo que as moleculas de ar nao se mexem inicialmente
    e = numpy.empty(nr)
    e.fill(1)
    p = numpy.empty(nr)
    p.fill(1)
    ln = numpy.empty(nr) ##
    jn = numpy.empty(nr)
    en= numpy.empty(nr)
    pn= numpy.empty(nr)
    
    #print r


    #l[0 : .1/dr]=0.2
    #j[.5/dr : .6/dr]=3
    #e[.5/dr : .6/dr]=4
    p[0]=10
    #p[1]=10

    for n in range(int(t)):  #iterate through time
        print "MUDEI DE NNNNNNNNNNNNNNNNNN"
        ln = l.copy() ##copy the existing values of u into un
        jn =j.copy()
        en =e.copy()
        pn = p.copy()
        for i in range(1,nr):
            print "r", i, r[i] , ln[i] , jn[i], en[i], pn[i]
            l[i] = ln[i] - dt/dr*(jn[i]-jn[i-1])-(2*dt*jn[i])/r[i]
            print "rho", ln[i], l[i]
            j[i] = jn[i] - dt/dr*(pn[i]-pn[i-1])-(2*dt*(jn[i])**2)/(ln[i]*r[i])
            print "j", jn[i], j[i]
            e[i] = en[i] - dt/dr*( (jn[i]*en[i])/ln[i]-(jn[i-1]*en[i-1])/ln[i-1]+  (jn[i]*pn[i])/ln[i]-(jn[i-1]*pn[i-1])/ln[i-1] )-(2*dt/r[i])*(jn[i]/ln[i])*(en[i]+pn[i])
            print "e", en[i], e[i]
            p[i]=(en[i]-(jn[i]*jn[i])/ln[i])*g
            print "p", pn[i], p[i]
    print p
    #print jn
    pyplot.figure(1)
    pyplot.xlabel('r')
    #pyplot.ylabel('Probability')
    pyplot.title('Density')
    pyplot.plot(r, l);
    #pyplot.xlimt([0,10])
    pyplot.figure(2)
    pyplot.xlabel('r')
    pyplot.title('Momentum Density')
    pyplot.plot(r, j);
    #pyplot.xlimt([0,10])
    pyplot.figure(3)
    pyplot.title('Energy Density')
    pyplot.xlabel('r')
    pyplot.plot(r, e);
    #pyplot.xlimt([0,10])
    pyplot.figure(4)
    pyplot.title('Pressure')
    pyplot.xlabel('r')
    pyplot.plot(r, p);
    #pyplot.xlimt([0,10])
    pyplot.show()
    

    


euler(1) #IC
#euler(10)
#burger(50)
#burger(100)
#burger(200)
#burger(400)
#burger(700)
#burger(2000)


