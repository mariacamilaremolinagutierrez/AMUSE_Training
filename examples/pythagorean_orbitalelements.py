"""
Calculates the orbital elements of the 3-body problem using dt=0.01
"""

import numpy
import time
import sys

from amuse.community.ph4.interface import ph4

from amuse.units import nbody_system
from amuse.units.quantities import AdaptingVectorQuantity
from amuse.datamodel import Particles
from amuse.ext.orbital_elements import orbital_elements_from_binary
from matplotlib import pyplot


def new_particles():
    particles = Particles(3)

    particles.mass=[3.,4.,5.] | nbody_system.mass
    particles.position = [
        [1, 3, 0],
        [-2, -1, 0],
        [1, -1, 0],
    ] | nbody_system.length

    particles.velocity = [0.,0.,0.] | nbody_system.speed
    particles.radius = 0 | nbody_system.length

    return particles

def run_pyth(te=100):

    dt = 0.01 | nbody_system.time
    t_end = te | nbody_system.time

    code = ph4()

    code.parameters.timestep_parameter = dt.value_in(nbody_system.time)
    code.particles.add_particles(new_particles())

    x = AdaptingVectorQuantity()
    y = AdaptingVectorQuantity()

    t = 0. | nbody_system.time

    eccents_12 = []
    eccents_23 = []
    eccents_31 = []

    semmajaxs_12 = []
    semmajaxs_23 = []
    semmajaxs_31 = []

    times = []

    while(t < t_end-dt/2):

        t=t+dt

        code.evolve_model(t)

        x.append(code.particles.x)
        y.append(code.particles.y)

        mass1, mass2, semimajor_axis_12, eccentricity_12, true_anomaly, inclination, long_asc_node, arg_per = orbital_elements_from_binary([code.particles[0], code.particles[1]])
        mass2, mass3, semimajor_axis_23, eccentricity_23, true_anomaly, inclination, long_asc_node, arg_per = orbital_elements_from_binary([code.particles[1], code.particles[2]])
        mass3, mass1, semimajor_axis_31, eccentricity_31, true_anomaly, inclination, long_asc_node, arg_per = orbital_elements_from_binary([code.particles[2], code.particles[0]])

        eccents_12.append(eccentricity_12)
        eccents_23.append(eccentricity_23)
        eccents_31.append(eccentricity_31)

        semmajaxs_12.append(semimajor_axis_12.value_in(nbody_system.length))
        semmajaxs_23.append(semimajor_axis_23.value_in(nbody_system.length))
        semmajaxs_31.append(semimajor_axis_31.value_in(nbody_system.length))

        times.append(t.value_in(nbody_system.time))

    code.stop()

    return x,y,times,semmajaxs_12,semmajaxs_23,semmajaxs_31,eccents_12,eccents_23,eccents_31


if __name__ in ('__main__', '__plot__'):

    if(len(sys.argv) != 2):
        print 'Please enter the time of integration'
        sys.exit()
    else:
        t_end = int(sys.argv[1])

    x, y, times, sma_12, sma_23, sma_31, ecc_12, ecc_23, ecc_31 = run_pyth(t_end)

    nrows = 2
    ncols = 4

    f=pyplot.figure(figsize=(36,18))

    subplot=f.add_subplot(nrows,ncols,1)

    x1 = x[:,0].value_in(nbody_system.length)
    x2 = x[:,1].value_in(nbody_system.length)
    x3 = x[:,2].value_in(nbody_system.length)

    y1 = y[:,0].value_in(nbody_system.length)
    y2 = y[:,1].value_in(nbody_system.length)
    y3 = y[:,2].value_in(nbody_system.length)

    subplot.plot(x1,y1,'r',label='Particle 1')
    subplot.plot(x2,y2,'b',label='Particle 2')
    subplot.plot(x3,y3,'g',label='Particle 3')

    subplot.set_xlabel("$x$", fontsize=20)
    subplot.set_ylabel("$y$", fontsize=20)
    subplot.set_xlim(-8,8)
    subplot.set_ylim(-6,6)
    subplot.set_title("Trajectory", fontsize=20)
    subplot.legend(loc='best', fontsize=20)

    semimajoraxes = [sma_12, sma_23, sma_31]
    eccentricities = [ecc_12,ecc_23,ecc_31]
    labels = ["1 and 2","2 and 3","3 and 1"]

    for i in range(2,5):

        subplot=f.add_subplot(nrows,ncols,i)

        subplot.plot(times, eccentricities[i-2], c='black', label='$e_{final} = $'+str(eccentricities[i-2][-1]))
        subplot.axhline(y=1, xmin=0, xmax=times[-1], c='m')

        subplot.set_title("Eccentricity for particles "+labels[i-2], fontsize=20)
        subplot.set_xlabel("$t$", fontsize=20)
        subplot.set_ylabel("$e$", fontsize=20)
        subplot.legend(loc='best', fontsize=20)
        #subplot.set_ylim(0,3000)

        subplot=f.add_subplot(nrows,ncols,i+4)

        subplot.plot(times, semimajoraxes[i-2], c='black', label='$a_{final} = $'+str(semimajoraxes[i-2][-1]))

        subplot.set_title("Semimajor axis for particles "+labels[i-2], fontsize=20)
        subplot.set_xlabel("$t$", fontsize=20)
        subplot.set_ylabel("$a$", fontsize=20)
        subplot.legend(loc='best', fontsize=20)

    pyplot.savefig('pythagorean_orbitalelements_tend'+str(t_end)+'.png')
    pyplot.close()
