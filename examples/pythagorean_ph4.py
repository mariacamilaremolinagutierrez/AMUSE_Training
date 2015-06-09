"""
Calculates the Pythagorean 3-body problem using different values for
the smoothing length in the n-body code.
"""

import numpy
import time

from amuse.community.ph4.interface import ph4

from amuse.units import nbody_system
from amuse.units.quantities import AdaptingVectorQuantity
from amuse.datamodel import Particles
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

def run_pyth(interface,tend=100,dt=0.125, label=""):

    code = interface()

    code.parameters.timestep_parameter = dt.value_in(nbody_system.time)
    code.particles.add_particles(new_particles())

    Ekin = code.kinetic_energy
    Epot = code.potential_energy
    Einit = Ekin + Epot

    Energy = [Einit.value_in(nbody_system.energy)]

    x = AdaptingVectorQuantity()
    y = AdaptingVectorQuantity()

    t = 0. | nbody_system.time

    Times = [t.value_in(nbody_system.time)]

    while(t < tend-dt/2):

        t=t+dt

        code.evolve_model(t)

        x.append(code.particles.x)
        y.append(code.particles.y)

        Ekin = code.kinetic_energy
        Epot = code.potential_energy
        Etot = Ekin + Epot

        Energy.append(Etot.value_in(nbody_system.energy))
        Times.append(t.value_in(nbody_system.time))

    code.stop()

    return x,y,Times,Energy


if __name__ in ('__main__', '__plot__'):

    parameters = [1.0,0.1,0.01,0.001]

    codes_to_run=[ ('ph4: dt = '+str(parameters[0]), ph4,  parameters[0]),
                 ('ph4: dt = '+str(parameters[1]), ph4,  parameters[1]),
                 ('ph4: dt = '+str(parameters[2]), ph4,  parameters[2]),
                 ('ph4: dt = '+str(parameters[3]), ph4,  parameters[3])]

    nrows = 2
    ncols = 4
    f=pyplot.figure(figsize=(32,16))

    for i,(label,interface,dts) in enumerate(codes_to_run):

        x,y,times,energies = run_pyth(interface, tend=100 | nbody_system.time, dt=dts | nbody_system.time, label=label)

        x = x.value_in(nbody_system.length)
        y = y.value_in(nbody_system.length)

        subplot=f.add_subplot(nrows,ncols,i+1)

        subplot.plot(x[:,0],y[:,0],'r')
        subplot.plot(x[:,1],y[:,1],'b')
        subplot.plot(x[:,2],y[:,2],'g')

        subplot.set_title(label)
        subplot.set_xlabel("$x$", fontsize=20)
        subplot.set_ylabel("$y$", fontsize=20)
        subplot.set_xlim(-8,8)
        subplot.set_ylim(-6,6)


        subplot=f.add_subplot(nrows,ncols,i+5)

        subplot.plot(times, (energies-energies[0])/energies[0])

        subplot.set_title("Total Change: "+str((energies[-1]-energies[0])/energies[0])+"\n with dt = "+str(parameters[i]))
        subplot.set_xlabel("$t$", fontsize=20)
        subplot.set_ylabel(r"$\Delta E/E$", fontsize=20)

    pyplot.savefig("pythagorean_ph4.png")
    pyplot.close()
