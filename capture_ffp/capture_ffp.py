import math
import numpy
import argparse
import os

from matplotlib.pyplot import *

from amuse.units import units, constants, nbody_system
from amuse.units.quantities import AdaptingVectorQuantity
from amuse.datamodel import Particles, ParticlesSuperset
from amuse.community.huayno.interface import Huayno
from amuse.community.kepler.interface import Kepler
from amuse.io import write_set_to_file
from amuse.ext.orbital_elements import orbital_elements_from_binary, new_binary_from_orbital_elements

def rm_file(file_name):
    """
    delete the file 'file_name', if it exists
    """
    if(os.path.isfile(file_name)):
        os.system('rm %s' % (file_name))

def initialize_code(bodies, code=Huayno, timestep_parameter=0.0169):
    """
    initialize gravity code for bodies
    """

    stars_gravity = code()
    stars_gravity.particles.add_particles(bodies)
    stars_gravity.commit_particles()
    stars_gravity.parameters.timestep_parameter = timestep_parameter  # default 0.0169 time

    return stars_gravity

def get_planets(m0,
                a_planets,
                m_planets,
                e_planets,
                phi=0.):

    # particle set with all planets
    planets = Particles()

    # initialize star--planet orbit
    for i,a_i in enumerate(a_planets):
        star_planet = new_binary_from_orbital_elements(m0,
                                                       m_planets[i],
                                                       a_planets[i],
                                                       e_planets[i],
                                                       true_anomaly=phi)
        # center on the star
        star_planet.position -= star_planet[0].position
        star_planet.velocity -= star_planet[0].velocity
        planets.add_particle(star_planet[1])

    #print " ** planets:"
    #print planets

    return planets

def get_ffp_in_orbit(m0,
                     b_ffp,
                     vinf_ffp,
                     m_ffp,
                     r_inf):
    """
    initialize attributes of the star and the ffp.
    """

    m0_and_ffp_in_orbit = Particles(2)

    # central star
    m0_and_ffp_in_orbit[0].mass = m0
    m0_and_ffp_in_orbit[0].position = (0,0,0) | nbody_system.length
    m0_and_ffp_in_orbit[0].velocity = (0,0,0) | nbody_system.speed

    # free-floating planet

    m0_and_ffp_in_orbit[1].mass = m_ffp
    m0_and_ffp_in_orbit[1].x = -r_inf
    m0_and_ffp_in_orbit[1].y = b_ffp
    m0_and_ffp_in_orbit[1].z = 0. | nbody_system.length
    m0_and_ffp_in_orbit[1].vx = vinf_ffp
    m0_and_ffp_in_orbit[1].vy = 0. | nbody_system.speed
    m0_and_ffp_in_orbit[1].vz = 0. | nbody_system.speed

    # orbital elements
    mass1, mass2, semimajor_axis, eccentricity, true_anomaly, inclination, long_asc_node, arg_per = \
    orbital_elements_from_binary(m0_and_ffp_in_orbit)
    print " ** star + FFP orbital elements: \n \t\t e =", eccentricity, "\n \t\t a =", semimajor_axis

    return m0_and_ffp_in_orbit

def energies_binaries(bodies, indexA, indexB):
    """
    function to calculate energy of a binary (particle set with two particles)
    """

    labels = ['Star','FFP', 'BP']

    particleA, particleB = bodies[indexA], bodies[indexB]

    m_A, m_B = particleA.mass, particleB.mass

    v_A = particleA.velocity.value_in(nbody_system.speed)
    vsquared_A = sum(v_A*v_A) | nbody_system.speed*nbody_system.speed

    v_B = particleB.velocity.value_in(nbody_system.speed)
    vsquared_B = sum(v_B*v_B) | nbody_system.speed*nbody_system.speed

    kinetic_energy = (m_A*vsquared_A + m_B*vsquared_B)/2.0

    r_AB = (particleA.position - particleB.position).value_in(nbody_system.length)
    rmag_AB = math.sqrt(sum(r_AB*r_AB)) | nbody_system.length

    potential_energy = -m_A*m_B/rmag_AB*(1|nbody_system.length**3 * nbody_system.time**(-2) / nbody_system.mass)

    binary_energy = kinetic_energy+potential_energy

    print '\t\t energy for',labels[indexA],'and',labels[indexB],':',binary_energy

    return binary_energy

def evolve_gravity(bodies,
                   fout,
                   converter,
                   t_end,
                   n_steps,
                   t_start=0.|nbody_system.time):

    # positions and velocities centered on the center of mass
    bodies.move_to_center()

    duration = t_end - t_start
    dt = duration / float(n_steps)
    time = t_start

    gravity = initialize_code(bodies, timestep_parameter = dt.value_in(nbody_system.time))

    channel_from_gr_to_framework = gravity.particles.new_channel_to(bodies)

    x = AdaptingVectorQuantity()
    y = AdaptingVectorQuantity()
    times = AdaptingVectorQuantity()

    #Order> 12, 23, 31
    eccentricities = []
    semimajoraxes = []

    Etot_init = gravity.kinetic_energy + gravity.potential_energy
    Etot = Etot_init

    print " ** evolving: \n \t\t t_start = ", converter.to_si(t_start).as_quantity_in(units.yr)
    print "\t\t t_end = ", converter.to_si(t_end).as_quantity_in(units.yr)
    print "\t\t dt = ", converter.to_si(dt).as_quantity_in(units.yr), "\n"
    print "\t\t E_tot = ", converter.to_si(Etot_init), "\n"
    #print " \t\t", "time", "\t\t", "dE/E", "\t\t\t", "[r_m0_FFP, r_m0_planets]"

    while time<=duration:

        gravity.evolve_model(time)
        channel_from_gr_to_framework.copy()

        bodies.collection_attributes.timestamp = time + t_start

        gravity.particles.collection_attributes.timestamp = time + t_start

        x.append(bodies.x)
        y.append(bodies.y)
        times.append(time)

        mass1, mass2, semimajor_axis_12, eccentricity_12, true_anomaly, inclination, long_asc_node, arg_per = orbital_elements_from_binary([bodies[0], bodies[1]])
        mass2, mass3, semimajor_axis_23, eccentricity_23, true_anomaly, inclination, long_asc_node, arg_per = orbital_elements_from_binary([bodies[1], bodies[2]])
        mass3, mass1, semimajor_axis_31, eccentricity_31, true_anomaly, inclination, long_asc_node, arg_per = orbital_elements_from_binary([bodies[2], bodies[0]])

        eccentricities.append([eccentricity_12,eccentricity_23,eccentricity_31])
        semimajoraxes.append([semimajor_axis_12.value_in(nbody_system.length),semimajor_axis_23.value_in(nbody_system.length),semimajor_axis_31.value_in(nbody_system.length)])

        Ekin = gravity.kinetic_energy
        Epot = gravity.potential_energy
        Etot = Ekin + Epot
        dE = Etot_init-Etot

        distances = (bodies.distances_squared(bodies[0])).sqrt()

        #print " \t\t", converter.to_si(time+t_start).as_quantity_in(units.yr), "\t", dE/Etot_init, "\t", converter.to_si(distances[1:].in_(nbody_system.length)).as_quantity_in(units.AU)

        if fout is not None:
            write_set_to_file(bodies, fout, "hdf5")

        time += dt

    print "\t\t (E_f - E_i)/E_tot = ", (Etot-Etot_init)/Etot_init, "\n"

    print " ** binary energies:"

    e12=energies_binaries(bodies, 0, 1).value_in(nbody_system.energy)
    e23=energies_binaries(bodies, 1, 2).value_in(nbody_system.energy)
    e31=energies_binaries(bodies, 2, 0).value_in(nbody_system.energy)

    results = ['flyby', 'exchange', 'temporary capture']

    if (e12>0 and e31<0 and e23>0):
        res = results[0]
    elif (e12<0 and e31>0 and e23>0):
        res = results[1]
    elif (e12<0 and e31<0 and e23>0):
        res = results[2]
    else:
        res = 'something that is not flyby, exchange or temporary capture has happened!'

    print "\n\t\t result =", res, "\n"

    gravity.stop()

    print " ** done"

    return x,y,times,numpy.array(eccentricities),numpy.array(semimajoraxes)

def plot_trayectory(x,y):

    f=figure(figsize=(35,15))

    x1 = x[:,0].value_in(nbody_system.length)
    x2 = x[:,1].value_in(nbody_system.length)
    x3 = x[:,2].value_in(nbody_system.length)

    y1 = y[:,0].value_in(nbody_system.length)
    y2 = y[:,1].value_in(nbody_system.length)
    y3 = y[:,2].value_in(nbody_system.length)

    plot(x1,y1,'y',label='Star')
    plot(x2,y2,'c',label='FFP')
    plot(x3,y3,'m',label='BP',alpha=0.5)

    title('Trajectory FFP (nbody units)')
    xlabel("$x$", fontsize=20)
    ylabel("$y$", fontsize=20)
    legend()

    scatter(x1[0],y1[0],c='black',marker='*')
    scatter(x2[0],y2[0],c='black')
    scatter(x3[0],y3[0],c='black')

    scatter(x1[-1],y1[-1],c='y',marker='*')
    scatter(x2[-1],y2[-1],c='c')
    scatter(x3[-1],y3[-1],c='m')

    axhline(y=0, xmin=-80, xmax=10, c='black', linestyle='--')
    axvline(x=0, ymin=-5, ymax=2, c='black', linestyle='--')

    savefig('trajectory.png')

    xlim(-3,1.1)
    ylim(-1,1)

    savefig('trajectory_zoom.png')

    close()

def plot_orbital_elements(times,eccentricities,semimajoraxes):

    nrows = 2
    ncols = 3

    f = figure(figsize=(35,15))
    labels = ['Star and FFP','FFP and BP', 'BP and Star']

    times = times.value_in(nbody_system.time)

    for i in range(1,4):

        subplot = f.add_subplot(nrows,ncols,i)

        eccs = eccentricities[:,i-1]
        smas = semimajoraxes[:,i-1]

        subplot.plot(times, eccs, c='black', label='$e_{final} = $'+str(eccs[-1]))
        subplot.axhline(y=1, xmin=0, xmax=times[-1], c='m')

        subplot.set_title("Eccentricity for: "+labels[i-1], fontsize=20)
        subplot.set_xlabel("$t$", fontsize=20)
        subplot.set_ylabel("$e$", fontsize=20)
        subplot.legend(loc='best', fontsize=20)
        #subplot.set_ylim(0,3000)

        subplot=f.add_subplot(nrows,ncols,i+3)

        subplot.plot(times, smas, c='black', label='$a_{final} = $'+str(smas[-1]))

        subplot.set_title("Semimajor axis for: "+labels[i-1], fontsize=20)
        subplot.set_xlabel("$t$", fontsize=20)
        subplot.set_ylabel("$a$", fontsize=20)
        subplot.legend(loc='best', fontsize=20)

    savefig('orbitalelements.png')
    close()

def new_option_parser():
    result = argparse.ArgumentParser()
    result.add_argument("--phi",
                      dest="phi", default = 0., type=float, action="store",
                      help="initial angle for of the bounded planet [%default]")
    result.add_argument("--fout",
                      dest="fout",  default="capture_ffp.hdf5", type=str, action="store",
                      help="output file [data.hdf5]")
    result.add_argument("--m0",
                      dest="m0", default = 1.0, type=float, action="store",
                      help="mass of the disk-central star in MSun [%default]")
    result.add_argument("--b",
                      dest="b", default = 5.0, type=float, action="store",
                      help="impact parameter of the FFP in AU [%default]")
    result.add_argument("--vinf",
                      dest="vinf", default = 3.0, type=float, action="store",
                      help="velocity of the FFP at a large distance (infinity) in km/s [%default]")
    result.add_argument("--m_ffp",
                      dest="m_ffp", default = 1.0, type=float, action="store",
                      help="mass of the FFP MJupiter [%default]")
    result.add_argument("--n_steps",
                      dest="n_steps", default = 10, type=int, action="store",
                      help="number of steps for the integrator [%default]")

    return result

if __name__ in ('__main__', '__plot__'):

    arg = new_option_parser().parse_args()

    converter = nbody_system.nbody_to_si(1 | units.MSun,  5 | units.AU)

    # single planet on a circular orbit
    a_planets = [converter.to_nbody(5. | units.AU)]
    m_planets = [converter.to_nbody(1. | units.MJupiter)]
    e_planets = [0.]  # circular orbit
    phi = arg.phi

    # initialize planets
    planets = get_planets(converter.to_nbody(arg.m0 | units.MSun),
                          a_planets,
                          m_planets,
                          e_planets,
                          phi)

    m0_units = converter.to_nbody(arg.m0 | units.MSun)
    b_units = converter.to_nbody(arg.b | units.AU)
    m_ffp_units = converter.to_nbody(arg.m_ffp | units.MJupiter)

    # r_inf -- initial distance to the planet (x-cordinate)
    #r_inf = min(1000., 40.*max(a_planets.value_in(nbody_system.length))) | nbody_system.length
    r_inf = 40.*max(a_planets)

    # set the velocity of FFP assuming parabolic orbit with respect to the star
    if arg.vinf is None:
        mu = m0_units*m_ffp_units / (m0_units + m_ffp_units)
        ep_m0_ffp = ((1.0 | nbody_system.length)**3 * nbody_system.mass**-1 * nbody_system.time**-2)*m0_units*m_ffp_units / ((b_ffp_units**2 + r_inf**2).sqrt())
        vx_ffp = (2.*ep_m0_ffp/mu).sqrt()
        vinf_units = vx_ffp
    else:
        vinf_units = converter.to_nbody(arg.vinf | units.kms)

    # initialize star and FFP
    star_and_ffp_in_orbit = get_ffp_in_orbit(m0_units,
                                             b_units,
                                             vinf_units,
                                             m_ffp_units,
                                             r_inf)

    # particle superset: star, FFP, planets
    bodies = ParticlesSuperset([star_and_ffp_in_orbit, planets])

    #paper = 650 yr
    t_end = converter.to_nbody(650. | units.yr)
    n_steps = arg.n_steps

    if arg.fout is not None:
        rm_file(arg.fout)

    x,y,times,eccentricities,semimajoraxes = evolve_gravity(bodies,
                                                           arg.fout,
                                                           converter,
                                                           t_end,
                                                           n_steps,
                                                           t_start= converter.to_nbody(0. | units.yr) )

    plot_trayectory(x,y)
    plot_orbital_elements(times,eccentricities,semimajoraxes)

##### TO DO:

# * create parameter phi for the planet in orbit -- so i can try with the values of b and phi stated in the paper (what are the units of phi in the plot? radians?)
# * fix the units of the output (i think some of them are still in nbody, not physical ones)


#### DONE:

# * units -- converter for units from the paper <--> physical units
# * plot of the experiment in the XY plane -- save and plot positions (xy) of the experiment and plot the result
# * orbital elements -- save eccentricity and semi-major axis of the star--FFP and star--planet pairs
#                    -- plot evolution of the orbital elements
# * energies -- function to calculate energy of a binary (particle set with two particles)
