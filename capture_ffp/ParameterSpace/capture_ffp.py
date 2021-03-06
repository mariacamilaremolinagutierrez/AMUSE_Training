import math, numpy, argparse, os, time

from matplotlib.pyplot import *

from amuse.units import units, constants, nbody_system
from amuse.units.quantities import AdaptingVectorQuantity
from amuse.datamodel import Particles, ParticlesSuperset
from amuse.community.smalln.interface import SmallN
from amuse.ext.orbital_elements import orbital_elements_from_binary, new_binary_from_orbital_elements

def initialize_code(bodies, code=SmallN, timestep_parameter=0.0169):
    """
    initialize gravity code for bodies
    """

    stars_gravity = code()
    stars_gravity.particles.add_particles(bodies)
    stars_gravity.commit_particles()
    stars_gravity.parameters.timestep_parameter = timestep_parameter  # default 0.0169 time

    return stars_gravity

def get_planets(m0, a_planets, m_planets, e_planets, phi=0.):

    #Particle set with all planets
    planets = Particles()

    #Initialize star-planet orbit
    for i,a_i in enumerate(a_planets):
        star_planet = new_binary_from_orbital_elements(m0,
                                                       m_planets[i],
                                                       a_planets[i],
                                                       e_planets[i],
                                                       true_anomaly=phi)
        #Center on the star
        star_planet.position -= star_planet[0].position
        star_planet.velocity -= star_planet[0].velocity
        planets.add_particle(star_planet[1])

    return planets

def get_ffp_in_orbit(m0, b_ffp, vinf_ffp, m_ffp, r_inf):
    """
    initialize attributes of the star and the ffp.
    """

    m0_and_ffp_in_orbit = Particles(2)

    #Central star
    m0_and_ffp_in_orbit[0].mass = m0
    m0_and_ffp_in_orbit[0].position = (0,0,0) | nbody_system.length
    m0_and_ffp_in_orbit[0].velocity = (0,0,0) | nbody_system.speed

    #Free-floating planet

    m0_and_ffp_in_orbit[1].mass = m_ffp
    m0_and_ffp_in_orbit[1].x = -r_inf
    m0_and_ffp_in_orbit[1].y = b_ffp
    m0_and_ffp_in_orbit[1].z = 0. | nbody_system.length
    m0_and_ffp_in_orbit[1].vx = vinf_ffp
    m0_and_ffp_in_orbit[1].vy = 0. | nbody_system.speed
    m0_and_ffp_in_orbit[1].vz = 0. | nbody_system.speed

    return m0_and_ffp_in_orbit

def energies_binaries(bodies, indexA, indexB):
    """
    function to calculate energy of a binary (particle set with two particles)
    """

    labels = ['Star','FFP', 'BP1', 'BP2', 'BP3', 'BP4'] #... temporary

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

    return binary_energy

def evolve_gravity(bodies, number_of_planets, converter, t_end, n_steps):

    #Positions and velocities centered on the center of mass
    bodies.move_to_center()

    time = 0. | nbody_system.time
    dt = t_end / float(n_steps)

    gravity = initialize_code(bodies, timestep_parameter = dt.value_in(nbody_system.time))
    channel_from_gr_to_framework = gravity.particles.new_channel_to(bodies)

    while time<=t_end:

        gravity.evolve_model(time)
        channel_from_gr_to_framework.copy()

        bodies.collection_attributes.timestamp = time

        gravity.particles.collection_attributes.timestamp = time

        time += dt

    results = [0,1,2,-1]

    planets_star_energies = []

    for i in range(0,number_of_planets):
        planet_star_energy = energies_binaries(bodies, 0, i+2).value_in(nbody_system.energy)
        planets_star_energies.append(planet_star_energy)

    planets_star_energy = sum(planets_star_energies)
    ffp_star_energy = energies_binaries(bodies, 0, 1).value_in(nbody_system.energy)

    if (planets_star_energy<0 and ffp_star_energy>0):
        res = results[0]
    elif (planets_star_energy<0 and ffp_star_energy<0):
        res = results[1]
    elif (planets_star_energy>0 and ffp_star_energy<0):
        res = results[2]
    else:
        res = results[3]

    print res

    gravity.stop()

def plot_trajectory(x,y,number_of_planets):

    colors = ['magenta', 'green', 'DarkOrange', 'red']

    f=figure(figsize=(35,15))

    x_star = x[:,0].value_in(nbody_system.length)
    x_ffp = x[:,1].value_in(nbody_system.length)

    y_star = y[:,0].value_in(nbody_system.length)
    y_ffp = y[:,1].value_in(nbody_system.length)

    plot(x_star,y_star,'y',label='Star')
    scatter(x_star[0],y_star[0],c='black',marker='*')
    scatter(x_star[-1],y_star[-1],c='y',marker='*')

    plot(x_ffp,y_ffp,'c',label='FFP')
    scatter(x_ffp[0],y_ffp[0],c='black')
    scatter(x_ffp[-1],y_ffp[-1],c='c')

    for i in range(0, number_of_planets):

        x_planet = x[:,i+2].value_in(nbody_system.length)
        y_planet = y[:,i+2].value_in(nbody_system.length)

        #color_planet = numpy.random.rand(3,1)
        color_planet = colors[i]

        plot(x_planet,y_planet,color=color_planet,label='BP',alpha=0.5)
        scatter(x_planet[0],y_planet[0],c='black')
        scatter(x_planet[-1],y_planet[-1],color=color_planet)

    axhline(y=0, xmin=-80, xmax=10, c='black', linestyle='--')
    axvline(x=0, ymin=-5, ymax=2, c='black', linestyle='--')

    title('Trajectory FFP (nbody units)')
    xlabel("$x$", fontsize=20)
    ylabel("$y$", fontsize=20)
    legend()

    savefig('trajectory.png')

    xlim(-3.1,1.1)
    ylim(-1,1)

    savefig('trajectory_zoom.png')
    close()

def plot_energy(times, energies):

    f=figure(figsize=(15,15))

    times = times.value_in(nbody_system.time)

    plot(times,energies)

    initial_energy = energies[0]
    energies = (energies-initial_energy)/initial_energy

    #plot(times,energies,color='black')

    #axhline(y=0, xmin=0, xmax=times[-1], c='m', linestyle='--')
    axhline(y=initial_energy, xmin=0, xmax=times[-1], c='m', linestyle='--')

    title('Total Energy of the System (nbody units)')
    xlabel("$t$", fontsize=20)
    ylabel("$|\Delta E| / E$", fontsize=20)

    savefig('energy.png')
    close()

def plot_orbital_elements(times,eccentricities,semimajoraxes,number_of_planets):

    nrows = 2
    ncols = 3

    f = figure(figsize=(35,15))
    labels = ['Star','FFP', 'BP']

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

    #Fixed parameters
    result.add_argument("--t_end",
                        dest="t_end", default = 650.0, type=float, action="store",
                        help="time of integration in years [%default]") #paper 650 yr
    result.add_argument("--m0",
                        dest="m0", default = 1.0, type=float, action="store",
                        help="mass of the disk-central star in MSun [%default]")
    result.add_argument("--m_ffp",
                        dest="m_ffp", default = 1.0, type=float, action="store",
                        help="mass of the FFP MJupiter [%default]")
    result.add_argument("--vinf",
                        dest="vinf", default = 3.0, type=float, action="store",
                        help="velocity of the FFP at a large distance (infinity) in km/s [%default]")
    result.add_argument("--m_planets",
                        dest = "m_planets", default = [1.0], type=float, action="store", nargs='*',
                        help="list of the masses of the bounded planets in MJupiter [%default]")
    result.add_argument("--a_planets",
                        dest = "a_planets", default = [5.0], type=float, action="store", nargs='*',
                        help="list of the semimajor axes of the bounded planets in AU [%default]")
    result.add_argument("--e_planets",
                        dest = "e_planets", default = [0.0], type=float, action="store", nargs='*',
                        help="list of the eccentricities of the bounded planets [%default]") #circular orbit

    #Steps of integrator
    result.add_argument("--n_steps",
                        dest="n_steps", default = 10000, type=int, action="store",
                        help="number of steps for the integrator [%default]")

    #Variable parameters
    result.add_argument("--phi",
                        dest="phi", default = 0., type=float, action="store",
                        help="initial angle for of the bounded planet in degrees[%default]")
    result.add_argument("--b",
                        dest="b", default = 5.0, type=float, action="store",
                        help="impact parameter of the FFP in AU [%default]")

    return result

def get_initial_parameters(converter):

    #Options
    arg = new_option_parser().parse_args()

    #Fixed parameters
    m0 = converter.to_nbody(arg.m0 | units.MSun)
    m_ffp = converter.to_nbody(arg.m_ffp | units.MJupiter)
    t_end = converter.to_nbody(arg.t_end | units.yr)
    n_steps = arg.n_steps

    #Planets masses and semimajor axes
    m_planets = converter.to_nbody( arg.m_planets | units.MJupiter )
    a_planets = converter.to_nbody( arg.a_planets | units.AU )
    e_planets = arg.e_planets
    number_of_planets  = len(m_planets)

    #Variable parameters
    phi = arg.phi
    b = converter.to_nbody(arg.b | units.AU)

    #Set the velocity of FFP assuming parabolic orbit with respect to the star
    if arg.vinf is None:
        mu = m0*m_ffp / (m0 + m_ffp)
        ep_m0_ffp = ((1.0 | nbody_system.length)**3 * nbody_system.mass**-1 * nbody_system.time**-2)*m0*m_ffp / ((b**2 + r_inf**2).sqrt())
        vx_ffp = (2.*ep_m0_ffp/mu).sqrt()
        vinf = vx_ffp
    else:
        vinf = converter.to_nbody(arg.vinf | units.kms)

    return t_end,n_steps,number_of_planets,m0,m_planets,a_planets,e_planets,m_ffp,phi,b,vinf


if __name__ in ('__main__', '__plot__'):

    #Converter used in this program
    converter = nbody_system.nbody_to_si(1 | units.MSun,  5 | units.AU)

    #Time starts
    start_time = time.time()

    #Get initial parameters for the functions above
    t_end,n_steps,number_of_planets,m0,m_planets,a_planets,e_planets,m_ffp,phi,b,vinf = get_initial_parameters(converter)

    #Initialize planets
    planets = get_planets(m0,a_planets,m_planets,e_planets,phi)

    #Initial distance to the planet (x-cordinate)
    r_inf = 40.*max(a_planets)

    #Initialize star and FFP
    star_and_ffp_in_orbit = get_ffp_in_orbit(m0,b,vinf,m_ffp,r_inf)

    #Particle superset: star, FFP, planets
    bodies = ParticlesSuperset([star_and_ffp_in_orbit, planets])

    #x,y,times,energies,eccentricities,semimajoraxes = evolve_gravity(bodies,number_of_planets,converter,t_end,n_steps)
    evolve_gravity(bodies,number_of_planets,converter,t_end,n_steps)
