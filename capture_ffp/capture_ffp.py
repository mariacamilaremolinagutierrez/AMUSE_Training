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

def get_ffp_in_orbit(m0, m, a, e, phi):
    """
    initialize attributes of the star and the ffp.
    """

    star_ffp = new_binary_from_orbital_elements(m0, m, a, e, true_anomaly=phi)
    star_ffp.position -= star_ffp[0].position
    star_ffp.velocity -= star_ffp[0].velocity

    print 'v_inf_orb =',(star_ffp[1].vx**2+star_ffp[1].vy**2+star_ffp[1].vz**2).sqrt()
    print star_ffp

    return star_ffp

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

    #distances = (bodies.distances_squared(bodies[0])).sqrt()
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

    x = AdaptingVectorQuantity()
    y = AdaptingVectorQuantity()
    times = AdaptingVectorQuantity()

    #Order: 12, 23, 31
    eccentricities = []
    semimajoraxes = []

    Etot_init = gravity.kinetic_energy + gravity.potential_energy
    Etot = Etot_init
    DeltaE_max = 0.0 | nbody_system.energy

    times.append(time)
    system_energies = [Etot.value_in(nbody_system.energy)]

    while time<=t_end:

        gravity.evolve_model(time)
        channel_from_gr_to_framework.copy()

        bodies.collection_attributes.timestamp = time

        gravity.particles.collection_attributes.timestamp = time

        x.append(bodies.x)
        y.append(bodies.y)
        times.append(time)

        for i in range(0,number_of_planets+2):
            for j in range(i+1,number_of_planets+2):

                binary = [bodies[i], bodies[j]]
                massA, massB, semimajor_axis, eccentricity, true_anomaly, inclination, long_asc_node, arg_per = orbital_elements_from_binary(binary)

                eccentricities.append(eccentricity)
                semimajoraxes.append(semimajor_axis.value_in(nbody_system.length))

        Etot = gravity.kinetic_energy + gravity.potential_energy
        DeltaE = abs(Etot-Etot_init)

        system_energies.append(Etot.value_in(nbody_system.energy))

        if ( DeltaE > DeltaE_max ):
            DeltaE_max = DeltaE

        time += dt

    print "Energy Change: max(|E_j - E_initial|)/E_initial = ", DeltaE_max/Etot_init

    results = ['flyby', 'temporary capture', 'exchange','nothing']
    #results = [0,1,2,-1]

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
        #res = 'something that is not flyby, exchange or temporary capture has happened!'
        res = results[3]

    print res

    gravity.stop()

    return x,y,times,numpy.array(system_energies),numpy.array(eccentricities),numpy.array(semimajoraxes)

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

    initial_energy = energies[0]
    energies = (energies-initial_energy)/initial_energy
    times = times.value_in(nbody_system.time)

    f=figure(figsize=(15,15))

    plot(times,energies,color='black')

    axhline(y=0, xmin=0, xmax=times[-1], c='m', linestyle='--')

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
        #subplot.set_ylim(0,3000)

        subplot=f.add_subplot(nrows,ncols,i+3)

        subplot.plot(times, smas, c='black', label='$a_{final} = $'+str(smas[-1]))

        subplot.set_title("Semimajor axis for: "+labels[i-1], fontsize=20)
        subplot.set_xlabel("$t$", fontsize=20)
        subplot.set_ylabel("$a$", fontsize=20)
        subplot.legend(loc='best', fontsize=20)

    savefig('orbitalelements.png')
    close()

def convert_units(converter, t_end_p, m0_p, m_ffp_p, m_planets_p, a_planets_p, e_planets_p, n_steps_p, phi_p, b_p):

    #time of integration in yr
    t_end = converter.to_nbody(t_end_p | units.yr)

    #mass of the disk-central star in MSun
    m0 = converter.to_nbody(m0_p | units.MSun)

    #mass of the FFP MJupiter
    m_ffp = converter.to_nbody(m_ffp_p | units.MJupiter)

    #list of the masses of the bounded planets in MJupiter
    m_planets = converter.to_nbody(m_planets_p | units.MJupiter)

    #list of the semimajor axes of the bounded planets in AU
    a_planets = converter.to_nbody(a_planets_p | units.AU)

    #list of the eccentricities of the bounded planets
    e_planets = e_planets_p

    #number of steps for the integrator
    n_steps = n_steps_p

    #initial angle for of the bounded planet in degrees
    phi = phi_p

    #impact parameter of the FFP in AU
    b = converter.to_nbody(b_p | units.AU)

    #Set the velocity of FFP assuming parabolic orbit with respect to the star
    r_inf = 40.*max(a_planets)
    mu = m0*m_ffp / (m0 + m_ffp)
    ep_m0_ffp = (1.0 | nbody_system.length**3 * nbody_system.time**-2 * nbody_system.mass**-1)*m0*m_ffp / ((b**2 + r_inf**2).sqrt())
    vx_ffp = (2.*ep_m0_ffp/mu).sqrt()
    print 'v_inf_theo =',vx_ffp

    return t_end, m0, m_ffp, m_planets, a_planets, e_planets, n_steps, phi, b

def run_capture(t_end_p=650.0, m0_p=1.0, m_ffp_p=1.0, m_planets_p=[1.0], a_planets_p=[5.0], e_planets_p=[0.0], n_steps_p=10000, phi_p=0.0, b_p=5.0):

    """
    Units: t_end(yr), m0(MSun), m_ffp(MJupiter), m_planets(MJupiter), a_planets(AU), e_planets(None), n_steps(None), phi(degrees), b(AU)
    """

    #Time starts
    start_time = time.time()

    #Converter used in this program
    converter = nbody_system.nbody_to_si(1 | units.MSun,  5 | units.AU)

    #Conversion of units
    t_end, m0, m_ffp, m_planets, a_planets, e_planets, n_steps, phi, b = convert_units(converter, t_end_p, m0_p, m_ffp_p, m_planets_p, a_planets_p, e_planets_p, n_steps_p, phi_p, b_p)

    #Initialize planets
    planets = get_planets(m0,a_planets,m_planets,e_planets,phi)

    #Number of planets
    number_of_planets = len(e_planets)
    #Initial distance to the planet (x-cordinate)
    r_inf = 40.*max(a_planets)

    #Initialize star and FFP
    # a_ffp = -r_inf**2
    # b_ffp = -b
    # e_ffp = (b_ffp/a_ffp)**2 + 1
    # phi_ffp = -176.5
    # phi_ffp = 180 + numpy.degrees(math.atan(r_0/r_inf))

    e_ffp = 1.000001
    a_ffp = -(r_inf**2+(b**2)/(e_ffp**2-1)).sqrt()
    phi_ffp = -179.5

    print 'a_ffp =',a_ffp

    #Set the parabolic orbit of the ffp around the star
    star_and_ffp_in_orbit = get_ffp_in_orbit(m0, m_ffp, a_ffp, e_ffp, phi_ffp)

    #Particle superset: star, FFP, planets
    bodies = ParticlesSuperset([star_and_ffp_in_orbit, planets])

    #Evolve time
    x,y,times,energies,eccentricities,semimajoraxes = evolve_gravity(bodies,number_of_planets,converter,t_end,n_steps)

    plot_trajectory(x,y,number_of_planets)
    plot_energy(times, energies)
    #plot_orbital_elements(times,eccentricities,semimajoraxes,number_of_planets)

    print '\nTime:', time.time()-start_time, 'seconds.'

if __name__ in ('__main__', '__plot__'):

    run_capture()
