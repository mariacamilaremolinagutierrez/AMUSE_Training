import os
import numpy as np

#Parameter Space

# -4r_0 < b < 4r_0
#bs = np.linspace(-4*5,4*5,10)

# 0 < phi < 0.7 (where 1=360 degrees)
#phis = np.linspace(0,360*0.7,10)

bs = [1.0]
phis = [360*0.1]

m_planets = '1.0 2.0'
a_planets = '5.0 8.0'
e_planets = '0.0 0.0'

n_steps = 4000

for b in bs:
    for phi in phis:
        long_command = 'amuse capture_ffp.py --b '+str(b)+' --phi '+str(phi)+' --m_planets '+m_planets+' --a_planets '+a_planets+' --e_planets '+e_planets+' --n_steps '+str(n_steps)+' > results.txt'
        #command = 'amuse capture_ffp.py --b '+str(b)+' --phi '+str(phi)+' > results.txt'
        #os.system(command)
        os.system(long_command)
