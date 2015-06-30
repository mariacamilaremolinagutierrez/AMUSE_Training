import os
import numpy as np
import matplotlib.pyplot as plt

def create_file(bs,phis):

    i=1

    par = open('./parameter_space_status.txt','r')
    num_runs = sum(1 for line in par)
    par.close()

    for b in bs:
        for phi in phis:
            if (i > num_runs):
                status_command = 'echo "'+str(i)+') phi = '+str(phi)+' b = '+str(b)+'" >> parameter_space_status.txt'
                os.system(status_command)
                command = 'amuse capture_ffp.py --b '+str(b)+' --phi '+str(phi)+' >> grid.txt'
                os.system(command)
            i += 1

def plot_grid(bs,phis):
    #File with results
    results = np.loadtxt('./grid.txt')

    #Plot
    fig = plt.figure(figsize=(10,6))

    i=0

    for b in bs:
        for phi in phis:

            try:
                r = results[i]
            except:
                break

            if(r==0):
                color = 'y'
            elif(r==1):
                color = 'c'
            elif(r==2):
                color = 'k'
            else:
                color = 'r'

            plt.scatter(phi,b,marker='s',s=100,color=color)

            i=i+1

    plt.xlabel('$\phi$')
    plt.ylabel('$b$')
    plt.xlim(0,360*0.7)
    plt.ylim(-4*5,4*5)
    plt.savefig('parameter_space.png')
    plt.close()

if __name__ in ('__main__', '__plot__'):

    # -4r_0 < b < 4r_0
    bs = np.linspace(-4*5,4*5,100)
    # 0 < phi < 0.7 (where 1=360 degrees)
    phis = np.linspace(0,360*0.7,100)

    create_file(bs,phis)
    plot_grid(bs,phis)
