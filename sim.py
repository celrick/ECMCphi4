import fft 
import matplotlib.pyplot as plt
import numpy as np
import sys



lattice = fft.FFT(2,32, mass=0.5)

steps = int(50000)
ones = []
twos = []
threes = []
fours = []
total_actions = []

def update_progress(progress):
            barLength = 50 # Modify this to change the length of the progress bar
            status = ""
            if isinstance(progress, int):
                progress = float(progress)
            if not isinstance(progress, float):
                progress = 0
                status = "error: progress var must be float\r\n"
            if progress < 0:
                progress = 0
                status = "Halt...\r\n"
            if progress >= 1:
                progress = 1
                status = "Done...\r\n"
            block = int(round(barLength*progress))
            text = "\rPercent: [{0}] {1}% {2}".format( "#"*block +\
                                   "-"*(barLength-block), progress*100, status)
            sys.stdout.write(text)
            sys.stdout.flush()
            
def std(data):
    return np.std(data)/(np.sqrt(len(data)))



for i in range(steps):
            
    update_progress(float(i)/float(steps))
    lattice.mcmc_move()

    if i % 10 == 0:
        ones.append(lattice.obs_1())
        twos.append(lattice.obs_2())
        threes.append(lattice.obs_3())
        fours.append(lattice.obs_4())
        total_actions.append(lattice.total_action())
print(len(total_actions))

observables = np.array([ones,twos,threes,fours, total_actions])
print(observables)
np.savetxt('observables.txt',observables.transpose())

"""
#Number of steps to use when figuring out averages
sta= -5
av_1, std_1= np.average(ones[sta:]), std(ones[sta:])
av_2, std_2= np.average(twos[sta:]), std(twos[sta:])
av_3, std_3= np.average(threes[sta:]), std(threes[sta:])
av_4, std_4= np.average(fours[sta:]), std(fours[sta:])
av_A, std_A =np.average(total_actions[sta:]), std(total_actions[sta:])

print("[O]The average action per site is % 5.5f +/ % 5.5f \n" %(av_A,std_A))
print("boogiewoo")

print("[O]The average observable 1 is % 5.5f +/ % 5.5f \n" %(av_1,std_1))
print("[O]The average observable 2 is % 5.5f +/ % 5.5f \n" %(av_2,std_2))
print("[O]The average observable 3 is % 5.5f +/ % 5.5f \n" %(av_3,std_3))
print("[O]The average observable 4 is % 5.5f +/ % 5.5f \n" %(av_4,std_4))

#Plotting of Results
xvalues = np.arange(0,len(twos))

plt.plot(xvalues, threes)
plt.show()


plt.plot(xvalues, total_actions)
plt.show()

print("\n\n\n [+] ECMC")


"""
Acceptance = 100.*lattice.acceptance/(steps*lattice.system.size*lattice.multihit)
print("\n\n")
print("[+] Simulation settings \n L = % d \n d = % d \n m = % 5.4f" \
      %(lattice.size[1], lattice.size[0], lattice.mass))
                                                           
print("[+] Moves \n Number of moves %d \n Acceptance = % 5.2f" %(steps,\
                                                                 Acceptance ))


"""
ones = []
twos = []
total_actions = []


for i in range(steps):
            
    update_progress(float(i)/float(steps))
    lattice.ecmc_move()

    if i % 100 == 0:
        ones.append(lattice.obs_1())
        twos.append(lattice.obs_2())
        total_actions.append(lattice.total_action())

#Number of steps to use when figuring out averages


av_1, std_1= np.average(ones[sta:]), std(ones[sta:])
av_2, std_2= np.average(twos[sta:]), std(twos[sta:])
av_A, std_A =np.average(total_actions[sta:]), std(total_actions[sta:])
Acceptance = 100.*lattice.acceptance/steps
print("\n\n")
print("[+] Simulation settings \n L = % d \n d = % d \n" %(lattice.size[1],\
                                                           lattice.size[0]))
print("[+] Moves \n Number of moves %d \n Acceptance = % 5.2f" %(steps,\
                                                                 Acceptance ))
print("[O]The average action per site is % 5.5f +/ % 5.5f \n" %(av_A,std_A))
print("[O]The average observable 1 is % 5.5f +/ % 5.5f \n" %(av_1,std_1))
print("[O]The average observable 2 is % 5.5f +/ % 5.5f \n" %(av_2,std_2))



plt.plot(xvalues, total_actions)
plt.show()
"""
