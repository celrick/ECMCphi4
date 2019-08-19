"""
FREE FIELD THEORY in n dimensions.

This is the main class for ecmc simulation of FFT, containing all the 
used methods. 

#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.


Created: 27/05/10
Created by Conor Elrick
"""
import numpy as np 
import matplotlib.pyplot as plt

class FFT():
    
    def __init__(self, n, L, start_type = 'hot', mass = 0.1, multihit = 3):
        self.size = (n,L)
        self.log_file = open("logfile_fft.out", 'w')
        self.system = self.system_setup(start_type)
        self.action = None
        self.mass = mass
        self.sigma = np.random.choice([-1,1]) # direction of ecmc movement
        self.chain_length = 125.0 #event chain length
        self.lifted_site = 0
        self.acceptance = 0.
        self.multihit = multihit
    
    
    def system_setup(self, start_type):
        """
        Initalises n dimensional array of points.
        """
        self.log_file.write("[I] Created an " + str(self.size[0]) + \
                            "-dimensional lattice of size " + \
                            str(self.size[1]) + ". \n") 
        if self.size[0] >= 4 or self.size[0] <1:
            print('[I] Error, dimension of lattice invalid')
            
        elif self.size[0] == 1:
            lattice = np.ndarray((1, self.size[1]))
            
        elif self.size[0] == 2:
            lattice = np.ndarray((self.size[1], self.size[1]))
            
        elif self.size[0] == 3:
            lattice = np.ndarray((self.size[1], self.size[1], self.size[1]))
            
        if start_type == 'hot':
            for i in range(lattice.size):
                lattice.flat[i] = self.rng_gauss()
                
        elif start_type == 'cold':
            for i in range(lattice.size):
                lattice.flat[i] = 0.
        self.log_file.write("[I] Initial setup has average of " + \
                           str(np.mean(lattice.flat)) + \
                           ", with a standard deviation of " + \
                           str(np.std(lattice.flat)) + ". \n" + \
                           "[I] Initial lattice is: \n ")
                           
        for data_slice in lattice:

            # The formatting string indicates that I'm writing out
            # the values in left-justified columns 7 characters in width
            # with 2 decimal places.  
            np.savetxt(self.log_file, data_slice, fmt='%-7.2f')

            # Writing out a break to indicate different slices...
            self.log_file.write('# New slice\n')
        self.log_file.write('\n\n\n')   
        return lattice
                
    def rng_gauss(self, centered = 0.):
        """
        Gaussian RNG
        """
        return np.random.normal(loc=centered,scale=2.0)
        #return np.random.randint(-5,5) #'For testing as a toy model
        
    def FT_lattice(self):
        "Returns the fourier transform of the lattice"
        return np.fft.fft(self.system)
        
    def nearest_neighbours(self, site):
        """
        Finds the nearest neighbours in a FLATTENED ARRAY of all sites. 
        Generalises to n dimensions
        """
        n_points = self.system.size
        
        neighbours = []
        L = self.size[1]
        if self.size[0] == 1:
            neighbours.append((site + 1) % n_points)
            neighbours.append((site - 1) % n_points)
        if self.size[0] == 2:
            tester = np.arange(0,n_points)
            tester = np.reshape(tester,(L,L))
            
            point = (site % L, int(site / L))
            
            neighbours.append(tester[(point[1]+1) % L ][point[0]])
            neighbours.append(tester[(point[1]-1) % L][point[0]])
            neighbours.append(tester[point[1]][(point[0]+1)% L])
            neighbours.append(tester[point[1]][(point[0]-1) % L])
        if self.size[0] >= 3:
            tester = np.arange(0,n_points)
            tester = np.reshape(tester,(L,L,L))
            
            point = (site % L, int(site / L) % L, int(site/(L**2)))
            
            neighbours.append(tester[point[2]][(point[1]+1) % L ]\
                              [point[0]])
            neighbours.append(tester[point[2]][(point[1]-1) % L]\
                              [point[0]])
            neighbours.append(tester[point[2]][point[1]]\
                              [(point[0]+1)% L])
            neighbours.append(tester[point[2]][point[1]]\
                              [(point[0]-1) % L])
            neighbours.append(tester[(point[2]+1) % L][point[1]]\
                              [point[0]])
            neighbours.append(tester[(point[2]-1) % L][point[1]]\
                              [point[0]])
            
        #return list(set(neighbours))
        return neighbours
    
    def action_self(self, site_value):
        """
        Self-interaction term in action. For FFT this is a quadratic.
        """
        return 0.5*(self.mass**2)*(site_value**2)
    
    def action_interaction(self, site_value, neighbour_value):
        """
        Interaction term in the action.
        """
        return  0.5*((site_value - neighbour_value)**2)
    
    def obs_1(self):
        """
        1 of 4 obervables. This is related to a term in the action.
        """
        total = 0.
        for site in range(0, self.system.size):
            differences = [(self.system.flat[site] -\
                           self.system.flat[neigh])**2 for neigh in\
                           self.nearest_neighbours(site)]
            total += np.sum(differences)
        total *= 1.0/(2*self.size[0]*self.system.size)
        return total
    
    def obs_2(self):
        """
        2 of 4 overvables. This is related to a term in the action. An extra
        sum sign is assumed.
        """
        total = 0.
        for site in range(0, self.system.size):
            total += self.system.flat[site]**2
        total *= 1.0/(2*self.system.size)
        return total
        
    def obs_3(self):
        """
        3 of 4 overvables. This is related to magnetisation
        """
        total = 0.
        for site in range(0, self.system.size):
            total += self.system.flat[site]
        total = total **2
        total /= (2.0 * self.system.size**2)
        return total
        
    def obs_4(self):
        """
        4 of 4 overvables. This is related to fourier transform. This is the 
        only one which is implemented in a certain number of dimensions for 
        checking.
        """
        total = 0.
        
        def ft(ind, phi):
            return np.exp(-2j*np.pi*ind / self.size[1]) * phi
         
        if self.size[0] == 1:
            phi_i = 0. + 0j
            for nm in range(self.size[1]):
                phi_i += ft(0, self.system[0][nm])
                #print(phi_i)
            total += phi_i * np.conjugate(phi_i) 
            
        elif self.size[0] == 2:
            phi_i_0 = 0. + 0j
            phi_i_1 = 0. + 0j
            
            #0th direction loop
            for ind in range(self.size[1]):
                for pos in range(self.size[1]):
                    
                    #print(ft(ind, self.system[ind][pos]))
                    phi_i_0 += ft(ind,self.system[ind][pos])
            #print(phi_i_0)
            #1st direction loop
            for ind in range(self.size[1]):
                for pos in range(self.size[1]):
                    phi_i_1 += ft(ind, self.system[pos][ind])
            #print(phi_i_1)
            total += phi_i_0 * np.conjugate(phi_i_0)\
                   + phi_i_1 * np.conjugate(phi_i_1)
            
        elif self.size[0] == 3:
            phi_i_0 = 0. + 0j
            phi_i_1 = 0. + 0j
            phi_i_2 = 0. + 0j
            
            #0th direction loop
            for ind in range(self.size[1]):
                for pos_y in range(self.size[1]):
                    for pos_z in range(self.size[1]):
                    
                    #print(ft(ind, self.system[ind][pos]))
                        phi_i_0 += ft(ind,self.system[ind][pos_y][pos_z])
                    
            #print(phi_i_0)
            #1st direction loop
            for ind in range(self.size[1]):
                for pos_x in range(self.size[1]):
                    for pos_z in range(self.size[1]):
                        phi_i_1 += ft(ind, self.system[pos_x][ind][pos_z])
            #2nd direction loop
            for ind in range(self.size[1]):
                for pos_x in range(self.size[1]):
                    for pos_y in range(self.size[1]):
                        phi_i_2 += ft(ind, self.system[pos_x][pos_y][ind])     
                        
            #print(phi_i_1)
            #print(phi_i_2)
            total += phi_i_0 * np.conjugate(phi_i_0)\
                   + phi_i_1 * np.conjugate(phi_i_1)\
                   + phi_i_2 * np.conjugate(phi_i_2)
            
        #print("h")
        total /= self.system.size
        total = np.real(total)
        return total
        
    def mc_local_step(self,site):
        """
        MCMC move given a site and an action. The neighbours and self are
        split in two cases.
        """
        
        old_site_value = self.system.flat[site]
        new_site_value = old_site_value + np.random.normal(loc=0.,scale=1.7)
        
        #True until proven otherwise
        move_acceptance = True
        
        old_action = 0.
        new_action = 0.
        
        #nearest neighbour term
        for neighbour in self.nearest_neighbours(site):
            
            old_action += self.action_interaction(old_site_value,\
                                                 self.system.flat[neighbour])
            
            new_action += self.action_interaction(new_site_value,\
                                                 self.system.flat[neighbour])
            
        #self term
        old_action += self.action_self(old_site_value)
        new_action += self.action_self(new_site_value)
        
        #Metropolis Step
        P_acceptance = min([1, np.exp(old_action - new_action)])
        gamma = np.random.uniform(0,1)
        
        if P_acceptance <= gamma:
            move_acceptance = False
        
        #Acceptance
        if move_acceptance == True:
            self.acceptance += 1.
            self.system.flat[site] = new_site_value
            #self.log_file.write("Move Accepted")
            #print(self.system)
        """
        print statements for checking
        """
        #print("site to change is " + str(site))
        #print('old sv = ' + str(old_site_value))
        #print('new sv = ' + str(new_site_value))
        #print("neighbour is " + str(neighbour))
        #print("acceptance = " + str(acceptance))
        #print("old-new action is" + str(old_action-new_action))
        #print("gamma is " + str(gamma))
        #print("is move accepted?" +str(move_acceptance))
        return None
        
    def mcmc_move(self):
        """
        local mcmc move of a site
        """
        for site in range(self.system.size):
            for hit in range(self.multihit):
        #site = np.random.randint(0, self.system.size)
                self.mc_local_step(site)
        return None
    
    def ecmc_t_calc(self, energy, delta_phi, mass):
        """
        Method for calculaing the t for a given change in energy. Formulae in
        notebook pg 2.8. Phi is either lattice site value or the delta given.
        """
        
        if self.sigma * delta_phi >= 0:
            #print("TRIGGER")
            return np.sqrt((2.0 * energy / (mass**2)) + delta_phi**2 ) -\
                   self.sigma*delta_phi
        else:
            return np.sqrt(2.0 * energy / (mass**2)) -\
                   self.sigma*delta_phi
    
    def ecmc_energy(self):
        """
        Generates the random change in energy. Note this needs to be called
        for each pair
        """
        return -1 * np.log(np.random.uniform(0,1))
     
    def ecmc_move(self):
        """
        Main loop for a chain of events
        """
        move_var = True
        chain_length = self.chain_length
        self.lifted_site = np.random.randint(0, self.system.size)
        while (move_var == True):
            
            #print("sigma is " + str(self.sigma))
            #print("lifted_site is " + str(self.lifted_site))
            #print("lifted_site value is " + str(self.system.flat[self.lifted_site]))
            neighbours = self.nearest_neighbours(self.lifted_site)
            #print("neighbours are " + str(neighbours))
            
            moves = np.zeros([len(neighbours) + 1,2]) # +1 for self interaction
            #print("moves are " + str(moves))
            
            #Self interaction
            delta_E = self.ecmc_energy()
            #print("self delta_E = " + str(delta_E))
            phi = self.system.flat[self.lifted_site]
            moves[0] = np.array([self.lifted_site, self.ecmc_t_calc(delta_E,\
                                phi, self.mass)])
            #print("moves are " + str(moves))
            
            #Pair interaction
            i = 1
            for neighbour in neighbours:
                delta_phi = self.system.flat[self.lifted_site] -\
                            self.system.flat[neighbour]
                #phi = self.system.flat[self.lifted_site]
                delta_E = self.ecmc_energy()
                #print("Neighbour is " + str(neighbour))
                #print("self delta_E = " + str(delta_E))
                moves[i] = np.array([neighbour, self.ecmc_t_calc(delta_E,\
                                    delta_phi, 1.0)])
                #print("move is " + str(moves[i]))                   
                i += 1
            t_move = moves[np.argmin(np.transpose(moves)[1], axis=0)]
            #print("moves are " + str(moves))
            #print("t_move " + str(t_move))
            
            if t_move[1] < self.chain_length:
                
                
                    
                self.chain_length -= t_move[1]
                
                self.system.flat[self.lifted_site] += self.sigma * t_move[1]
                    
                if t_move[0] == self.lifted_site:    
                
                    self.sigma *= -1
                #print("new value at site is " + str(self.system.flat[self.lifted_site]))
                self.lifted_site = int(t_move[0])
                #print("Chain length is " +str(self.chain_length))
                #print(self.system)
            else:
                self.system.flat[self.lifted_site] += self.sigma *\
                                                (self.chain_length)
                if t_move[0] == self.lifted_site:
                    self.sigma *= -1
                #print("lifted_site value is " + str(self.system.flat[self.lifted_site]))
                move_var = False
                self.chain_length = chain_length
                #print("WOOOOOO" + str(self.chain_length))
           
        return None
    
    def total_action(self):
        """
        Method for calculating the total action of the system. Should be 
        1/2 kbt per degree of freedom.
        """
        action = 0.
        for site_I in range(0,self.system.size):
            n_neighbours = self.nearest_neighbours(site_I)
            #print(site_I)
            #print("N")
            for neigh in n_neighbours:
                #0.5 for double counting
                
                action += 0.5*self.action_interaction(\
                             self.system.flat[site_I], self.system.flat[neigh])
                #print("interaction is "+ str(action))
            action += self.action_self(self.system.flat[site_I])
           
        action = action/(self.system.size)
        return action
        
        
        
        
if __name__ == '__main__':
    obsers1 = []
    obsers2 = []
    obsers3 = []
    obsers4 = []
    action = []
    sizes = [8]
    print("sizes are:")
    print(sizes) 
    
    for j in sizes:
        
        import time, sys

    # update_progress() : Displays or updates a console progress bar
    ## Accepts a float between 0 and 1. Any int will be converted to a float.
    ## A value under 0 represents a 'halt'.
    ## A value at 1 or bigger represents 100%
    ### Stolen from stack exchange
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
        
        
        
        
        lattice = FFT(2, j)
        print(lattice.system)

        total_action = []
        
        ones = []
        twos = []
        threes = [] 
        fours = []
        
        steps = int(50000)
        #print(steps)
        
        #print(lattice.system.size)
        for i in range(steps):
            
            #update_progress(float(i)/float(steps))
            #print("\n")
            
            lattice.mcmc_move()
            #print("\n")
            #print(lattice.system)
            if i% 10000 == 0:
                print(lattice.system)
            if i % 10 == 0:
                ones.append(lattice.obs_1())
                twos.append(lattice.obs_2())
                threes.append(lattice.obs_3())
                fours.append(lattice.obs_4())
                total_action.append(lattice.total_action())
                #print(str(lattice.total_action())+"\n")
        #print(lattice.system)
        xvalues = np.arange(0,len(twos))
        print("Acceptance is " + str(100.*lattice.acceptance/steps) +"%\n")
        #print('ones ' + str(np.average(ones[5:])))
        obsers1.append(np.average(ones[-200:]))
        obsers2.append(np.average(twos[-200:]))
        obsers3.append(np.average(threes[-200:]))
        obsers4.append(np.average(fours[-200:]))
        action.append(np.average(total_action[-1000:]))
        #print(obsers3)        #print(fours)
        #print('threes ' + str(np.average(threes[5:])))
        #print('fours' + str(np.average(fours[5:])))
        #print(lattice.system)
        plt.plot(xvalues[0:], twos[0:], label="Observable_2")
        plt.savefig("Observable_2 for L =" + str(j) + ".png")
        plt.clf()
        plt.plot(xvalues[0:], ones[0:], label="Observable_1")
        plt.savefig("Observable_1 for L =" + str(j) + ".png")
        plt.clf()
        plt.plot(xvalues[0:], threes[0:], label="Observable_3")
        plt.savefig("Observable_3 for L =" + str(j) + ".png")
        plt.clf()
        plt.plot(xvalues[0:], fours[0:], label="Observable_4")
        plt.savefig("Observable_4 for L =" + str(j) + ".png")
        plt.clf()
        plt.plot(xvalues[0:], total_action[0:], label="Total action")
        plt.savefig("total_action for L =" + str(j) + ".png")
        plt.clf()
        #plt.plot(xvalues,threes, label="Observable_3")
       
        #plt.legend()
        #plt.title("Observables 1 and 2 taken every 10 moves")
        #plt.savefig("Observables12.png")
        #plt.show()
        #plt.plot(xvalues,threes)
        #plt.savefig("Observable3.png")
        #plt.show()
        #plt.plot(xvalues,fours)
        #plt.show()
       # 
        #with open("observables.out",'w') as f:
            #for item in xvalues:
                #f.write("1- %s\n" % ones[item])
                #f.write("2- %s\n" % twos[item])
                #f.write("3- %s\n" % threes[item])

        #for data_slice in lattice.system:

            # The formatting string indicates that I'm writing out
            # the values in left-justified columns 7 characters in width
            # with 2 decimal places.  
            #np.savetxt(lattice.log_file, data_slice, fmt='%-7.2f')
            # Writing out a break to indicate different slices...
            #lattice.log_file.write('# New slice\n')
        print("Finshed " + str(j))
    
    print("Observable_1 ")
    print(obsers1)
    print("Observable_2 ")
    print(obsers2)
    print("Observable_3 ")
    print(obsers3)
    print("Observable_4 ")
    print(obsers4)
    print("Action ")
    print(action)
    
    #plt.plot(sizes, obsers3, color = 'black',label='simulation')
    #plt.plot(sizes, 9.0/(4*(np.array(sizes))), color = 'red',label='9/4')
    #plt.plot(sizes, 9.0/((np.square(sizes))), color = 'yellow')
    #plt.plot(sizes, 9.0/((np.array(sizes))), color='blue')
    #plt.plot(sizes, 4.0/((9.0*np.array(sizes))), color = 'blue', linestyle='--',label='4/9')
    #plt.plot(sizes, (9.0/8.0)*np.square(np.array(sizes)))
    #plt.legend()
    #plt.show()
    
    
