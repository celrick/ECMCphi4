"""
Phi^4 THEORY in n dimensions.

This is the main class for ecmc simulation of Phi-4, containing all the 
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


Created: 04/07/10
Created by Conor Elrick
"""
import fft
import numpy as np
import matplotlib.pyplot as plt
import cmath

np.seterr(all='raise')

class phi4(fft.FFT):
    """
     
    """
    
    def __init__(self, n, L, start_type = 'hot', lamda = 0.3,\
                 multihit = 3,  mass = 1.0j):
        self.size = (n,L)
        self.log_file = open("logfile_phi4.out", 'w')
        self.system = self.system_setup(start_type)
        self.sigma = np.random.choice([-1,1]) # direction of ecmc movement
        self.chain_length = L**2 #event chain length
        self.lifted_site = 0
        self.acceptance = 0.
        self.multihit = multihit
        #self.kappa = kappa
        self.lambda_0 = lamda
        self.mass = mass
        #try:
        #self.omega = np.sqrt(-6*self.mass**2/self.lambda_0)
        #print(self.omega)
        #except:
        self.omega = None
        #    print("[O] Omega is complex. omega -> None")
        
        
    def action_self(self, site_value):
        """
        Self-interaction term in action. Note the differing normalisation to FFT
        """
       
        #return site_value**2 + self.lambda_0*(site_value**2-1.0)**2
        return 0.5*(self.mass**2)*(site_value**2) + self.lambda_0*(site_value**4)/24.0
        
    def action_interaction(self, site_value, neighbour_value):
        """
        Interaction term in the action. Note the differing normalisation to FFT.
        """
        
        #return  -2.0*self.kappa*(site_value*neighbour_value)
        return  0.5*((site_value - neighbour_value)**2)
    
    def ecmc_t_calc(self, energy, phi):
        """
        Case B
        """
        #if self.lambda_0 *self.mass**2 >0:
        if self.omega == None:
            #Only one root
            psi_2 = - self.sigma*phi
            if psi_2 <= 0:
                
                return self.ecmc_0_solve(energy,phi)
            else:
                
                return self.ecmc_1_solve(energy,phi,psi_2)
        else:
            # 3 roots
            psi_2 = -self.sigma*phi
            psi_1 = psi_2 - self.omega
            psi_3 = psi_2 + self.omega
            if psi_3 < 0:
                #print("One")
                return self.ecmc_0_solve(energy,phi)
                
            elif psi_2 <0:
                #print("two")
                return self.ecmc_2_solve(energy,phi)
            elif psi_1 <0:
                energy_subs = np.real(-(1.0/24)*self.lambda_0*phi**4-0.5*self.mass**2*phi**2)
                #print(energy_subs)
                if energy_subs <= energy:
                    assert(self.ecmc_2_solve(energy-energy_subs,phi) > psi_3)
                    return self.ecmc_2_solve(energy-energy_subs,phi)
                else:
                    #print("four")
                    #print(energy)
                    #print(energy_subs)
                    #print(self.lambda_0)
                    #print(self.mass)
                    #print(phi)
                    #print(self.sigma)
                    #print(psi_2)
                    #print(self.ecmc_0_solve(energy,phi))
                    assert(self.ecmc_0_solve2(energy,phi)<psi_2)
                    assert(self.ecmc_0_solve2(energy,phi)>0)
                    return self.ecmc_0_solve2(energy,phi)
                    
            elif psi_1 >0:
                energy_subs = np.real(3.0*self.mass**4/(2.0*self.lambda_0))
                #print(energy_subs)
                if energy_subs <= energy:
                    #print("five")
                    assert(self.ecmc_2_solve(energy-energy_subs,phi) > psi_3)
                    return self.ecmc_2_solve(energy-energy_subs,phi)
                else:
                    #print(energy)
                    #print(energy_subs)
                    #print(self.lambda_0)
                    #print(self.mass)
                    #print(phi)
                    #print(self.sigma)
                    #print(psi_2)
                    #print(self.ecmc_2_solve2(energy,phi))
                    #assert(self.ecmc_0_solve2(energy,phi)<psi_2)
                    assert(self.ecmc_2_solve2(energy,phi)<psi_2)
                    assert(self.ecmc_2_solve2(energy,phi)>psi_1)
                    return self.ecmc_2_solve2(energy,phi)
                    
            else:
                raise FloatingPointError
                
    


            
    def ecmc_t_calc_neighbour(self, energy, delta_phi,mass):

        if self.sigma * delta_phi >= 0:
            
            return np.sqrt((2.0 * energy / (mass**2)) + delta_phi**2 ) -\
                   self.sigma*delta_phi
        else:
            return np.sqrt(2.0 * energy / (mass**2)) -\
                   self.sigma*delta_phi

    def ecmc_move(self):
        """
        Main loop for a chain of events
        """
        move_var = True
        chain_length = self.chain_length
        self.lifted_site = np.random.randint(0, self.system.size)
        while (move_var == True):
            
            neighbours = self.nearest_neighbours(self.lifted_site)

            
            moves = np.zeros([len(neighbours) + 1,2]) # +1 for self interaction

            
            #Self interaction
            delta_E = self.ecmc_energy()
            phi = self.system.flat[self.lifted_site]
            moves[0] = np.array([self.lifted_site, self.ecmc_t_calc(delta_E,\
                                phi)])

            
            #Pair interaction
            i = 1
            #   #
            for neighbour in neighbours:
                delta_phi = self.system.flat[self.lifted_site] -\
                            self.system.flat[neighbour]
                
                delta_E = self.ecmc_energy()
                
                moves[i] = np.array([neighbour, self.ecmc_t_calc_neighbour(delta_E,\
                                    delta_phi, 1.0)])
                                 
                i += 1
            #print(moves)    
            t_move = moves[np.argmin(np.transpose(moves)[1], axis=0)]
            assert(t_move[1] > 0)
                
            if t_move[1] < self.chain_length:
                
                
                    
                self.chain_length -= t_move[1]
                
                self.system.flat[self.lifted_site] += self.sigma * t_move[1]
                    
                if t_move[0] == self.lifted_site:    
                
                    self.sigma *= -1

                self.lifted_site = int(t_move[0])

            else:
                self.system.flat[self.lifted_site] += self.sigma *\
                                                (self.chain_length)
                if t_move[0] == self.lifted_site:
                    self.sigma *= -1

                move_var = False
                self.chain_length = chain_length

           
        return None
    
    def ecmc_0_solve(self, energy, phi):
        L = self.lambda_0
        M = self.mass**2
        S = self.sigma
        return -S*phi + (1.0/L)*np.sqrt(-6*L*M +L*np.sqrt(L**2*phi**4 + 12*M*L*phi**2 + 36*M**2 + 24*energy*L))
    
    def ecmc_0_solve2(self, energy, phi):
        L = self.lambda_0
        M = self.mass**2
        S = self.sigma
        return -S*phi - (1.0/L)*np.sqrt(-6*L*M -L*np.sqrt(L**2*phi**4 + 12*M*L*phi**2 + 36*M**2 + 24*energy*L))
    
    def ecmc_1_solve(self, energy, phi, w):
        L = self.lambda_0
        M = self.mass**2
        S = self.sigma
        return -S*phi + (1.0/L)*np.sqrt(-6*L*M +2*L*np.sqrt(9*M**2 + 6*energy*L))
        
    def ecmc_2_solve(self, energy, phi):
        L = self.lambda_0
        M = self.mass**2
        S = self.sigma
        return -S*phi + (1.0/L)*np.sqrt(2*L*np.sqrt(6*energy*L)-6*L*M)
     
    def ecmc_2_solve2(self, energy, phi):
        L = self.lambda_0
        M = self.mass**2
        S = self.sigma
        return -S*phi - (1.0/L)*np.sqrt(-2*L*np.sqrt(6*energy*L)-6*L*M)
        
    def mc_local_step(self,site):
        """
        MCMC move given a site and an action. The neighbours and self are
        split in two cases.
        """
        
        old_site_value = self.system.flat[site]
        new_site_value = old_site_value + np.random.normal(loc=0.,scale=0.8)
        
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
        try:
            P_acceptance = min([1, np.exp(old_action - new_action)])
        except FloatingPointError as error:
            #print(error)
            P_acceptance = 0.

        gamma = np.random.uniform(0,1)
        
        if P_acceptance <= gamma:
            move_acceptance = False
        
        #Acceptance
        if move_acceptance == True:
            self.acceptance += 1.
            self.system.flat[site] = new_site_value
           
        
        return None
    def newobs(self):
        tot = 0.
        for site in self.system.flat:
            tot += site
        return tot
       
    
        
        
if __name__ == '__main__':

        lattice = phi4(2,32, lamda = 100,mass=9.5j)
        #max 3.2j
        #min 2.5j
        steps = int(500)
        newobs = []
        for i in range(steps):
      
            lattice.mcmc_move()
            newobs.append(lattice.newobs())
        print(np.average(newobs))
        
        
        
        lattice2 = phi4(2,32, lamda = 100,mass=9.5j)
        #max 3.2j
        #min 2.5j
        steps = int(500)
        newobs2 = []
        for i in range(steps):
      
            lattice2.mcmc_move()
            newobs2.append(lattice2.newobs())
        print(np.average(newobs2))
        
        lattice3 = phi4(2,32, lamda = 100,mass=9.5j)
        #max 3.2j
        #min 2.5j
        steps = int(500)
        newobs3 = []
        for i in range(steps):
      
            lattice3.mcmc_move()
            newobs3.append(lattice3.newobs())
        print(np.average(newobs3))
        
        plt.plot(np.arange(steps), newobs)
        plt.plot(np.arange(steps), newobs2)
        plt.plot(np.arange(steps), newobs3)
        plt.show()
