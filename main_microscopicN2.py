import numpy as np
#import matplotlib.pyplot as plt
#import scipy as sp
# Logging module
import logging
#logging.basicConfig(filename='example.log', encoding='utf-8', level=logging.DEBUG)
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
#==============================================================================
# Initial Conditions Class
class init:
    def __init__(self, N, electronNumber, muonNumber):
        self.N= N
        self.electronNumber= electronNumber
        self.muonNumber= muonNumber
        
        self.nue_state= np.array([[1], [0]])
        self.numu_state= np.array([[0], [1]])
        
        self.mu0= 1
        
        logging.info(f"Number of electron neutrino = {self.electronNumber}")
        logging.info(f"Number of muon neutrino = {self.muonNumber}")
        
        # Check if the number of particles is correct
        if self.electronNumber+ self.muonNumber!= self.N:
            # Throw an error
            logging.error("Total number of particles is not equal to N. N = nu_e + nu_mu")
            exit()
        
        # Check Total Number of Particles
        if self.N != 2:
            # Throw an error
            logging.error("This script is only for 2 particles.")
            exit()
#==============================================================================
# Constants Class
class const:
    sigma_mat= np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]])
    gamma_matrices= np.zeros((4, 4, 4), dtype= complex)
    gamma_matrices[0]= [[1,0,0,0], [0,1,0,0], [0,0,-1,0], [0,0,0,-1]]
    for it in range(3):
        # Combine [[0 , sigma_mat],[-sigma_mat, 0]] the 0 is 2x2 and sigma_mat is 2x2
        gamma_matrices[it+1]= np.concatenate((np.concatenate((np.zeros((2, 2), dtype= complex), sigma_mat[it]), axis= 1),\
                                              np.concatenate((-sigma_mat[it], np.zeros((2, 2), dtype= complex)), axis= 1)), axis= 0)
#==============================================================================
# Create Initial States of Levels
def state_level1(init):
    # Create a state vector
    state= np.zeros((init.N, 1), dtype= complex)
    state= np.kron(init.nue_state, init.numu_state)
    return state
#==============================================================================
def hamiltonian_level1(init):
    # Hamiltonian is sigma_A1 (dot) sigma_B1
    
    hamilt= np.zeros((4, 4), dtype= complex)
    a= np.zeros((4, 4), dtype= complex)
    # https://sci-hub.se/https://www.sciencedirect.com/science/article/abs/pii/037026939291887F?via%3Dihub
    for it in range(3):
        hamilt= hamilt+ np.kron(const.sigma_mat[it], np.eye(2))\
                      @ np.kron(np.eye(2), const.sigma_mat[it])
        a= a+ np.kron(const.sigma_mat[it], const.sigma_mat[it])
    print("sigma_x (x) sigma_x + sigma_y (x) sigma_y + sigma_z (x) sigma_z")
    print(a)
    print("="*5)
    print("[sigma_x (x) I] * [I (x) sigma_x] + [sigma_y (x) I] * [I (x) sigma_y] + [sigma_z (x) I] * [I (x) sigma_z]")
    print(hamilt)
    exit()
    return hamilt
#==============================================================================
if __name__ == '__main__':
    # Initial Conditions
    init= init(N=2, electronNumber=1, muonNumber=1)
    state= state_level1(init)
    print(hamiltonian_level1(init))
#==============================================================================