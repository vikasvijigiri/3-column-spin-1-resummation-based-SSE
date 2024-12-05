# ***********************************************************************************
#                              Exact diagonalization 
# Description: (Heisenberg + three-spin) model on square lattice for RSSE
# Author: Vikas Vijigiri
# Date: October, 2023
# ***********************************************************************************
import numpy as np
from   scipy.sparse import csr_matrix
import scipy as sp
from   scipy.sparse.linalg import eigsh
import itertools as it
import pandas as pd
import concurrent.futures


H_BONDS = 4
Q_BONDS = 4
BQ_BONDS = 2



# ***********************************************************************************
#                               Hamiltonian part 
# ***********************************************************************************


# basis_rep
def basis_no(basis_state, L):
    return sum((3**j) * ((basis_state[j] + 1)) for j in range(L))


# spin-1/2 action
def action(ii, Hbond):
    i, j = Hbond
    state_list = []
    coef_list = []
    
    # Heisenberg AFM only model
    # diagonal
    de = - (1 - (ii[i] * ii[j]))
    if abs(de) > 1e-6:
      coef_list.append(de)
      state_list.append(ii)

    # off-diagonal
    
    # S_+i*s_-j
    m, local_state = pmflip(ii[i], ii[j])
    states = flipper(Hbond, local_state, ii)
    if abs(m) > 1e-6:
      coef_list.append(m)
      state_list.append(states)
    
    # S_-i*s_+j
    m, local_state = mpflip(ii[i], ii[j])
    states = flipper(Hbond, local_state, ii)
    if abs(m) > 1e-6:
      coef_list.append(m)
      state_list.append(states)

    return state_list, coef_list 


# you can edit here
def pmflip(l0, l1):  # this is OK
    return [1, [l0+1, l1-1]] if l0 != 1 and l1 != -1  else [0, [l0,l1]]

def mpflip(l0, l1):  # this is OK
    return [1, [l0-1, l1+1]] if l0 != -1 and l1 != 1  else [0, [l0,l1]]   
 
    
def flipper(bond, loc_state, ii):
    jj = ii                
    jj = list(jj)
    jj[bond[0]]=loc_state[0] 
    jj[bond[1]]=loc_state[1] 
    jj = tuple(jj)
    return jj



# Spin-1
def H_3spin_ED(L, JH, J3, Hbonds, threeS_Spinbnds, threeS_Hbnds, threeS_NNHbnds):

    data = []
    rows = []
    cols = []


    for ii in it.product([-1, 0, 1], repeat=L): 

        # Heisenberg interactions
        if abs(JH) > 1e-6:
          for i, bond in enumerate(Hbonds):
              #print(i, bond, HHsgn[i])
              state_list, coef_list = action(ii, bond)
              for jc, jj in enumerate(state_list):
                  data.append(JH * coef_list[jc])
                  rows.append(basis_no(ii, L))
                  cols.append(basis_no(jj, L)) 

        # 3-spin interactions
        if abs(J3) > 1e-6:
          for i, (bond1, bond2) in enumerate(threeS_Spinbnds):
              state_list, coef_list = action(ii, bond1)

              for jc, jj in enumerate(state_list):
                  f_state_list, f_coef_list = action(jj, bond2)

                  for kc, kk in enumerate(f_state_list):
                      data.append(-J3 * coef_list[jc] * f_coef_list[kc])
                      rows.append(basis_no(ii, L))
                      cols.append(basis_no(kk, L)) 


        # additional H-bonds from three-spin interactions (N and NN)
        if abs(J3) > 1e-6:
          for i, bond in enumerate(threeS_Hbnds + threeS_NNHbnds):
              #print(i, bond, HHsgn[i])
              state_list, coef_list = action(ii, bond)
              for jc, jj in enumerate(state_list):
                  data.append(J3 * coef_list[jc])
                  rows.append(basis_no(ii, L))
                  cols.append(basis_no(jj, L)) 

    Hsp = csr_matrix((data, (rows, cols)), shape=(3**L, 3**L), dtype=np.float64)     
    return Hsp 






def set_H_Bonds(bonds, baseIndex, s1, s2, offset, N, shift, sgn, lamda): 
    bonds[baseIndex + offset][0] = s1 + N * shift[0]
    bonds[baseIndex + offset][1] = s2 + N * shift[1]
    sgn[baseIndex + offset] = lamda


def set_Bi_Bonds(bonds, baseIndex, s1, s2, offset, N, shift, sgn, lamda):
    bonds[baseIndex + offset, 0, 0] = s1 + N * shift[0]
    bonds[baseIndex + offset, 0, 1] = s2 + N * shift[1]
    bonds[baseIndex + offset, 1, 0] = s1 + N * shift[2]
    bonds[baseIndex + offset, 1, 1] = s2 + N * shift[3]
    sgn[baseIndex + offset] = lamda



# comment: it's a 4-site chain problem
def gen_bonds(lx, ly):


    Ns = lx*ly
    JHsites = np.zeros((Ns, 2), dtype=int)
    J3sites = np.zeros((Ns, 2, 2), dtype=int)
    J3Hbonds = np.zeros((2 * Ns, 2), dtype=int)
    J3NNHbonds = np.zeros((Ns, 2), dtype=int)

    Hsgn = np.zeros((Ns), dtype=np.float128)
    Bsgn = np.zeros((Ns), dtype=np.float128)


    JHsites = [[0, 1], [1, 2], [2, 3], [3, 0]]
    J3sites = [[[0, 1], [1, 2]], [[1, 2], [2, 3]], [[2, 3], [3, 0]], [[3, 0], [0, 1]]]
    J3Hbonds = [[0, 1], [1, 2], [2,3], [3,0], [0,1], [1,2], [2,3], [3,0]]
    J3NNHbonds = [[0, 2], [1, 3], [2, 0], [3, 1]]

#    for y1 in range(ly):
#        for x1 in range(lx):
#            # ***************************************************************************************
#            #       Heisenberg interaction bonds
#            # ***************************************************************************************

#            # Horizontal bonds
#            s1 = x1 + y1 * lx
#            x2 = (x1 + 1) % lx
#            s2 = x2 + y1 * lx

#            #set_H_Bond(JHsites, H_BONDS * s1, s1, s2, Ns)
#            baseIndexHorizontal = s1;

#            shiftValues = [0, 0]; # lambda = 1 
#            set_H_Bonds(JHsites, baseIndexHorizontal, s1, s2, 0, Ns, shiftValues, Hsgn, 1.);            



#            # vertical bonds
#            x2 = x1
#            y2 = (y1 + 1) % ly
#            s2 = x2 + y2 * lx

#            #set_H_Bond(JHsites, H_BONDS * (Ns + s1), s1, s2, Ns)
#            baseIndexVertical = Ns + s1;

#            shiftValues = [0, 0]; # lambda = 1
#            set_H_Bonds(JHsites, baseIndexVertical, s1, s2, 0, Ns, shiftValues, Hsgn, 1.);


#            # ***************************************************************************************
#            #       3-spin interaction bonds
#            # ***************************************************************************************
#            // Horizontal bonds
#            s1 = x1 + y1 * lx;
#            x2 = (x1 + 1) % lx;
#            s2 = x2 + y1 * lx;

#            s3 = s2;
#            x3 = (x1 + 2) % lx;
#            s4 = x3 + y1 * lx;

#            baseIndexHorizontal = 2 * s1

#            shiftValues = [0, 0, 1, 1]
#            set_Bi_Bonds(J3sites, baseIndexHorizontal, s1, s2, 0, Ns, shiftValues, Bsgn, 1.);                 

#            shiftValues = [0, 1, 1, 0]
#            set_Bi_Bonds(J3sites, baseIndexHorizontal, s1, s2, 1, Ns, shiftValues, Bsgn, lmbda)

#            # Vertical bonds
#            x2 = x1
#            y2 = (y1 + 1) % ly
#            s2 = x2 + y2 * lx

#            baseIndexVertical = 2 * Ns + 2 * s1

#            shiftValues = [0, 0, 1, 1]
#            set_Bi_Bonds(J3sites, baseIndexVertical, s1, s2, 0, Ns, shiftValues, Bsgn, 1.);

#            shiftValues = [0, 1, 1, 0]
#            set_Bi_Bonds(J3sites, baseIndexVertical, s1, s2, 1, Ns, shiftValues, Bsgn, lmbda);

    return JHsites, J3sites, J3Hbonds, J3NNHbonds #, Hsgn, Bsgn, QQsgn






def Heisenberg_threeSpin_ED(Lx, Ly, JH, J3, Hbonds, threeS_Spinbnds, threeS_Hbnds, threeS_NNHbnds, n, sparse=False):
    #if n < 30: sparse = True
    Hsp = H_3spin_ED(Lx*Ly, JH, J3, Hbonds, threeS_Spinbnds, threeS_Hbnds, threeS_NNHbnds)
    print("Hamiltonian Built -> success!")
    if not sparse: Wsp, Vsp = sp.linalg.eigh(Hsp.toarray())
    else: Wsp, Vsp = eigsh(Hsp, k=n, which='SA',return_eigenvectors=True)
    print("Diagonalization -> success!")
    return Wsp+J3, Vsp


# ***********************************************************************************
#                               Measurements 
# ***********************************************************************************

class Site:
    def __init__(self, x, y):
        self.x = x
        self.y = y

def calculate_position_vectors(Lx, Ly):
    positions = []

    # Reference position for site 0
    ref_x, ref_y = 0, 0

    # Calculate position vectors
    for i in range(Lx * Ly):
        rel_x = i % Lx
        rel_y = i // Lx
        x = rel_x - ref_x
        y = rel_y - ref_y

        # Apply periodic boundary conditions along x
        if x > Lx / 2:
            x -= Lx
        if x < -Lx / 2:
            x += Lx

        # Apply periodic boundary conditions along y
        if y > Ly / 2:
            y -= Ly
        if y < -Ly / 2:
            y += Ly

        positions.append(Site(x, y))

    return positions


def zMag_square(Lx, Ly, beta, basis_list, W, V):
    smag_square = smag_four = 0.;

    z = np.sum(np.exp(-np.float128(beta * W)))
  
    for ii, w in enumerate(W):
        state_vector = V[:,ii]
        weight = np.exp(-np.float128(beta * w))
        for ii, basis in enumerate(basis_list):
            si = state_vector[ii] ** 2
            mg = 0. # layer 1
            for x in range(Lx*Ly):
                mg += 0.5 * basis[x] * (-1) ** ((x // Lx) + (x % Lx))

            smag_square += ((mg/Lx/Ly)**2) * si * weight/z
            smag_four += ((mg/Lx/Ly)**4) * si * weight/z

    return smag_square, smag_four







def estimate_all_observables(Lx, Ly, JH, J3, Wsp, Vsp, n, beta):
    # Partition function
    z = np.sum(np.exp(-np.float128(beta * Wsp)))

    # Energy
    emin = Wsp.min()/Lx/Ly  
    enrg = np.sum((Wsp/Lx/Ly)*np.exp(-np.float128(beta * Wsp)))/z  
    enrg2 = np.sum((Wsp/Lx/Ly)**2 *np.exp(-np.float128(beta * Wsp)))/z  
    enrg4 = np.sum((Wsp/Lx/Ly)**4 *np.exp(-np.float128(beta * Wsp)))/z  



    basis_list = list(it.product([-1, 1], repeat=Lx*Ly))
    
    
    # Staggered Magnetization, Staggered Magnetization^2, Staggered Magnetization^4
    SMag_square, SMag_four = zMag_square(Lx, Ly, beta, basis_list, Wsp, Vsp)

    # Printing the observables in one go using pandas Dataframe (With clear headers)
    headers = [ "Lx", "Ly", "beta", "J_H", "J3", 
                "enrg", "enrg2", "enrg4", 
                "SMag_square", "SMag_four"
              ]

    values  = [ Lx, Ly, beta, JH, J3,
                enrg, enrg2, enrg4,
                SMag_square, SMag_four
              ]
    # Reshape values to have a single column
    values_reshaped = [[val] for val in values]

    return pd.DataFrame(values_reshaped, columns=['Values'], index=headers)
