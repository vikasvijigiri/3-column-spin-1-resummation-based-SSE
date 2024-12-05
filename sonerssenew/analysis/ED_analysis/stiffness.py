import numpy as np
from   scipy.sparse import csr_matrix
import scipy as sp
from   scipy.sparse.linalg import eigsh
import itertools as it


H_BONDS = 4
Q_BONDS = 4
BQ_BONDS = 2



def basis_no(i, L):
    return sum((2**j) * ((i[j] + 1) // 2) for j in range(L))


def action(ii, Hbond, theta):
    i, j = Hbond
    state_list = []
    coef_list = []
    
    # Heisenberg AFM only model
    # diagonal
    de = -0.5**2 * (1 - (ii[i] * ii[j])) * np.cos(theta)**2
    if abs(de) > 1e-6:
      coef_list.append(de)
      state_list.append(ii)

    # off-diagonal
    if (ii[i]*ii[j] == -1):
      jj = list(ii)
      jj[i] = -ii[i] # flip site i
      jj[j] = -ii[j] # flip site j
      jj = tuple(jj)
      coef_list.append(0.5  * (1 + np.cos(theta))/2.)
      state_list.append(jj)

    return state_list, coef_list   



# Spin-1/2
def Bilayer_Bi_QQ_JH(L, JH, JB, JQ, theta, Hbonds, Bplaqs, QQplaqs):

    data = []
    rows = []
    cols = []

    for ii in it.product([-1, 1], repeat=L): 

        # Heisenberg interactions
        if abs(JH) > 1e-6:
          for bond in Hbonds:
              state_list, coef_list = action(ii, bond, theta)
              for jc, jj in enumerate(state_list):
                  data.append(JH * coef_list[jc])
                  rows.append(basis_no(ii, L))
                  cols.append(basis_no(jj, L)) 

        # Biquad interactions
        if abs(JB) > 1e-6:
          for bond1, bond2 in Bplaqs:
              state_list, coef_list = action(ii, bond1, theta)

              for jc, jj in enumerate(state_list):
                  f_state_list, f_coef_list = action(jj, bond2, theta)

                  for kc, kk in enumerate(f_state_list):
                      data.append(-JB * 2. * coef_list[jc] * f_coef_list[kc])
                      rows.append(basis_no(ii, L))
                      cols.append(basis_no(kk, L)) 

        # QQ interactions
        if abs(JQ) > 1e-6:
          for bond1, bond2, bond3, bond4 in QQplaqs:
              state_list, coef_list = action(ii, bond1, theta)

              for jc, jj in enumerate(state_list):
                  j_state_list, j_coef_list = action(jj, bond2, theta)

                  for kc, kk in enumerate(j_state_list):
                      k_state_list, k_coef_list = action(kk, bond3, theta)

                      for lc, ll in enumerate(k_state_list):
                          l_state_list, l_coef_list = action(ll, bond4, theta)

                          for mc, mm in enumerate(l_state_list):
                              data.append(-JQ * (2.**2) * coef_list[jc] * j_coef_list[kc] * k_coef_list[lc] * l_coef_list[mc])
                              rows.append(basis_no(ii, L))
                              cols.append(basis_no(mm, L)) 

    
    Hsp = csr_matrix((data, (rows, cols)), shape=(2**L, 2**L), dtype=np.float64)     
    return Hsp 






def set_H_Bond(bonds, index, s1, s2, N):
    bonds[index][0] = s1
    bonds[index][1] = s2

    bonds[index + 1][0] = s1 + N
    bonds[index + 1][1] = s2

    bonds[index + 2][0] = s1
    bonds[index + 2][1] = s2 + N

    bonds[index + 3][0] = s1 + N
    bonds[index + 3][1] = s2 + N


def set_QQ_Bonds(bonds, baseIndex, s1, s2, s3, s4, offset, N, shift):
    bonds[baseIndex + offset, 0, 0] = s1 + N * shift[0]
    bonds[baseIndex + offset, 0, 1] = s2 + N * shift[1]
    bonds[baseIndex + offset, 1, 0] = s1 + N * shift[2]
    bonds[baseIndex + offset, 1, 1] = s2 + N * shift[3]

    bonds[baseIndex + offset, 2, 0] = s3 + N * shift[4]
    bonds[baseIndex + offset, 2, 1] = s4 + N * shift[5]
    bonds[baseIndex + offset, 3, 0] = s3 + N * shift[6]
    bonds[baseIndex + offset, 3, 1] = s4 + N * shift[7]



def set_Bi_Bonds(bonds, baseIndex, s1, s2, offset, N, shift):
    bonds[baseIndex + offset, 0, 0] = s1 + N * shift[0]
    bonds[baseIndex + offset, 0, 1] = s2 + N * shift[1]
    bonds[baseIndex + offset, 1, 0] = s1 + N * shift[2]
    bonds[baseIndex + offset, 1, 1] = s2 + N * shift[3]


def gen_bonds(lx, ly):

    Ns = lx*ly
    JHsites = np.zeros((8 * Ns, 2), dtype=int)
    JQsites = np.zeros((8 * Ns, 4, 2), dtype=int)
    JBsites = np.zeros((4 * Ns, 2, 2), dtype=int)

    for y1 in range(ly):
        for x1 in range(lx):
            # ***************************************************************************************
            #       Heisenberg interaction bonds
            # ***************************************************************************************

            # Horizontal bonds
            s1 = x1 + y1 * lx
            x2 = (x1 + 1) % lx
            s2 = x2 + y1 * lx

            set_H_Bond(JHsites, H_BONDS * s1, s1, s2, Ns)

            # vertical bonds
            x2 = x1
            y2 = (y1 + 1) % ly
            s2 = x2 + y2 * lx

            set_H_Bond(JHsites, H_BONDS * (Ns + s1), s1, s2, Ns)

            # ***************************************************************************************
            #       QoQ interaction bonds
            # ***************************************************************************************
            # Horizontal bonds
            s1 = x1 + y1 * lx
            x2 = (x1 + 1) % lx
            s2 = x2 + y1 * lx

            s3 = (s1 + lx) % Ns
            s4 = (s2 + lx) % Ns

            baseIndexHorizontal = 4 * s1

            shiftValues = [0, 0, 1, 1, 0, 0, 1, 1]
            set_QQ_Bonds(JQsites, baseIndexHorizontal, s1, s2, s3, s4, 0, Ns, shiftValues)

            shiftValues = [0, 0, 1, 1, 0, 1, 1, 0]
            set_QQ_Bonds(JQsites, baseIndexHorizontal, s1, s2, s3, s4, 1, Ns, shiftValues)

            shiftValues = [0, 1, 1, 0, 0, 0, 1, 1]
            set_QQ_Bonds(JQsites, baseIndexHorizontal, s1, s2, s3, s4, 2, Ns, shiftValues)

            shiftValues = [0, 1, 1, 0, 0, 1, 1, 0]
            set_QQ_Bonds(JQsites, baseIndexHorizontal, s1, s2, s3, s4, 3, Ns, shiftValues)

            # Vertical bonds
            s1 = x1 + y1 * lx
            s2 = (s1 + lx) % Ns

            s3 = x2 + y1 * lx
            s4 = (s3 + lx) % Ns

            # Vertical bonds
            baseIndexVertical = 4 * Ns + 4 * s1

            shiftValues = [0, 0, 1, 1, 0, 0, 1, 1]
            set_QQ_Bonds(JQsites, baseIndexVertical, s1, s2, s3, s4, 0, Ns, shiftValues)

            shiftValues = [0, 0, 1, 1, 0, 1, 1, 0]
            set_QQ_Bonds(JQsites, baseIndexVertical, s1, s2, s3, s4, 1, Ns, shiftValues)

            shiftValues = [0, 1, 1, 0, 0, 0, 1, 1]
            set_QQ_Bonds(JQsites, baseIndexVertical, s1, s2, s3, s4, 2, Ns, shiftValues)

            shiftValues = [0, 1, 1, 0, 0, 1, 1, 0]
            set_QQ_Bonds(JQsites, baseIndexVertical, s1, s2, s3, s4, 3, Ns, shiftValues)

            # ***************************************************************************************
            #       Bi interaction bonds
            # ***************************************************************************************
            # Horizontal bonds
            s1 = x1 + y1 * lx
            x2 = (x1 + 1) % lx
            s2 = x2 + y1 * lx

            baseIndexHorizontal = 2 * s1

            shiftValues = [0, 0, 1, 1]
            set_Bi_Bonds(JBsites, baseIndexHorizontal, s1, s2, 0, Ns, shiftValues)

            shiftValues = [0, 1, 1, 0]
            set_Bi_Bonds(JBsites, baseIndexHorizontal, s1, s2, 1, Ns, shiftValues)

            x2 = x1
            y2 = (y1 + 1) % ly
            s2 = x2 + y2 * lx

            baseIndexVertical = 2 * Ns + 2 * s1

            shiftValues = [0, 0, 1, 1]
            set_Bi_Bonds(JBsites, baseIndexVertical, s1, s2, 0, Ns, shiftValues)

            shiftValues = [0, 1, 1, 0]
            set_Bi_Bonds(JBsites, baseIndexVertical, s1, s2, 1, Ns, shiftValues)

    return JHsites, JQsites, JBsites






def Bilayer_JH_Bi_QQ_ED(Lx, Ly, Nl, JH, JB, JQ, theta, Hbnds, Bplqs, QQplqs, n, sparse=False):
    #if n < 30: sparse = True
    Hsp = Bilayer_Bi_QQ_JH(Lx*Ly*Nl, JH, JB, JQ, theta, Hbnds, Bplqs, QQplqs)
    print("Hamiltonian Built success")
    if not sparse: Wsp, Vsp = sp.linalg.eigh(Hsp.toarray())
    else: Wsp = eigsh(Hsp, k=n, which='SA',return_eigenvectors=False)
    print("Diagonalization -> success!")
    return Wsp



# *********************************** Measurements **********************************

def zMag_square(Lx, Lz, Nl, basis_list, state_vector):
    smag = 0.; smag2 = 0.; smag4 = 0.;
    for i, basis in enumerate(basis_list):
        mg = 0.        
        for x in range(Lx*Lz):
            mg += 0.5 * (basis[x] + basis[x + Lx*Lz]) * (-1) ** ((x // Lx) + (x % Lx))
        smag += (mg)*(state_vector[i] ** 2)
        smag2 += (mg**2)*(state_vector[i] ** 2)
        smag4 += (mg**4)*(state_vector[i] ** 2)
    return smag/(Lx*Lz*Nl), smag2/(Lx*Lz*Nl)**2, smag4/(Lx*Lz*Nl)**4        
