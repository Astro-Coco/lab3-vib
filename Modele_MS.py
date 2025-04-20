import scipy.integrate as integrate
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import math as mt
import scipy.linalg as la
import sys
from scipy.optimize import fsolve
from Xi_module import Xi, Xij, Xij_dd

# use LaTeX fonts in the plot
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Dimension de l'aile et propriétés  (%%% A compléter %%%)
Longueur = 10
Largeur = 2
Epaisseur = 0.3
Paroi = 0.05
E = 70e9
RHO = 2700
A = (Largeur*Epaisseur) - (Largeur-2*Paroi)*(Epaisseur-2*Paroi)
I = 2*( (1/12)*(Largeur*Paroi**3) + Largeur*Paroi*(Epaisseur/2-Paroi/2)**2 ) + 2*((1/12)*Paroi*(Epaisseur-2*Paroi)**3)

# Moteur
LM = 4 # Emplacement du moteur en m
MM = 1000. # masse du moteur en kg

# Carburant
LC1 = 2 # Début de la zone de carburant en m
LC2 = 6 # Fin de la zone de carburant en m
MC = 304. # Masse par unité de longueur (kg/m)

# Calcul des fréquences propres et modes propres théorique (sans moteur ni carburant)
# Nombre de modes à calculer
N = 10
BL = np.zeros(N)
# Calcul de Betak.L, un vecteur [1 x N] (%%% A compléter %%%)
for i in range(N):
    BL[i] = (i+0.5)*np.pi

def func(x):
    return np.cos(x)+1/np.cosh(x)

for i in range(N):
    BL[i] = BL[i] = fsolve(func, BL[i])[0]


# Discrétisation de l'aile en 100 points
Nb_Points = 100
x = np.linspace(0,Longueur,Nb_Points)

# Calul de mij et kij      (%%% A compléter %%%)
m = RHO*A
mij = np.zeros((N,N))
kij = np.zeros((N,N))
for i in range(N):
    for j in range(i,N):
        # Masse ponctuelle du moteur à LM
        xi_LM = Xi([LM], BL[i], Longueur)
        xj_LM = Xi([LM], BL[j], Longueur)
        mij[i,j] += MM * xi_LM[0] * xj_LM[0]
        
        # masse de l’aile
        mij[i,j], _ = integrate.quad(lambda xx: m * Xij(np.array([xx]), BL[i], BL[j], Longueur), 0, Longueur, limit = 200)

        # masse carburant
        mij[i,j] += integrate.quad(lambda xx: MC * Xij(np.array([xx]), BL[i], BL[j], Longueur), LC1, LC2, limit = 200)[0]

        # raideur
        kij[i,j], _ = integrate.quad(lambda xx: E * I * Xij_dd(np.array([xx]), BL[i], BL[j], Longueur), 0, Longueur, limit = 200)

        # mij[i,j] = mij[i,j] si j ne va que de i-1 à N-1: range(i,N)
        # kij[i,j] = kij[i,j] si j ne va que de i-1 à N-1: range(i,N)

        # Symétrie de la matrice, donc j in range(i,N)

# Calul des fréquences propres    (%%% A compléter %%%)
LAMBDA, PHI = la.eig(kij, mij)
# print("LAMBDA : Les valeurs propres")
# print("PHI normalisée (norme euclidienne)")
OMEGA = sorted(np.sqrt(np.real(LAMBDA)))
for i in range(N):
    #LAR = PHI[i,:]
    PHI[i,:] = PHI[i,:][np.sqrt(np.real(LAMBDA)).argsort()]

def Modes_MS(Q,x,Longueur,BL):
    from Xi_module import Xi
    # Description
    #
    # Cette fonction permet de calculer les modes propres de votre modèle à
    # partir de la matrice de la matrice Q obtenue par la méthode des modes
    # supposées
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # %% Variables d'entrée
    # %
    # % Q: Matrice de coeffcients obtenue par la méthode des modes supposées
    # % x: Coordonnées des points de discrétisation de l'aile
    # % Longueur: Longueur de l'aile
    # % BkL: Vecteur BL
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Variables de sortie
    #
    # v_AM: Matrice [Nb_Points x N]. Chaque colonne de la matrice
    # correspond à un mode propre en fonction de x.
    #
    # Nombre de modes calculés
    N = len(Q)
    # Calcul du nombre de points de discrétisation
    Nb_Points = len(x)
    # Initialisation de la matrice v_AM
    v_AM = np.zeros((Nb_Points,N));
    # Calcul des vecteur propres de l'aile
    for i in range(N):
        for n in range(N):
            v_AM[:,i] = v_AM[:,i] + Q[n,i]*Xi(x,BL[n],Longueur)
        # Normalisation du vecteur propre
        v_AM[:,i] = abs(v_AM[2,i])/(v_AM[2,i])*(v_AM[:,i]/np.linalg.norm(v_AM[:,i]))
    return v_AM

# Calcul de modes propres    (%%% A compléter %%%)
v_AM = Modes_MS(PHI, x, Longueur, BL)

# Tracage des 5 premiers Modes et impression des fréquences en rad/s
Nb_Modes = 5

print("Les fréquences propres en rad/s sont")
for i in range(N):
    print(OMEGA[i])
    
couleur = iter(plt.cm.rainbow(np.linspace(0,0.3,Nb_Modes)))

for i in range(Nb_Modes):
    c=next(couleur)
    plt.plot(x[:],v_AM[:,i],'-',label='mode '+str(i+1), color=c)

plt.xlabel(r'x/L', fontsize=11)
plt.ylabel(r'Amplitude normée', fontsize=11)
plt.legend(loc='best')
plt.grid(True)
plt.savefig('modes_supposes.png', dpi = 300, bbox_inches='tight')
plt.show()
