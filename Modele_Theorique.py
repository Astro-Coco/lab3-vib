import scipy.integrate as integrate
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import math as mt
import sys
from scipy.optimize import fsolve
from Xi_module import Xi

# lint =  integrate.quad(integrand,0., 20. ,args=(a,r,xp,xl))[0]

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
I = 2*( (1/12)*(Largeur*Paroi**3) + Largeur*Paroi*(Epaisseur/2 - Paroi/2)**2 ) + 2*((1/12)*Paroi*(Epaisseur-2*Paroi)**3)


# Calcul des fréquences propres et modes propres théorique (sans moteur ni carburant)
# Nombre de modes à calculer
N = 10
BL = np.zeros(N)
# Calcul de Betak.L, un vecteur [1 x N] (%%% A compléter %%%)
# Trouver les beta.L (BL) en utilisant fsolve et une estimée initiale des BL

def eq_beta(betaL):
    return np.cos(betaL) + 1/np.cosh(betaL)

# Discrétisation de l'aile en 100 points
Nb_Points = 100
x = np.linspace(0,Longueur,Nb_Points)

# Initialisation du vecteur OMEGA_theorique
OMEGA_theorique = np.zeros((N,1))
# Initialisation de la matrice v_theorique
v_theorique = np.zeros((Nb_Points,N))

# Calcul des modes propres théoriques (%%% A compléter %%%)
for i in range(N):
    guess = (i + 0.5) * np.pi
    BL[i] = fsolve(eq_beta, guess)
    beta = BL[i] / Longueur
    print('BL[',i,'] : ',BL[i])
    OMEGA_theorique[i] = (beta**2) * np.sqrt(E * I / (RHO * A))
    print('OMEGA_theorique[',i,'] : ',OMEGA_theorique[i], '\n')
    v_theorique[:,i] = Xi(x, BL[i], Longueur) # appel à la routine Xi dans le module Xi_module.py
    # Normalisation du vecteur propre
    v_theorique[:,i] = abs(v_theorique[2,i])/(v_theorique[2,i])*(v_theorique[:,i]/np.linalg.norm(v_theorique[:,i]))

# Tracage des 5 premiers Modes
Nb_Modes = 5
couleur = iter(plt.cm.rainbow(np.linspace(0,0.3,Nb_Modes)))

for i in range(Nb_Modes):
    c=next(couleur)
    plt.plot(x[:],v_theorique[:,i],'-',label='mode '+str(i+1), color=c)

plt.xlabel(r'x/L', fontsize=11)
plt.ylabel(r'Amplitude normée', fontsize=11)
plt.legend(loc='best')
plt.grid(True)
plt.savefig('modes_theoriques.png', dpi = 300, bbox_inches='tight')
plt.show()
exit(0)
