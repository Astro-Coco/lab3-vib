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
Longueur = 
Largeur = 
Epaisseur = 
Paroi = 
E = 
RHO = 
A = 
I = 

# Calcul des fréquences propres et modes propres théorique (sans moteur ni carburant)
# Nombre de modes à calculer
N = 10
BL = np.zeros(N)
# Calcul de Betak.L, un vecteur [1 x N] (%%% A compléter %%%)
# Trouver les beta.L (BL) en utilisant fsolve et une estimée initiale des BL







# Discrétisation de l'aile en 100 points
Nb_Points = 100
x = np.linspace(0,Longueur,Nb_Points)

# Initialisation du vecteur OMEGA_theorique
OMEGA_theorique = np.zeros((N,1))
# Initialisation de la matrice v_theorique
v_theorique = np.zeros((Nb_Points,N))

# Calcul des modes propres théoriques (%%% A compléter %%%)
for i in range(N):
    OMEGA_theorique[i] = 
    print('OMEGA_theorique[',i,']',OMEGA_theorique[i])
    v_theorique[:,i] =  # appel à la routine Xi dans le module Xi_module.py (àcompléter aussi)
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
