import scipy.integrate as integrate
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import math as mt
import scipy.linalg as la
import sys
from scipy.optimize import fsolve
from Modes_EF_module import Modes_EF
from Chargement_module import Chargement

# #########################################################################
# Dimension de l'aile et propriétés %%% À compléter %%%
Longueur = 10.0
Largeur = 2
Epaisseur = 0.3
Paroi = 0.05
E = 70e9
RHO = 2700
A = (Largeur*Epaisseur) - (Largeur-2*Paroi)*(Epaisseur-2*Paroi)
I = 2*( (1/12)*(Largeur*Paroi**3) + Largeur*Paroi*(Epaisseur/2-Paroi/2)**2 ) + 2*((1/12)*Paroi*(Epaisseur-2*Paroi)**3)
# Moteur
LM =  4 # Emplacement du moteur en m
MM =  1000 # Masse du moteur en kg

# Carburant
LC1 =  2 # Début de la zone de carburant en m
LC2 =  6 # Fin de la zone de carburant en m
MC = 304  # Masse par unité de longueur en kg/m

# #########################################################################
# Maillage
Nb_Elements = 10
L = Longueur/Nb_Elements
DDL = 2*(Nb_Elements+1)

# Matrice de rigidité et de masse d'une poutre %%% À compléter %%%
K_poutre = (E*I/L**3)*np.array([[12, 6*L, -12, 6*L],[6*L, 4*L**2, -6*L, 2*L**2],[-12, -6*L, 12, -6*L],[6*L, 2*L**2, -6*L, 4*L**2]])
M_poutre = (L/420)*np.array([[156, 22*L, 54, -13*L],[22*L, 4*L**2, 13*L, -3*L**2],[54, 13*L, 156, -22*L],[-13*L, -3*L**2, -22*L, 4*L**2]])

# Matrice de rigidité et de masse globale du système %%% À compléter %%%
K = np.zeros((DDL,DDL))
M = np.zeros((DDL,DDL))
for i in range(0,DDL-3,2):
    K[i:i+4,i:i+4] += K_poutre
    M[i:i+4,i:i+4] += RHO*A*M_poutre

# Ajout de la masse du moteur %%% À compléter %%%
Noeud = round(LM/L)
M[2*Noeud,2*Noeud] += MM

# Ajout de la masse du carburant
Noeud1 = round(LC1/L)
Noeud2 = round(LC2/L)
# print(Noeud,Noeud1,Noeud2) %%% À compléter %%%
for i in range(2*Noeud1+2,2*Noeud2-1,2):
     M[i:i+4,i:i+4] += MC*M_poutre 

#Efforts tranchants suivis du moments fléchissants
M1 = M[:2,2:]
K1 = K[:2,2:]

# Conditions aux rives
K = K[2:,2:]
M = M[2:,2:]
Nb_DDL = len(K)
# print(DDL,Nb_DDL)

# #########################################################################
# Calcul des fréquences propres (numériquement) %%% À compléter %%%
LAMBDA, PHI = la.eig(K,M)
# print("PHI est normalisée par la.eig (norme euclidienne)")
OMEGA = sorted(np.sqrt(np.real(LAMBDA)))
print("OMEGA: Les 5 fréquences naturelles les plus basses en rad/s sont")
Mm=np.diag(np.diag(PHI.transpose()@M@PHI))
PHI=PHI@np.sqrt(np.linalg.inv(Mm))
# Verification: Mm doit être la matrice unitaire I
#print("Matrice masse modale")
#Mm = PHI.transpose()@M@PHI
for i in range(Nb_DDL):
    PHI[i,:] = PHI[i,:][np.sqrt(np.real(LAMBDA)).argsort()]

PHI = np.block([[np.zeros((2,Nb_DDL))],[PHI]])

# #########################################################################
# Calcul des modes propres (numériquement)
# Discrétisation de l'aile en 100 points
Nb_Points = 100
x = np.linspace(0,Longueur,Nb_Points)

v_FEM = np.zeros((Nb_Points,DDL))
# Calcul des modes propres
v_FEM = Modes_EF(PHI,x,Longueur,Nb_Elements)

# #########################################################################
# Impression des résultats
# Tracage des 5 premier Modes
Nb_Modes = 5
for i in range(Nb_Modes):
    print(OMEGA[i])

couleur = iter(plt.cm.rainbow(np.linspace(0,0.4,Nb_Modes)))
# use LaTeX fonts in the plot
plt.rc('text', usetex=False)
plt.rc('font', family='serif')

for i in range(Nb_Modes):
    c=next(couleur)
    plt.plot(x[:],v_FEM[:,i],'-',label='mode '+str(i+1), color=c)

plt.xlabel(r'x/L', fontsize=11)
plt.ylabel(r'Amplitude normée', fontsize=11)
plt.legend(loc='best')
plt.grid(True)
plt.savefig('modes_EF.png', dpi = 300, bbox_inches='tight')
plt.show()


# #########################################################################
# Intégration temporelle
# Temps de simulation
#Faire à Temps = 0.1, 0.2 1.0, 3.0
Temps = 3
# Nombre de pas de temps
n = 5000
# visualisation des N_vis derniers pas de temps
N_vis = n-1
# time points
t = np.linspace(0,Temps,n)
# Type de chargement (1 pour Balourd et 2 pour pour Rafale de vent)
Type = 2

# Amortissement
C = 0.8*M

# fonction qui retourne la dérivée temporelle (dz/dt)
def model(t,z,invm,invmc,invmk,n,Type):
    xl = np.zeros(n)
    yl = np.zeros(n)
    dzdt = np.zeros(2*n)
    u = np.zeros(Nb_DDL)
    if Type == 1:
        # Terme de forçage de type Balourd à l'emplacement du moteur
        w = 10000*np.pi/30
        mb = 0.5
        e = 0.3
        u[2*Noeud-2] = mb*e*w**2*np.sin(w*t)
    if Type == 3:
        u[2*Noeud-2] = 1000.
    for i in range(n):
        xl[i] = z[2*i]
        yl[i] = z[2*i+1]
    dxdt = yl
    dydt = - invmc @ yl - invmk @ xl + invm @ u
    for i in range(n):
        dzdt[2*i]   = dxdt[i]
        dzdt[2*i+1] = dydt[i]
    return dzdt

# Calculs préliminaires des matrices du système pour alléger le coût d'évaluation de la solution
invM  = np.linalg.inv(M)
invMC = np.linalg.inv(M)@C
invMK = np.linalg.inv(M)@K
# plage de temps de résolution
tspan = [t[0],t[n-1]]

# initial condition
z0 = np.zeros(2*Nb_DDL)
if Type == 2:
    # Rafale de vent Conditions initiales
    Impulsion = np.zeros((Nb_DDL))
    for i in range(0,Nb_DDL,2):
        Impulsion[i] = 3100
    Impulsion[0] = 3100/2.
    Impulsion[1] = 3100/12.
    Impulsion[Nb_DDL-2] = 3100/2.
    Impulsion[Nb_DDL-1] = -3100/12.
    # Définition des conditions initiales (vitesse non nulle - impulsion)
    ydep = np.zeros((Nb_DDL))
    Impulsion = np.linalg.inv(M)@Impulsion
    for i in range(Nb_DDL):
        z0[2*i+1] = Impulsion[i]

# résolution du système d'EDO couplées
z = solve_ivp(model,tspan,z0,method='Radau',args=(invM,invMC,invMK,Nb_DDL,Type),t_eval=t[n-N_vis-1:n-1])

# enregistrement de la solution pour visualisation
x_sol = np.zeros((N_vis,Nb_DDL))
v_sol = np.zeros((N_vis,Nb_DDL))
a_sol = np.zeros((N_vis,Nb_DDL))
F_sol = np.zeros((N_vis))
M_sol = np.zeros((N_vis))






# DÉPLACEMENT DE L’EXTRÉMITÉ DE L’AILE EN FONCTION DU TEMPS

x_extremite = np.array([Longueur])  
v_extremite = Modes_EF(PHI, x_extremite, Longueur, Nb_Elements) @ z.y[::2, :]  


plt.figure()
plt.plot(z.t, v_extremite[0, :])
plt.xlabel("Temps (s)")
plt.ylabel("Déplacement vertical à x = L (m)")
plt.title("Déplacement de l’extrémité de l’aile dans le temps")
plt.grid(True)
plt.savefig("deplacement_extremite.png", dpi=300, bbox_inches='tight')
plt.show()





k = 0
u = np.zeros(Nb_DDL)
for i in range(n-N_vis-1,n-1):
    if Type == 1:
        # Terme de forçage de type Balourd à l'emplacement du moteur
        w = 10000*np.pi/30
        mb = 0.5
        e = 0.3
    if Type == 3:
        u[2*Noeud-2] = 1000.
    for j in range(Nb_DDL):
        x_sol[k,j] = z.y[2*j,k]
        v_sol[k,j] = z.y[2*j+1,k]
    a_sol[k,:] = - invMC @ v_sol[k,:] - invMK @ x_sol[k,:] + invM @ u[:]
    F_sol[k] = M1[0,:] @ (a_sol[k,:]+0.8*v_sol[k,:]) + K1[0,:] @ x_sol[k,:]
    M_sol[k] = M1[1,:] @ (a_sol[k,:]+0.8*v_sol[k,:]) + K1[1,:] @ x_sol[k,:]
    k += 1

#Affichage du graphique
plt.plot(t[n-N_vis-1:n-1],F_sol[:],'-',label='Effort Tranchant en N', color='r', linewidth=1)
plt.plot(t[n-N_vis-1:n-1],M_sol[:],'-',label='Moment Fléchissant en N/m', color='g', linewidth=1)
plt.xlabel(r'Temps (s)', fontsize=11)
plt.grid(True)    # plt.savefig('solution_temporelle_forcee.png', dpi = 300, bbox_inches='tight')
plt.legend()
plt.show()

# Calcul de la flèche de la poutre en fonction du temps
Chargement(x_sol,Longueur,Nb_Elements)



# ################################################################
# Analyse de convergence 

# Maillage
Nb_Elements1 = 2*Nb_Elements
L1 = Longueur/Nb_Elements1
DDL1 = 2*(Nb_Elements1+1)

# Matrice de rigidité et de masse d'une poutre    %%% A compléter %%%
K_poutre1 = (E*I/L1**3)*np.array([[12, 6*L1, -12, 6*L1],[6*L1,4*L1**2,-6*L1,2*L1**2],[-12,-6*L1,12,-6*L1],[6*L1,2*L1**2,-6*L1,4*L1**2]])
M_poutre1 = (L1/420)*np.array([[156,22*L1,54,-13*L1],[22*L1,4*L1**2,13*L1,-3*L1**2],[54,13*L1,156,-22*L1],[-13*L1,-3*L1**2,-22*L1,4*L1**2]])

# Matrice de rigidité et de masse globale du système    %%% A compléter %%%
K1 = np.zeros((DDL1,DDL1))
M1 = np.zeros((DDL1,DDL1)) #Bloc 4x4 décalé de 2 ddl, on ne doit pas écraser les termes élémentaires, on ajoute par dessus 
for i in range(0,DDL1-3,2):
    K1[i:i+4,i:i+4] = K1[i:i+4,i:i+4] + K_poutre1
    M1[i:i+4,i:i+4] = M1[i:i+4,i:i+4] + RHO*A*M_poutre1

# Ajout de la masse du moteur   %%% A compléter %%%
Noeud3 = round(LM/L1)
M1[2*Noeud3,2*Noeud3] = M1[2*Noeud3,2*Noeud3] + MM

# Ajout de la masse du carburant
Noeud4 = round(LC1/L1)
Noeud5 = round(LC2/L1)
for i in range(2*Noeud4+2,2*Noeud5-1,2):
    M1[i:i+4,i:i+4] = M1[i:i+4,i:i+4] + MC*M_poutre1

#% Conditions aux rives
K1 = K1[2:,2:]
M1 = M1[2:,2:]
Nb_DDL = len(K1)
# print(DDL,Nb_DDL)








# Analyse de convergence 
# a) Fréquences propres 

def analyse_convergence_frequences_EF_plot(epsilon=0.01):
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.linalg import eig

    def assemble_K_M_EF(Nb_Elements):
        L = Longueur / Nb_Elements
        DDL = 2 * (Nb_Elements + 1)
        K_poutre = (E*I/L**3)*np.array([[12, 6*L, -12, 6*L], [6*L, 4*L**2, -6*L, 2*L**2],
                                        [-12, -6*L, 12, -6*L], [6*L, 2*L**2, -6*L, 4*L**2]])
        M_poutre = (L/420)*np.array([[156, 22*L, 54, -13*L], [22*L, 4*L**2, 13*L, -3*L**2],
                                     [54, 13*L, 156, -22*L], [-13*L, -3*L**2, -22*L, 4*L**2]])
        K = np.zeros((DDL, DDL))
        M = np.zeros((DDL, DDL))
        for i in range(0, DDL - 3, 2):
            K[i:i+4, i:i+4] += K_poutre
            M[i:i+4, i:i+4] += RHO * A * M_poutre
        # Masse moteur
        noeud_moteur = round(LM / L)
        M[2*noeud_moteur, 2*noeud_moteur] += MM
        # Masse carburant
        noeud1 = round(LC1 / L)
        noeud2 = round(LC2 / L)
        for i in range(2*noeud1+2, 2*noeud2-1, 2):
            M[i:i+4, i:i+4] += MC * M_poutre
        return K[2:, 2:], M[2:, 2:]

    def get_frequencies(K, M):
        LAMBDA, _ = eig(K, M)
        OMEGA = sorted(np.sqrt(np.real(LAMBDA)))
        return np.array(OMEGA[:5])

    print("\n======= ANALYSE DE CONVERGENCE - FRÉQUENCES PROPRES EF =======")
    N_vals = []
    epsilons = []



    N = 5
    while N <= 320:
        K1, M1 = assemble_K_M_EF(N)
        K2, M2 = assemble_K_M_EF(2*N)
        omega1 = get_frequencies(K1, M1)
        omega2 = get_frequencies(K2, M2)
        
        erreurs_f = np.abs((omega2 - omega1) / omega1)
        epsilon_local = np.max(erreurs_f)
        
        print("\n======= ANALYSE DE CONVERGENCE (Fréquences propres EF) =======")
        print(f"Erreur max = {epsilon_local*100:.6f} % pour N = {N}")
        
        N_vals.append(N)
        epsilons.append(epsilon_local)
        N *= 2


    # Tracé 1 : axes réguliers
    plt.figure(figsize=(8, 5))
    plt.plot(N_vals, [e * 100 for e in epsilons], 'o-', color='orange', linewidth=2, label="Données")
    plt.xlabel("Nombre d'éléments EF (N)")
    plt.ylabel("Erreur maximale (%)")
    plt.title("Erreur maximale vs N (échelle normale)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("Convergence_fréquences_EF.png", dpi=400, bbox_inches='tight')
    
    # Tracé 2 : log-log
    from scipy.optimize import curve_fit
    import numpy as np
    
    def exp_func(x, a, b):
        return a * np.power(x, b)
    
    x_fit = np.array(N_vals)
    y_fit = np.array(epsilons) * 100  # pour affichage en %
    
    params, _ = curve_fit(exp_func, x_fit, y_fit)
    a, b = params
    y_model = exp_func(x_fit, a, b)
    
    plt.figure(figsize=(8, 5))
    plt.loglog(N_vals, y_fit, 'o-', color='orange', linewidth=2, label='Données')
    plt.loglog(N_vals, y_model, 'r--', label=fr"$y = {a:.2f} \cdot x^{{{b:.2f}}}$")
    plt.xlabel("N (log)")
    plt.ylabel("Erreur maximale (%)")
    plt.title("Erreur maximale vs N (log-log)")
    plt.grid(True, which="both", ls="--")
    plt.legend()
    plt.tight_layout()
    plt.savefig("Convergence_fréquences_log_EF.png", dpi=400, bbox_inches='tight')


analyse_convergence_frequences_EF_plot(epsilon=0.01)


# b) Modes propres


from scipy.linalg import eig
from numpy import linspace, sum, zeros, block
import matplotlib.pyplot as plt




def assemble_K_M_EF(Nb_Elements):
    L = Longueur / Nb_Elements
    DDL = 2 * (Nb_Elements + 1)
    K_poutre = (E*I/L**3)*np.array([[12,6*L,-12,6*L],[6*L,4*L**2,-6*L,2*L**2],[-12,-6*L,12,-6*L],[6*L,2*L**2,-6*L,4*L**2]])
    M_poutre = (L/420)*np.array([[156,22*L,54,-13*L],[22*L,4*L**2,13*L,-3*L**2],[54,13*L,156,-22*L],[-13*L,-3*L**2,-22*L,4*L**2]])
    K = np.zeros((DDL, DDL))
    M = np.zeros((DDL, DDL))
    for i in range(0, DDL - 3, 2):
        K[i:i+4, i:i+4] += K_poutre
        M[i:i+4, i:i+4] += RHO * A * M_poutre
    # masses ajoutées
    noeud_moteur = round(LM / L)
    M[2*noeud_moteur, 2*noeud_moteur] += MM
    noeud1 = round(LC1 / L)
    noeud2 = round(LC2 / L)
    for i in range(2*noeud1+2, 2*noeud2-1, 2):
        M[i:i+4, i:i+4] += MC * M_poutre
    return K[2:, 2:], M[2:, 2:]

# Récupérer les modes propres
from Modes_EF_module import Modes_EF

def get_modes(K, M, N_elements):
    _, PHI = eig(K, M)
    idx = np.argsort(np.sqrt(np.real(_)))
    PHI = PHI[:, idx]
    PHI = block([[zeros((2, PHI.shape[1]))], [PHI]])
    x = linspace(0, Longueur, 100)
    v = Modes_EF(PHI, x, Longueur, N_elements)
    return v[:, :5]


# Lancer l'analyse sur différentes tailles de maillage
N_vals = []
delta_vals = []

N = 5
while N <= 320:
    vN = get_modes(*assemble_K_M_EF(N), N)
    v2N = get_modes(*assemble_K_M_EF(2*N), 2*N)
    erreurs_modes = []

    for i in range(5):
        num = sum((v2N[:, i] - vN[:, i])**2)
        denom = sum((vN[:, i])**2)
        erreur_i = num / denom
        erreurs_modes.append(erreur_i)

    delta = max(erreurs_modes)

    print("\n======= ANALYSE DE CONVERGENCE (Modes propres EF) =======")
    print(f"Erreur max = {delta*100:.6f} % pour N = {N}")

    N_vals.append(N)
    delta_vals.append(delta)
    N *= 2

    

# Tracé 1 : axes réguliers
plt.figure(figsize=(8, 5))
plt.plot(N_vals, [d * 100 for d in delta_vals], 'o-', color='orange', linewidth=2, label="Données")
plt.xlabel("Nombre d'éléments EF (N)")
plt.ylabel("Erreur quadratique maximale (%)")
plt.title("Erreur maximale vs N (échelle normale)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("Convergence_modes_EF.png", dpi=400, bbox_inches='tight')

# Tracé 2 : axes log-log avec régression exponentielle
from scipy.optimize import curve_fit
import numpy as np

def exp_func(x, a, b):
    return a * np.power(x, b)

x_fit = np.array(N_vals)
y_fit = np.array(delta_vals) * 100
params, _ = curve_fit(exp_func, x_fit, y_fit)
a, b = params
y_model = exp_func(x_fit, a, b)

plt.figure(figsize=(8, 5))
plt.loglog(N_vals, y_fit, 'o-', color='orange', linewidth=2, label='Données')
plt.loglog(N_vals, y_model, 'r--', label=fr"$y = {a:.2f} \cdot x^{{{b:.2f}}}$")
plt.xlabel("N (log)")
plt.ylabel("Erreur quadratique maximale (%)")
plt.title("Erreur maximale vs N (log-log)")
plt.grid(True, which="both", linestyle="--")
plt.legend()
plt.tight_layout()
plt.savefig("Convergence_modes_log_EF.png", dpi=400, bbox_inches='tight')











