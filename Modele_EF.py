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
E = 
RHO = 
A = 
I = 
# Moteur
LM =  # Emplacement du moteur en m
MM =  # Masse du moteur en kg

# Carburant
LC1 =   # Début de la zone de carburant en m
LC2 =   # Fin de la zone de carburant en m
MC =   # Masse par unité de longueur en kg/m

# #########################################################################
# Maillage
Nb_Elements = 160
L = Longueur/Nb_Elements
DDL = 2*(Nb_Elements+1)

# Matrice de rigidité et de masse d'une poutre %%% À compléter %%%
K_poutre = 
M_poutre = 

# Matrice de rigidité et de masse globale du système %%% À compléter %%%
K = np.zeros((DDL,DDL))
M = np.zeros((DDL,DDL))
for i in range(0,DDL-3,2):
    K[i:i+4,i:i+4] = 
    M[i:i+4,i:i+4] =  

# Ajout de la masse du moteur %%% À compléter %%%
Noeud = round(LM/L)
M[2*Noeud,2*Noeud] =  

# Ajout de la masse du carburant
Noeud1 = round(LC1/L)
Noeud2 = round(LC2/L)
# print(Noeud,Noeud1,Noeud2) %%% À compléter %%%
for i in range(2*Noeud1+2,2*Noeud2-1,2):
     M[i:i+4,i:i+4] =  

# Conditions aux rives
K = K[2:,2:]
M = M[2:,2:]
Nb_DDL = len(K)
# print(DDL,Nb_DDL)

# #########################################################################
# Calcul des fréquences propres (numériquement) %%% À compléter %%%
LAMBDA, PHI =  
# print("PHI est normalisée par la.eig (norme euclidienne)")
OMEGA = sorted(np.sqrt(np.real(LAMBDA)))
# print("OMEGA: Les 5 fréquences naturelles les plus basses en rad/s sont")
Mm=np.diag(np.diag(PHI.transpose()@M@PHI))
PHI=PHI@np.sqrt(np.linalg.inv(Mm))
# Verification: Mm doit être la matrice unitaire I
# print("Matrice masse modale")
# Mm = PHI.transpose()@M@PHI
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

# # #########################################################################
# # Impression des résultats
# # Tracage des 5 premier Modes
# Nb_Modes = 5
# for i in range(Nb_Modes):
#     print(OMEGA[i])
#
# couleur = iter(plt.cm.rainbow(np.linspace(0,0.4,Nb_Modes)))
# # use LaTeX fonts in the plot
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')
#
# for i in range(Nb_Modes):
#     c=next(couleur)
#     plt.plot(x[:],v_FEM[:,i],'-',label='mode '+str(i+1), color=c)
#
# plt.xlabel(r'x/L', fontsize=11)
# plt.ylabel(r'Amplitude normée', fontsize=11)
# plt.legend(loc='best')
# plt.grid(True)
# plt.savefig('modes_EF.png', dpi = 300, bbox_inches='tight')
# plt.show()
#
#
# # #########################################################################
# # Intégration temporelle
# # Temps de simulation
# Temps = 4.
# # Nombre de pas de temps
# n = 5000
# # time points
# t = np.linspace(0,Temps,n)
# # Type de chargement (1 pour Balourd et 2 pour pour Rafale de vent)
# Type = 2
#
# # Amortissement
# C = 0.8*M
#
# # fonction qui retourne la dérivée temporelle (dz/dt)
# def model(t,z,invm,invmc,invmk,n,Type):
#     xl = np.zeros(n)
#     yl = np.zeros(n)
#     dzdt = np.zeros(2*n)
#     u = np.zeros(Nb_DDL)
#     if Type == 1:
#         # Terme de forçage de type Balourd à l'emplacement du moteur
#         w = 10000*np.pi/30
#         mb = 0.5
#         e = 0.3
#         u[2*Noeud-2] = mb*e*w**2*np.sin(w*t)
#     for i in range(n):
#         xl[i] = z[2*i]
#         yl[i] = z[2*i+1]
#     dxdt = yl
#     dydt = - invmc @ yl - invmk @ xl + invm @ u
#     for i in range(n):
#         dzdt[2*i]   = dxdt[i]
#         dzdt[2*i+1] = dydt[i]
#     return dzdt
#
# # Calculs préliminaires des matrices du système pour alléger le coût d'évaluation de la solution
# invM  = np.linalg.inv(M)
# invMC = np.linalg.inv(M)@C
# invMK = np.linalg.inv(M)@K
# # plage de temps de résolution
# tspan = [t[0],t[n-1]]
# # visualisation des N_vis derniers pas de temps
# N_vis = 1000
#
# # initial condition
# z0 = np.zeros(2*Nb_DDL)
# if Type == 2:
#     # Rafale de vent Conditions initiales
#     Impulsion = np.zeros((Nb_DDL))
#     for i in range(0,Nb_DDL,2):
#         Impulsion[i] = 3100
#     Impulsion[0] = 3100/2.
#     Impulsion[1] = 3100/12.
#     Impulsion[Nb_DDL-2] = 3100/2.
#     Impulsion[Nb_DDL-1] = -3100/12.
#     # Définition des conditions initiales (vitesse non nulle - impulsion)
#     ydep = np.zeros((Nb_DDL))
#     Impulsion = np.linalg.inv(M)@Impulsion
#     for i in range(Nb_DDL):
#         z0[2*i+1] = Impulsion[i]
#
# # résolution du système d'EDO couplées
# z = solve_ivp(model,tspan,z0,method='Radau',args=(invM,invMC,invMK,Nb_DDL,Type),t_eval=t[n-N_vis-1:n-1])
#
# # enregistrement de la solution pour visualisation
# x_sol = np.zeros((N_vis,Nb_DDL))
# v_sol = np.zeros((N_vis,Nb_DDL))
#
# k = 0
# for i in range(n-N_vis-1,n-1):
#     for j in range(Nb_DDL):
#         x_sol[k,j] = z.y[2*j,k]
#         v_sol[k,j] = z.y[2*j+1,k]
#     k += 1
#
# # Calcul de la flèche de la poutre en fonction du temps
# Chargement(x_sol,Longueur,Nb_Elements)

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

# #########################################################################
# Calcul des fréquences propres (numériquement)
LAMBDA1, PHI1 = la.eig(K1,M1)
# print("PHI est normalisée par la.eig (norme euclidienne)")
OMEGA1 = sorted(np.sqrt(np.real(LAMBDA1)))
# print("OMEGA: Les 5 fréquences naturelles les plus basses en rad/s sont")
Mm1=np.diag(np.diag(PHI1.transpose()@M1@PHI1))
PHI1=PHI1@np.sqrt(np.linalg.inv(Mm1))
# Verification: Mm doit être la matrice unitaire I
# print("Matrice masse modale")
# Mm = PHI.transpose()@M@PHI
for i in range(Nb_DDL):
    PHI1[i,:] = PHI1[i,:][np.sqrt(np.real(LAMBDA1)).argsort()]

PHI1 = np.block([[np.zeros((2,Nb_DDL))],[PHI1]])

print(OMEGA[0:5])
print(OMEGA1[0:5])

# eps = np.max(abs((OMEGA1[0:4]-OMEGA[0:4])/OMEGA[0:4]))
eps =np.zeros(5)
for i in range(5):
    eps[i] = abs((OMEGA1[i]-OMEGA[i])/OMEGA[i])
print('epsilon - convergence fréquences')
print(max(eps))

# #########################################################################
# Calcul des modes propres (numériquement)
# Discrétisation de l'aile en 100 points
Nb_Points = 100
x = np.linspace(0,Longueur,Nb_Points)

v_FEM1 = np.zeros((Nb_Points,DDL))
# Calcul des modes propres
v_FEM1 = Modes_EF(PHI1,x,Longueur,Nb_Elements1)
#Analyse de convergence des modes propres

delta = 0.
for i in range(5):
    delta = delta + np.max(sum((v_FEM1[1:100,i]-v_FEM[1:100,i])**2)/sum((v_FEM[1:100,i])**2))
print('delta - convergence vecteurs propres')
print(delta)
