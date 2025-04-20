def Modes_EF(PHI,x,Longueur,Nb_Elements):
    import numpy as np
    # Description
    #
    # Cette fonction permet de calculer les modes propres de votre modèle à
    # partir de la matrice des vecteurs propres PHI. Les fonctions de formes
    # Ni(x) sont utilisées.
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # % Variables d'entrée
    # %
    # % PHI: Matrice des vecteurs propres obtenus par EF
    # % x: Coordonnées des points de discrétisation de l'aile
    # % Longueur: Longueur de l'aile
    # % Nb_Elements: Nombre d'éléments dans le modèle EF
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Variables de sortie
    #
    # v_FEM: Matrice [Nb_Points x Nb_DDL]. Chaque colonne de la matrice
    # correspond à un mode propre en fonction de x.
    #
    # Calcul de la longueur de chaque élément et du nombre de DDL actifs dans
    # le modèle EF
    L = Longueur/Nb_Elements
    Nb_DDL = np.shape(PHI)[1]
    # Calcul du nombre de points de discrétisation
    Nb_Points = len(x)
    # Initialisation de la matrice v_FEM
    v_FEM = np.zeros((Nb_Points,Nb_DDL))
    # Identification du numéro d'élément auquel appartient le ième point de l'aile
    Numero_Element = np.ceil(x/L)
    Numero_Element[0] = 1
    for i in range(Nb_Points):
        # Calcul de x local dans l'élément poutre
        x_local = x[i] - (Numero_Element[i]-1)*L
        # Identification des 4 DDLs
        j = np.int(2*Numero_Element[i]-2)
        # Fonctions de forme %%% À compléter %%%
        N1 =  
        N2 =  
        N3 =  
        N4 =  
        # Calcul de la flèche à la position i, et ce pour chaque mode propre. %%% À compléter %%%
        v_FEM[i,:] =  
    # Normalisation des vecteurs propres
    for i in range(Nb_DDL):
        v_FEM[:,i] = abs(v_FEM[1,i])/(v_FEM[1,i])*(v_FEM[:,i]/np.linalg.norm(v_FEM[:,i]))
    return v_FEM