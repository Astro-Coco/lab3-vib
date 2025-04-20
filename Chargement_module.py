def Chargement(Deplacement,Longueur,Nb_Elements):
    import numpy as np
    import matplotlib.pyplot as plt
    #  Description
    # 
    #  Cette fonction permet de calculer la flèche de l'aile en tout point x
    #  suite au chargement imposé (Balourd ou Rafale de vent)
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #  Variables d'entrée
    # 
    #  Deplacement: Matrice de déplacement obtenu par Simulink
    #  Longueur: Longueur de l'aile
    #  Nb_Elements: Nombre d'éléments dans le modèle EF
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #  Variables de sortie
    #  v_chargement: Matrice [Nb_Points x Nb_Temps]. Chaque colonne de la
    #  matrice correspond à la flèche de la poutre à un instant donné.
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # Ajouter des 0 à la matrice Déplacement pour tenir compte des DDLs à
    # l'encastrement
    Nb_Temps = np.shape(Deplacement)[0]
    # print(np.shape(Deplacement))
    Deplacement = np.block([np.zeros((Nb_Temps,2)),Deplacement])
    # print(np.shape(Deplacement))
    #
    # Calcul de la longueur de chaque élément et du nombre de DDL actifs dans
    # le modèle EF
    L = Longueur/Nb_Elements

    # Discrétisation de l'aile en 100 points
    Nb_Points = 100
    x = np.linspace(0,Longueur,Nb_Points)

    # Initialisation de la matrice v_chargement
    v_chargement = np.zeros((Nb_Points,Nb_Temps))

    # Identification du numéro d'élément auquel appartient le ième point de l'aile
    Numero_Element = np.ceil(x/L)
    Numero_Element[0] = 1

    # print(np.shape(v_chargement))
    for i in range(Nb_Points):
        # Calcul de x local dans l'élément poutre
        x_local = x[i] - (Numero_Element[i]-1)*L
        # %Identification des 4 DDLs
        j = int(2*Numero_Element[i]-2)
        # %Fonctions de forme
        N1 = (2/L**3)*(x_local**3) - (3/L**2)*(x_local**2) + 1
        N2 = (1/L**2)*(x_local**3) - (2/L)*(x_local**2) + x_local
        N3 = (-2/L**3)*(x_local**3) + (3/L**2)*(x_local**2)
        N4 = (1/L**2)*(x_local**3) - (1/L)*(x_local**2)
        # %Calcul de la flèche en un point donné, pour tous les temps simulés
        v_chargement[i,:] = N1*Deplacement[:,j] + N2*Deplacement[:,j+1] + N3*Deplacement[:,j+2] + N4*Deplacement[:,j+3]

    # use LaTeX fonts in the plot
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    # Tracage des déformées de l'aile sur les n pas de temps
    couleur = iter(plt.cm.rainbow(np.linspace(0,0.4,Nb_Temps)))
    for i in range(Nb_Temps):
        c=next(couleur)
        plt.plot(x[:],v_chargement[:,i],'-',label='mode '+str(i+1), color=c, linewidth=0.25)

    plt.xlabel(r'x/L', fontsize=11)
    plt.ylabel(r'Amplitude normée', fontsize=11)
    # plt.legend(loc='best')
    plt.grid(True)
    # plt.savefig('solution_temporelle_forcee.png', dpi = 300, bbox_inches='tight')
    plt.savefig('solution_temporelle_libre.png', dpi = 300, bbox_inches='tight')
    plt.show()

    return