def Xi(x,BiL,Longueur):
    import numpy as np
    #  Description
    # 
    #  Cette fonction permet de calculer le ième mode propre théorique de
    #  l'aile (sans moteur ni carburant)
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #  Variables d'entrée
    # 
    #  x: Coordonnées des points de discrétisation de l'aile
    #  BiL: ième valeur de betak L
    #  Longueur: Longueur de l'aile
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #  Variables de sortie
    #
    #  xi: Mode propre en fonction de la position x.
    #
    # (%%% A compléter %%%)
    ci = 
    di = 
    ai = 

    xi = ai * np.sin(beta * x) - np.cos(beta * x) + ci * np.exp(-beta * x) + di * np.exp(beta * (x - L))
    ... \
    ... \
    ...

    return xi
    
def Xij(x,BiL,BjL,Longueur):
    #  Description
    # 
    #  Cette fonction permet de calculer le produit phii(x)*phij(x)
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # %% Variables d'entrée
    # %    #
    # % x: Coordonnées des points de discrétisation de l'aile
    # % BiL: ième valeur de betak L
    # % BjL: jème valeur de betak L
    # % Longueur: Longueur de l'aile
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #  Variables de sortie
    #
    #  Xij: Valeur du produit Xi(x)*Xj(x) en fonction de x
    #
    # (%%% A compléter %%%)

    xi = # appel à Xi pour BiL...
    xj = # appel à Xi pour BjL...
    xij = xi * xj # produite terme à terme
    return xij

def Xij_dd(x,BiL,BjL,Longueur):
    import numpy as np
    # Description
    # 
    #  Cette fonction permet de calculer le produit phii''(x)*phij''(x)
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # %% Variables d'entrée
    # %
    # % x: Coordonnées des points de discrétisation de l'aile
    # % BiL: ième valeur de betak L
    # % BjL: jème valeur de betak L
    # % Longueur: Longueur de l'aile
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #  Variables de sortie
    #
    #  xij_dd: Valeur du produit xi''(x)*xj''(x) en fonction de x
    #
    # (%%% A compléter %%%)
    ci = 
    di = 
    ai = 

    xi_dd = ... \
    ... \
    ... \
    ...

    cj = 
    dj = 
    aj = 

    xj_dd = ... \
    ... \
    ... \
    ...

    xij_dd = xi_dd * xj_dd
    return xij_dd
