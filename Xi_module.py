def Xi(x,BiL,Longueur, coeff = False):
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
    x = np.array(x)
    
    Beta = BiL / Longueur
    denominateur = 1 + 2*np.sin(BiL)*np.exp(-BiL)-np.exp(-2*BiL)

    ci = ((np.cos(BiL) + np.sin(BiL)) * np.exp(-BiL) + 1) / denominateur
    di = (-np.cos(BiL) + np.sin(BiL) - np.exp(-BiL)) / denominateur
    ai = ci - di * np.exp(-BiL)

    if coeff:
        return ci, di, ai
    else:
        xi = ai * np.sin(Beta * x) - np.cos(Beta * x) + ci * np.exp(-Beta * x) + di * np.exp(Beta * (x - Longueur))
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

    xi = Xi(x, BiL, Longueur)
    xj = Xi(x, BjL, Longueur)
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
    ci , di, ai = Xi(x, BiL, Longueur, coeff = True)
    beta_i = BiL / Longueur
    xi_dd = -ai * beta_i**2 * np.sin(beta_i * x) \
            - beta_i**2 * np.cos(beta_i * x) \
            + ci * beta_i**2 * np.exp(-beta_i * x) \
            + di * beta_i**2 * np.exp(beta_i * (x - Longueur))

    cj, dj, aj = Xi(x, BjL, Longueur, coeff = True)
    beta_j = BjL / Longueur


    xj_dd = -aj * beta_j**2 * np.sin(beta_j * x) \
            - beta_j**2 * np.cos(beta_j * x) \
            + cj * beta_j**2 * np.exp(-beta_j * x) \
            + dj * beta_j**2 * np.exp(beta_j * (x - Longueur))

    xij_dd = xi_dd * xj_dd
    return xij_dd
