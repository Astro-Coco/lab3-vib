
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, curve_fit
import scipy.integrate as integrate
import scipy.linalg as la
from Xi_module import Xi, Xij, Xij_dd

# === Paramètres physiques ===
Longueur = 10
Largeur = 2
Epaisseur = 0.3
Paroi = 0.05
E = 70e9
RHO = 2700
MC = 304.
MM = 1000.
LM = 4
LC1 = 2
LC2 = 6
A = (Largeur * Epaisseur) - (Largeur - 2*Paroi)*(Epaisseur - 2*Paroi)
I = 2*((1/12)*(Largeur*Paroi**3) + Largeur*Paroi*(Epaisseur/2 - Paroi/2)**2) + 2*((1/12)*Paroi*(Epaisseur-2*Paroi)**3)
m = RHO * A
x = np.linspace(0, Longueur, 100)

# === A) Étude de convergence des fréquences propres ===
def calcul_omega(N):
    BL = np.array([fsolve(lambda b: np.cos(b) + 1/np.cosh(b), [(i + 0.5)*np.pi])[0] for i in range(N)])
    M = np.zeros((N, N))
    K = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            M[i, j] = m * integrate.quad(Xij, 0, Longueur, args=(BL[i], BL[j], Longueur))[0]
            M[i, j] += MC * integrate.quad(Xij, LC1, LC2, args=(BL[i], BL[j], Longueur))[0]
            M[i, j] += MM * Xij(LM, BL[i], BL[j], Longueur)
            K[i, j] = E * I * integrate.quad(Xij_dd, 0, Longueur, args=(BL[i], BL[j], Longueur))[0]
    LAMBDA, _ = la.eig(K, M)
    omega = np.sort(np.sqrt(np.real(LAMBDA)))
    return omega[:5]

N_vals = list(range(5, 16))
erreurs_freq = []
print("\n📈 Erreurs relatives maximales sur les fréquences propres:")
for N in N_vals:
    omega_N = calcul_omega(N)
    omega_2N = calcul_omega(2*N)
    eps = np.abs((omega_2N - omega_N) / omega_N)
    err_max = np.max(eps) * 100
    erreurs_freq.append(err_max)
    print(f"N = {N}  →  erreur max = {err_max:.4f} %")

# === Fit exponentiel (fréquences) ===
def expo(x, a, b):
    return a * np.exp(b * x)

N_array = np.array(N_vals)
err_freq_array = np.array(erreurs_freq)
params_freq, _ = curve_fit(expo, N_array, err_freq_array)
a_freq, b_freq = params_freq
xfit = np.linspace(min(N_array), max(N_array), 100)
yfit = expo(xfit, a_freq, b_freq)

# === B) Étude de convergence des modes propres ===
def calcul_modes(N):
    BL = np.array([fsolve(lambda b: np.cos(b) + 1/np.cosh(b), [(i + 0.5)*np.pi])[0] for i in range(N)])
    M = np.zeros((N, N))
    K = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            M[i, j] = m * integrate.quad(Xij, 0, Longueur, args=(BL[i], BL[j], Longueur))[0]
            M[i, j] += MC * integrate.quad(Xij, LC1, LC2, args=(BL[i], BL[j], Longueur))[0]
            M[i, j] += MM * Xij(LM, BL[i], BL[j], Longueur)
            K[i, j] = E * I * integrate.quad(Xij_dd, 0, Longueur, args=(BL[i], BL[j], Longueur))[0]
    _, PHI = la.eig(K, M)
    idx = np.argsort(np.real(np.diag(PHI.T @ K @ PHI)))
    PHI = PHI[:, idx]
    v = np.zeros((len(x), N))
    for i in range(N):
        for n in range(N):
            v[:, i] += PHI[n, i] * Xi(x, BL[n], Longueur)
        v[:, i] *= np.sign(v[2, i]) / np.linalg.norm(v[:, i])
    return v[:, :5]

print("\n📊 Erreurs relatives maximales sur les modes propres:")
err_modes = []
for N in N_vals:
    vN = calcul_modes(N)
    v2N = calcul_modes(2 * N)
    deltas = [np.sum((v2N[:, i] - vN[:, i])**2) / np.sum(vN[:, i]**2) for i in range(5)]
    delta_max = np.max(deltas) * 100
    err_modes.append(delta_max)
    print(f"N = {N}  →  erreur max = {delta_max:.4f} %")

# === Fit exponentiel (modes) ===
err_modes_array = np.array(err_modes)
params_modes, _ = curve_fit(expo, N_array, err_modes_array)
a_modes, b_modes = params_modes
yfit_modes = expo(xfit, a_modes, b_modes)

# === Affichage des résultats ===
plt.figure(figsize=(14, 6))

plt.subplot(1, 2, 1)
plt.plot(N_array, err_freq_array, 'o-', color='royalblue', label="Fréquences")
plt.plot(xfit, yfit, 'r--', label=fr"$y = {a_freq:.2f} e^{{{b_freq:.2f}N}}$")
plt.xlabel("Nombre de modes (N)")
plt.ylabel("Erreur max fréquence (%)")
plt.title("Convergence des fréquences naturelles")
plt.grid(True)
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(N_array, err_modes_array, 'o-', color='darkorange', label="Modes propres")
plt.plot(xfit, yfit_modes, 'r--', label=fr"$y = {a_modes:.2f} e^{{{b_modes:.2f}N}}$")
plt.xlabel("Nombre de modes (N)")
plt.ylabel("Erreur max forme modale (%)")
plt.title("Convergence des formes modales")
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.savefig("convergence_frequences_et_modes.png", dpi=300)
plt.show()
