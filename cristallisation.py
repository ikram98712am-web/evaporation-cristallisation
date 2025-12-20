"""
cristallisation_complete.py
---------------------------
Modèle complet de cristallisation batch du saccharose
avec toutes les fonctionnalités du projet.
"""

import numpy as np
from scipy.integrate import odeint, solve_ivp
from scipy.optimize import minimize, fsolve
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import warnings
warnings.filterwarnings('ignore')

# ============================================
# 1. CONSTANTES ET PARAMÈTRES
# ============================================

R = 8.31446261815324  # J/mol/K

# Constantes cinétiques (PDF)
kb = 1.5e10           # noyaux/(m³·s)
b_nuc = 2.5           # exposant de sursaturation (nucléation)
j_nuc = 0.5           # exposant de la masse (nucléation)

kg = 2.8e-7           # m/s (pré-exponentiel croissance)
g_gr = 1.5            # exposant de sursaturation (croissance)
Eg = 45e3             # J/mol (énergie d'activation croissance)

# Propriétés physiques
rho_cristal = 1580.0  # kg/m³ (masse volumique saccharose)
rho_solution = 1300.0  # kg/m³ (solution sucrée approx)
k_v = np.pi / 6.0     # facteur de forme volume (sphère)
k_a = np.pi          # facteur de forme surface (sphère)

# Conditions opératoires
T0 = 70.0 + 273.15    # Température initiale [K]
Tf = 35.0 + 273.15    # Température finale [K]
t_total = 4 * 3600.0  # Temps total [s]

# ============================================
# 2. FONCTIONS THERMODYNAMIQUES
# ============================================

def solubilite_saccharose(T):
    """
    Solubilité C* (g saccharose / 100 g solution) en fonction de T (K).
    
    Formule du sujet :
    C* = 64.18 + 0.1337 T + 5.52e-3 T^2 - 9.73e-6 T^3
    avec T en °C.
    
    Returns
    -------
    float
        Solubilité [g/100g solution]
    """
    T_C = T - 273.15  # conversion K -> °C
    C_star = (
        64.18
        + 0.1337 * T_C
        + 5.52e-3 * T_C**2
        - 9.73e-6 * T_C**3
    )
    return max(C_star, 0.0)


def concentration_to_masse(C):
    """
    Convertit la concentration [g/100g] en fraction massique.
    """
    return C / 100.0


def sursaturation(C, T):
    """
    Sursaturation relative S = (C - C*) / C*.
    
    Returns
    -------
    float
        Sursaturation [-]
    """
    C_star = solubilite_saccharose(T)
    if C_star <= 1e-8:
        return 0.0
    return max((C - C_star) / C_star, 0.0)


def concentration_from_sursaturation(S, T):
    """
    Calcule la concentration à partir de la sursaturation.
    """
    C_star = solubilite_saccharose(T)
    return C_star * (1.0 + S)


# ============================================
# 3. CINÉTIQUES DE CRISTALLISATION
# ============================================

def taux_nucleation(S, m_solide=0.0):
    """
    Taux de nucléation B = kb S^b m^j.
    
    Parameters
    ----------
    S : float
        Sursaturation [-]
    m_solide : float, optional
        Fraction massique de solide [-]
        
    Returns
    -------
    float
        Taux de nucléation [noyaux/(m³·s)]
    """
    if S <= 0:
        return 0.0
    return kb * (S ** b_nuc) * (max(m_solide, 1e-6) ** j_nuc)


def taux_croissance(S, T):
    """
    Taux de croissance G = kg S^g exp(-Eg / (R T)).
    
    Returns
    -------
    float
        Taux de croissance [m/s]
    """
    if S <= 0:
        return 0.0
    return kg * (S ** g_gr) * np.exp(-Eg / (R * T))


# ============================================
# 4. PROFILS DE REFROIDISSEMENT
# ============================================

def profil_refroidissement_lineaire(t, T0, Tf, t_total):
    """
    Profil linéaire : T(t) = T0 - αt
    """
    alpha = (T0 - Tf) / t_total
    return T0 - alpha * t


def profil_refroidissement_exponentiel(t, T0, Tf, t_total, beta=None):
    """
    Profil exponentiel : T(t) = Tf + (T0 - Tf) exp(-βt)
    """
    if beta is None:
        beta = 3.0 / t_total  # Constante pour arriver à ~95% à t_total
    return Tf + (T0 - Tf) * np.exp(-beta * t)


def profil_refroidissement_optimal(t, T0, S_cible=0.05, C0=65.0):
    """
    Profil optimal pour maintenir S constant.
    Résolution numérique de l'équation différentielle.
    """
    # Fonction pour trouver T qui donne S = S_cible à chaque instant
    def equation_temperature(T, C_current):
        C_star = solubilite_saccharose(T)
        S_current = (C_current - C_star) / C_star if C_star > 0 else 0
        return S_current - S_cible
    
    # Approche simplifiée : on génère un profil pré-calculé
    # En réalité, il faudrait résoudre simultanément avec les bilans
    times_opt = np.linspace(0, t_total, 100)
    T_opt = np.zeros_like(times_opt)
    T_opt[0] = T0
    
    # Approximation : refroidissement plus lent au début
    for i in range(1, len(times_opt)):
        # Diminution progressive de la pente
        fraction = times_opt[i] / t_total
        pente = 1.0 - 0.7 * fraction  # Diminue la pente avec le temps
        T_opt[i] = T0 + (Tf - T0) * (fraction ** pente)
    
    # Interpolation pour tout t
    return np.interp(t, times_opt, T_opt)


# ============================================
# 5. BILAN DE POPULATION (MÉTHODE DES CLASSES)
# ============================================

class CristalliseurBatch:
    """
    Modèle de cristalliseur batch avec bilan de population.
    """
    
    def __init__(self, n_classes=100, L_min=1e-6, L_max=2000e-6):
        """
        Initialise le cristalliseur.
        
        Parameters
        ----------
        n_classes : int
            Nombre de classes de taille
        L_min : float
            Taille minimale [m]
        L_max : float
            Taille maximale [m]
        """
        self.n_classes = n_classes
        self.L_min = L_min
        self.L_max = L_max
        
        # Classes de taille (centres)
        self.dL = (L_max - L_min) / n_classes
        self.L = np.linspace(L_min + self.dL/2, L_max - self.dL/2, n_classes)
        
        # Variables d'état
        self.C = None      # Concentration [g/100g]
        self.n = None      # Distribution de taille [#/m³/m]
        self.T = None      # Température [K]
        self.t = 0.0       # Temps [s]
        
        # Historique
        self.history = {
            'time': [],
            'C': [],
            'T': [],
            'S': [],
            'B': [],
            'G': [],
            'moments': [],
            'n_distrib': []
        }
    
    def set_conditions_initiales(self, C0, T0):
        """
        Définit les conditions initiales.
        """
        self.C = C0
        self.T = T0
        self.n = np.zeros(self.n_classes)  # Pas de cristaux initialement
        self.t = 0.0
        
        # Réinitialise l'historique
        for key in self.history:
            self.history[key] = []
        
        self._enregistrer_etat()
    
    def _calculer_moments(self, n):
        """
        Calcule les moments de la distribution.
        
        Returns
        -------
        dict
            Moments m0, m1, m2, m3
        """
        m0 = np.sum(n) * self.dL                     # Moment 0 [#/m³]
        m1 = np.sum(n * self.L) * self.dL            # Moment 1 [m/m³]
        m2 = np.sum(n * self.L**2) * self.dL         # Moment 2 [m²/m³]
        m3 = np.sum(n * self.L**3) * self.dL         # Moment 3 [m³/m³]
        
        # Masse de cristaux
        masse_cristaux = k_v * rho_cristal * m3      # [kg/m³]
        
        return {
            'm0': m0,
            'm1': m1,
            'm2': m2,
            'm3': m3,
            'masse': masse_cristaux
        }
    
    def _calculer_fraction_massique(self):
        """
        Calcule la fraction massique de cristaux.
        """
        moments = self._calculer_moments(self.n)
        masse_cristaux = moments['masse']
        
        # Masse totale de saccharose dans 1 m³ de solution
        C_frac = self.C / 100.0  # Conversion g/100g -> fraction
        masse_saccharose_totale = rho_solution * C_frac  # [kg/m³]
        
        if masse_saccharose_totale > 0:
            return masse_cristaux / masse_saccharose_totale
        else:
            return 0.0
    
    def _calculer_derivees(self, t, y, T_func):
        """
        Calcule les dérivées pour l'intégration.
        """
        # Extraction des variables
        C = y[0]
        n = y[1:self.n_classes+1]
        
        # Température à l'instant t
        T = T_func(t)
        
        # Calcul des grandeurs cinétiques
        m_solide = self._calculer_fraction_massique()
        S = sursaturation(C, T)
        B = taux_nucleation(S, m_solide)
        G = taux_croissance(S, T)
        
        # Dérivée de la concentration (bilan matière)
        moments = self._calculer_moments(n)
        surface_totale = k_a * moments['m2']  # Surface totale des cristaux [m²/m³]
        
        # Masse cristallisée par unité de temps
        dm_dt = 3.0 * k_v * rho_cristal * G * moments['m2']  # [kg/(m³·s)]
        
        # Conversion en dérivée de concentration
        # dm_dt = - (masse solution) * dC/dt
        # Mais attention aux unités : C en g/100g
        dC_dt = -dm_dt / rho_solution * 100.0  # [g/100g·s]
        
        # Dérivée de la distribution (bilan population)
        dn_dt = np.zeros_like(n)
        
        # Terme de croissance (advection)
        if G > 0:
            # Schéma upwind pour la stabilité
            for i in range(1, self.n_classes):
                flux = G * n[i-1]
                dn_dt[i] += flux / self.dL
                dn_dt[i-1] -= flux / self.dL
        
        # Terme de nucléation (source à L=0)
        if G > 0:
            # Les nouveaux cristaux apparaissent dans la première classe
            dn_dt[0] += B / (G * self.dL)
        
        return np.concatenate([[dC_dt], dn_dt])
    
    def simuler(self, T_profile, dt=10.0, t_final=None):
        """
        Simule la cristallisation batch.
        
        Parameters
        ----------
        T_profile : function
            Fonction T(t) retournant la température [K]
        dt : float
            Pas de temps [s]
        t_final : float
            Temps final [s]
            
        Returns
        -------
        dict
            Résultats de la simulation
        """
        if t_final is None:
            t_final = t_total
        
        # Conditions initiales
        y0 = np.concatenate([[self.C], self.n])
        
        # Fonction pour l'intégrateur
        def dydt(t, y):
            return self._calculer_derivees(t, y, T_profile)
        
        # Intégration
        t_eval = np.arange(0, t_final + dt, dt)
        sol = solve_ivp(dydt, [0, t_final], y0, 
                       t_eval=t_eval, method='BDF', rtol=1e-6, atol=1e-9)
        
        # Extraction des résultats
        times = sol.t
        C_history = sol.y[0, :]
        n_history = sol.y[1:, :].T
        
        # Calcul des grandeurs supplémentaires
        S_history = []
        B_history = []
        G_history = []
        moments_history = []
        f_massique_history = []
        
        for i, t in enumerate(times):
            T = T_profile(t)
            C = C_history[i]
            n = n_history[i, :]
            
            m_solide = self._calculer_fraction_massique()
            S = sursaturation(C, T)
            B = taux_nucleation(S, m_solide)
            G = taux_croissance(S, T)
            moments = self._calculer_moments(n)
            
            S_history.append(S)
            B_history.append(B)
            G_history.append(G)
            moments_history.append(moments)
            
            # Fraction massique
            if C_history[0] > 0:
                f_mass = 1.0 - C / C_history[0]
            else:
                f_mass = 0.0
            f_massique_history.append(f_mass)
        
        return {
            'time': times,
            'C': C_history,
            'T': [T_profile(t) for t in times],
            'S': np.array(S_history),
            'B': np.array(B_history),
            'G': np.array(G_history),
            'n': n_history,
            'L': self.L,
            'moments': moments_history,
            'f_massique': np.array(f_massique_history)
        }


# ============================================
# 6. ANALYSE DES DISTRIBUTIONS DE TAILLE
# ============================================

def analyser_distribution(L, n_distribution):
    """
    Analyse une distribution de taille.
    
    Parameters
    ----------
    L : array
        Classes de taille [m]
    n_distribution : array
        Distribution [#/m³/m]
        
    Returns
    -------
    dict
        Indicateurs de qualité
    """
    if len(n_distribution) == 0 or np.sum(n_distribution) == 0:
        return {
            'L50': 0.0,
            'L_mean': 0.0,
            'L_std': 0.0,
            'CV': 0.0,
            'uniformite': 0.0
        }
    
    dL = L[1] - L[0]
    
    # Distribution cumulée
    dist_cum = np.cumsum(n_distribution) * dL
    dist_cum_norm = dist_cum / dist_cum[-1]
    
    # Taille médiane L50
    idx_50 = np.searchsorted(dist_cum_norm, 0.5)
    L50 = L[min(idx_50, len(L)-1)]
    
    # Moments
    m0 = np.sum(n_distribution) * dL
    m1 = np.sum(n_distribution * L) * dL
    m2 = np.sum(n_distribution * L**2) * dL
    
    # Statistiques
    L_mean = m1 / m0 if m0 > 0 else 0.0
    variance = m2 / m0 - L_mean**2 if m0 > 0 else 0.0
    L_std = np.sqrt(max(variance, 0))
    
    # Coefficient de variation
    CV = (L_std / L_mean * 100) if L_mean > 0 else 0.0
    
    # Indice d'uniformité
    L10_idx = np.searchsorted(dist_cum_norm, 0.1)
    L90_idx = np.searchsorted(dist_cum_norm, 0.9)
    L10 = L[min(L10_idx, len(L)-1)]
    L90 = L[min(L90_idx, len(L)-1)]
    
    uniformite = (L90 - L10) / L50 if L50 > 0 else 0.0
    
    return {
        'L50': L50 * 1e6,      # en µm
        'L_mean': L_mean * 1e6, # en µm
        'L_std': L_std * 1e6,   # en µm
        'CV': CV,               # en %
        'uniformite': uniformite,
        'L10': L10 * 1e6,
        'L90': L90 * 1e6
    }


# ============================================
# 7. COMPARAISON DES PROFILS DE REFROIDISSEMENT
# ============================================

def comparer_profils_refroidissement(C0=65.0, t_total=4*3600.0):
    """
    Compare les trois profils de refroidissement.
    
    Returns
    -------
    dict
        Résultats de comparaison
    """
    print("Simulation du profil linéaire...")
    
    # Création du cristalliseur
    cristalliseur = CristalliseurBatch(n_classes=80, L_min=1e-7, L_max=1000e-6)
    cristalliseur.set_conditions_initiales(C0, T0)
    
    # 1. Profil linéaire
    def T_lin(t):
        return profil_refroidissement_lineaire(t, T0, Tf, t_total)
    
    results_lin = cristalliseur.simuler(T_lin, dt=60.0, t_final=t_total)
    
    print("Simulation du profil exponentiel...")
    cristalliseur.set_conditions_initiales(C0, T0)
    
    # 2. Profil exponentiel
    def T_exp(t):
        return profil_refroidissement_exponentiel(t, T0, Tf, t_total)
    
    results_exp = cristalliseur.simuler(T_exp, dt=60.0, t_final=t_total)
    
    print("Simulation du profil optimal...")
    cristalliseur.set_conditions_initiales(C0, T0)
    
    # 3. Profil optimal (approximation)
    def T_opt(t):
        return profil_refroidissement_optimal(t, T0, S_cible=0.05, C0=C0)
    
    results_opt = cristalliseur.simuler(T_opt, dt=60.0, t_final=t_total)
    
    # Analyse des distributions finales
    n_final_lin = results_lin['n'][-1, :]
    n_final_exp = results_exp['n'][-1, :]
    n_final_opt = results_opt['n'][-1, :]
    
    analyse_lin = analyser_distribution(results_lin['L'], n_final_lin)
    analyse_exp = analyser_distribution(results_exp['L'], n_final_exp)
    analyse_opt = analyser_distribution(results_opt['L'], n_final_opt)
    
    return {
        'lin': {
            'results': results_lin,
            'analyse': analyse_lin,
            'rendement': results_lin['f_massique'][-1] * 100
        },
        'exp': {
            'results': results_exp,
            'analyse': analyse_exp,
            'rendement': results_exp['f_massique'][-1] * 100
        },
        'opt': {
            'results': results_opt,
            'analyse': analyse_opt,
            'rendement': results_opt['f_massique'][-1] * 100
        }
    }


# ============================================
# 8. DIMENSIONNEMENT DU CRISTALLISEUR
# ============================================

def dimensionner_cristalliseur_complet(
    production_batch_kg=5000.0,
    temps_batch_h=4.0,
    C0=65.0,
    f_crist_finale=0.85,
    rapport_H_D=1.5,
    fraction_remplissage=0.7,
    puissance_agitation_W_m3=1000.0,
    U_echange=500.0
):
    """
    Dimensionne un cristalliseur batch.
    
    Parameters
    ----------
    production_batch_kg : float
        Production de sucre par batch [kg]
    temps_batch_h : float
        Temps de batch [h]
    C0 : float
        Concentration initiale [g/100g]
    f_crist_finale : float
        Fraction cristallisée finale [-]
    rapport_H_D : float
        Rapport Hauteur/Diamètre
    fraction_remplissage : float
        Fraction de remplissage [-]
    puissance_agitation_W_m3 : float
        Puissance d'agitation spécifique [W/m³]
    U_echange : float
        Coefficient d'échange global [W/m²·K]
        
    Returns
    -------
    dict
        Caractéristiques du cristalliseur
    """
    # 1. Calcul des masses
    masse_sucre_total_kg = production_batch_kg / f_crist_finale
    masse_solution_kg = masse_sucre_total_kg / (C0 / 100.0)
    
    # 2. Volume de solution
    volume_solution_m3 = masse_solution_kg / rho_solution
    
    # 3. Volume du cristalliseur
    volume_cristalliseur_m3 = volume_solution_m3 / fraction_remplissage
    
    # 4. Dimensions géométriques
    diametre_m = (4 * volume_cristalliseur_m3 / (np.pi * rapport_H_D)) ** (1/3)
    hauteur_m = rapport_H_D * diametre_m
    
    # Vérification du volume
    volume_calc_m3 = np.pi * (diametre_m/2)**2 * hauteur_m
    
    # 5. Agitation
    puissance_agitation_W = puissance_agitation_W_m3 * volume_solution_m3
    puissance_agitation_kW = puissance_agitation_W / 1000.0
    
    # 6. Échange thermique
    # Chaleur à évacuer : refroidissement de T0 à Tf
    delta_T_mean = (T0 + Tf)/2 - 273.15  # Différence moyenne en °C
    Cp_solution = 3.5e3  # J/kg·K (capacité calorifique)
    
    charge_thermique_J = masse_solution_kg * Cp_solution * (T0 - Tf)
    puissance_thermique_W = charge_thermique_J / (temps_batch_h * 3600)
    
    # Surface d'échange nécessaire
    delta_T_log = (T0 - Tf) / np.log((T0 - 273.15) / (Tf - 273.15))  # Approximation
    surface_echange_m2 = puissance_thermique_W / (U_echange * delta_T_log)
    
    # 7. Taux de cisaillement (pour agitation)
    viscosite = 0.1  # Pa·s (approximation solution sucrée)
    vitesse_tip = 3.0  # m/s (vitesse en bout de pale)
    diametre_agitateur = diametre_m * 0.4  # 40% du diamètre du tank
    taux_cisaillement = vitesse_tip / (diametre_agitateur/2)  # 1/s
    
    return {
        'production_batch_kg': production_batch_kg,
        'temps_batch_h': temps_batch_h,
        'masse_solution_kg': masse_solution_kg,
        'masse_sucre_total_kg': masse_sucre_total_kg,
        'volume_solution_m3': volume_solution_m3,
        'volume_cristalliseur_m3': volume_cristalliseur_m3,
        'volume_calc_m3': volume_calc_m3,
        'diametre_m': diametre_m,
        'hauteur_m': hauteur_m,
        'rapport_H_D': rapport_H_D,
        'puissance_agitation_kW': puissance_agitation_kW,
        'charge_thermique_MJ': charge_thermique_J / 1e6,
        'puissance_thermique_kW': puissance_thermique_W / 1000,
        'surface_echange_m2': surface_echange_m2,
        'taux_cisaillement_1_s': taux_cisaillement,
        'nombre_reynolds': rho_solution * vitesse_tip * diametre_agitateur / viscosite
    }


# ============================================
# 9. VISUALISATION DES RÉSULTATS
# ============================================

def visualiser_resultats_complets(comparaison_results):
    """
    Visualise les résultats de la comparaison.
    """
    fig = plt.figure(figsize=(16, 12))
    gs = GridSpec(3, 3, figure=fig, hspace=0.3, wspace=0.3)
    
    # 1. Températures au cours du temps
    ax1 = fig.add_subplot(gs[0, 0])
    results_lin = comparaison_results['lin']['results']
    results_exp = comparaison_results['exp']['results']
    results_opt = comparaison_results['opt']['results']
    
    time_h = results_lin['time'] / 3600.0
    
    ax1.plot(time_h, np.array(results_lin['T']) - 273.15, 'b-', label='Linéaire', linewidth=2)
    ax1.plot(time_h, np.array(results_exp['T']) - 273.15, 'r--', label='Exponentiel', linewidth=2)
    ax1.plot(time_h, np.array(results_opt['T']) - 273.15, 'g-.', label='Optimal', linewidth=2)
    ax1.set_xlabel('Temps [h]')
    ax1.set_ylabel('Température [°C]')
    ax1.set_title('Profils de refroidissement')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. Sursaturation
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.plot(time_h, results_lin['S'], 'b-', label='Linéaire', linewidth=2)
    ax2.plot(time_h, results_exp['S'], 'r--', label='Exponentiel', linewidth=2)
    ax2.plot(time_h, results_opt['S'], 'g-.', label='Optimal', linewidth=2)
    ax2.axhline(y=0.05, color='k', linestyle=':', alpha=0.5, label='S cible=0.05')
    ax2.set_xlabel('Temps [h]')
    ax2.set_ylabel('Sursaturation S [-]')
    ax2.set_title('Évolution de la sursaturation')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 3. Rendement
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.plot(time_h, results_lin['f_massique']*100, 'b-', label='Linéaire', linewidth=2)
    ax3.plot(time_h, results_exp['f_massique']*100, 'r--', label='Exponentiel', linewidth=2)
    ax3.plot(time_h, results_opt['f_massique']*100, 'g-.', label='Optimal', linewidth=2)
    ax3.set_xlabel('Temps [h]')
    ax3.set_ylabel('Rendement [%]')
    ax3.set_title('Évolution du rendement')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # 4. Distributions de taille finales
    ax4 = fig.add_subplot(gs[1, :])
    
    L_um = results_lin['L'] * 1e6
    
    # Normalisation des distributions
    n_lin = comparaison_results['lin']['results']['n'][-1, :]
    n_exp = comparaison_results['exp']['results']['n'][-1, :]
    n_opt = comparaison_results['opt']['results']['n'][-1, :]
    
    if np.sum(n_lin) > 0:
        n_lin_norm = n_lin / np.max(n_lin)
        ax4.plot(L_um, n_lin_norm, 'b-', label=f"Linéaire (L50={comparaison_results['lin']['analyse']['L50']:.1f} µm)", linewidth=2)
    
    if np.sum(n_exp) > 0:
        n_exp_norm = n_exp / np.max(n_exp)
        ax4.plot(L_um, n_exp_norm, 'r--', label=f"Exponentiel (L50={comparaison_results['exp']['analyse']['L50']:.1f} µm)", linewidth=2)
    
    if np.sum(n_opt) > 0:
        n_opt_norm = n_opt / np.max(n_opt)
        ax4.plot(L_um, n_opt_norm, 'g-.', label=f"Optimal (L50={comparaison_results['opt']['analyse']['L50']:.1f} µm)", linewidth=2)
    
    ax4.set_xlabel('Taille des cristaux [µm]')
    ax4.set_ylabel('Distribution normalisée [-]')
    ax4.set_title('Distributions de taille finales')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    # 5. Tableau comparatif
    ax5 = fig.add_subplot(gs[2, :])
    ax5.axis('tight')
    ax5.axis('off')
    
    # Données du tableau
    data = [
        ['Profil', 'L50 [µm]', 'CV [%]', 'Rendement [%]', 'L_mean [µm]', 'Uniformité'],
        ['Linéaire', 
         f"{comparaison_results['lin']['analyse']['L50']:.1f}",
         f"{comparaison_results['lin']['analyse']['CV']:.1f}",
         f"{comparaison_results['lin']['rendement']:.1f}",
         f"{comparaison_results['lin']['analyse']['L_mean']:.1f}",
         f"{comparaison_results['lin']['analyse']['uniformite']:.2f}"],
        ['Exponentiel',
         f"{comparaison_results['exp']['analyse']['L50']:.1f}",
         f"{comparaison_results['exp']['analyse']['CV']:.1f}",
         f"{comparaison_results['exp']['rendement']:.1f}",
         f"{comparaison_results['exp']['analyse']['L_mean']:.1f}",
         f"{comparaison_results['exp']['analyse']['uniformite']:.2f}"],
        ['Optimal',
         f"{comparaison_results['opt']['analyse']['L50']:.1f}",
         f"{comparaison_results['opt']['analyse']['CV']:.1f}",
         f"{comparaison_results['opt']['rendement']:.1f}",
         f"{comparaison_results['opt']['analyse']['L_mean']:.1f}",
         f"{comparaison_results['opt']['analyse']['uniformite']:.2f}"]
    ]
    
    table = ax5.table(cellText=data, loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)
    
    plt.suptitle('Comparaison des profils de refroidissement - Cristallisation batch', fontsize=14, y=0.98)
    plt.tight_layout()
    plt.show()


def visualiser_dimensionnement(dimensionnement):
    """
    Visualise les résultats du dimensionnement.
    """
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # 1. Diagramme du cristalliseur
    ax1 = axes[0, 0]
    D = dimensionnement['diametre_m']
    H = dimensionnement['hauteur_m']
    
    # Dessin simplifié
    circle = plt.Circle((0.5, 0.3), D/10, color='lightblue', alpha=0.5)
    ax1.add_patch(circle)
    rect = plt.Rectangle((0.5 - D/10, 0.3), D/5, H/10, color='lightblue', alpha=0.5)
    ax1.add_patch(rect)
    
    ax1.text(0.5, 0.25, f"D = {D:.2f} m", ha='center', va='center')
    ax1.text(0.5, 0.35 + H/20, f"H = {H:.2f} m", ha='center', va='center')
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    ax1.set_aspect('equal')
    ax1.axis('off')
    ax1.set_title('Schéma du cristalliseur')
    
    # 2. Données clés
    ax2 = axes[0, 1]
    ax2.axis('tight')
    ax2.axis('off')
    
    data_dim = [
        ['Paramètre', 'Valeur', 'Unité'],
        ['Volume total', f"{dimensionnement['volume_cristalliseur_m3']:.1f}", 'm³'],
        ['Volume solution', f"{dimensionnement['volume_solution_m3']:.1f}", 'm³'],
        ['Diamètre', f"{dimensionnement['diametre_m']:.2f}", 'm'],
        ['Hauteur', f"{dimensionnement['hauteur_m']:.2f}", 'm'],
        ['Rapport H/D', f"{dimensionnement['rapport_H_D']:.1f}", '-'],
    ]
    
    table_dim = ax2.table(cellText=data_dim, loc='center', cellLoc='center')
    table_dim.auto_set_font_size(False)
    table_dim.set_fontsize(9)
    table_dim.scale(1, 1.5)
    
    # 3. Données thermiques
    ax3 = axes[0, 2]
    ax3.axis('tight')
    ax3.axis('off')
    
    data_therm = [
        ['Paramètre', 'Valeur', 'Unité'],
        ['Puissance agitation', f"{dimensionnement['puissance_agitation_kW']:.1f}", 'kW'],
        ['Charge thermique', f"{dimensionnement['charge_thermique_MJ']:.1f}", 'MJ'],
        ['Puissance thermique', f"{dimensionnement['puissance_thermique_kW']:.1f}", 'kW'],
        ['Surface échange', f"{dimensionnement['surface_echange_m2']:.1f}", 'm²'],
        ['Taux cisaillement', f"{dimensionnement['taux_cisaillement_1_s']:.0f}", '1/s'],
    ]
    
    table_therm = ax3.table(cellText=data_therm, loc='center', cellLoc='center')
    table_therm.auto_set_font_size(False)
    table_therm.set_fontsize(9)
    table_therm.scale(1, 1.5)
    
    # 4. Production
    ax4 = axes[1, 0]
    categories = ['Solution', 'Sucre total', 'Sucre cristallisé']
    valeurs = [
        dimensionnement['masse_solution_kg'] / 1000,  # en tonnes
        dimensionnement['masse_sucre_total_kg'] / 1000,
        dimensionnement['production_batch_kg'] / 1000
    ]
    
    bars = ax4.bar(categories, valeurs, color=['lightblue', 'lightgreen', 'gold'])
    ax4.set_ylabel('Masse [tonnes]')
    ax4.set_title('Masses par batch')
    ax4.grid(True, alpha=0.3, axis='y')
    
    for bar, val in zip(bars, valeurs):
        height = bar.get_height()
        ax4.text(bar.get_x() + bar.get_width()/2., height,
                f'{val:.1f} t', ha='center', va='bottom')
    
    # 5. Volumes
    ax5 = axes[1, 1]
    labels_vol = ['Volume cristalliseur', 'Volume solution']
    sizes_vol = [dimensionnement['volume_cristalliseur_m3'], 
                dimensionnement['volume_solution_m3']]
    colors_vol = ['lightcoral', 'lightskyblue']
    
    ax5.pie(sizes_vol, labels=labels_vol, colors=colors_vol, autopct='%1.1f%%',
           startangle=90)
    ax5.set_title('Répartition des volumes')
    
    # 6. Puissances
    ax6 = axes[1, 2]
    categories_puiss = ['Agitation', 'Refroidissement']
    valeurs_puiss = [dimensionnement['puissance_agitation_kW'],
                    dimensionnement['puissance_thermique_kW']]
    
    bars_puiss = ax6.bar(categories_puiss, valeurs_puiss, color=['orange', 'lightseagreen'])
    ax6.set_ylabel('Puissance [kW]')
    ax6.set_title('Besoins en puissance')
    ax6.grid(True, alpha=0.3, axis='y')
    
    for bar, val in zip(bars_puiss, valeurs_puiss):
        height = bar.get_height()
        ax6.text(bar.get_x() + bar.get_width()/2., height,
                f'{val:.1f} kW', ha='center', va='bottom')
    
    plt.suptitle('Dimensionnement du cristalliseur batch', fontsize=14, y=0.98)
    plt.tight_layout()
    plt.show()


# ============================================
# 10. FONCTION PRINCIPALE DE DÉMONSTRATION
# ============================================

def demo_cristallisation_complete():
    """
    Démonstration complète du module de cristallisation.
    """
    print("=" * 80)
    print("DÉMONSTRATION COMPLÈTE - MODULE DE CRISTALLISATION")
    print("=" * 80)
    
    # Partie 1 : Thermodynamique
    print("\n1. VÉRIFICATION THERMODYNAMIQUE")
    print("-" * 40)
    
    T_test = np.array([20, 40, 60, 80]) + 273.15
    for T in T_test:
        C_star = solubilite_saccharose(T)
        print(f"  Température {T-273.15:3.0f}°C : Solubilité = {C_star:6.1f} g/100g")
    
    # Partie 2 : Comparaison des profils
    print("\n2. COMPARAISON DES PROFILS DE REFROIDISSEMENT")
    print("-" * 40)
    print("  Simulation en cours... (cela peut prendre quelques minutes)")
    
    comparaison = comparer_profils_refroidissement(C0=65.0, t_total=4*3600.0)
    
    print("\n  Résultats finaux :")
    print("  " + "-" * 60)
    print(f"  | Profil        | L50 [µm] | CV [%] | Rendement [%] |")
    print("  " + "-" * 60)
    print(f"  | Linéaire      | {comparaison['lin']['analyse']['L50']:8.1f} | {comparaison['lin']['analyse']['CV']:6.1f} | {comparaison['lin']['rendement']:13.1f} |")
    print(f"  | Exponentiel   | {comparaison['exp']['analyse']['L50']:8.1f} | {comparaison['exp']['analyse']['CV']:6.1f} | {comparaison['exp']['rendement']:13.1f} |")
    print(f"  | Optimal       | {comparaison['opt']['analyse']['L50']:8.1f} | {comparaison['opt']['analyse']['CV']:6.1f} | {comparaison['opt']['rendement']:13.1f} |")
    print("  " + "-" * 60)
    
    # Partie 3 : Dimensionnement
    print("\n3. DIMENSIONNEMENT DU CRISTALLISEUR")
    print("-" * 40)
    
    dimensionnement = dimensionner_cristalliseur_complet(
        production_batch_kg=5000.0,
        temps_batch_h=4.0,
        C0=65.0,
        f_crist_finale=0.85
    )
    
    print(f"  Production par batch : {dimensionnement['production_batch_kg']:.0f} kg")
    print(f"  Volume cristalliseur : {dimensionnement['volume_cristalliseur_m3']:.1f} m³")
    print(f"  Dimensions : D = {dimensionnement['diametre_m']:.2f} m, H = {dimensionnement['hauteur_m']:.2f} m")
    print(f"  Puissance agitation  : {dimensionnement['puissance_agitation_kW']:.1f} kW")
    print(f"  Surface échange      : {dimensionnement['surface_echange_m2']:.1f} m²")
    
    # Partie 4 : Visualisation
    print("\n4. GÉNÉRATION DES GRAPHIQUES")
    print("-" * 40)
    print("  Génération des visualisations...")
    
    visualiser_resultats_complets(comparaison)
    visualiser_dimensionnement(dimensionnement)
    
    print("\n" + "=" * 80)
    print("DÉMONSTRATION TERMINÉE AVEC SUCCÈS")
    print("=" * 80)
    
    return comparaison, dimensionnement


# ============================================
# 11. EXÉCUTION
# ============================================

if __name__ == "__main__":
    # Exécute la démonstration complète
    comparaison, dimensionnement = demo_cristallisation_complete()
    
    # Option : sauvegarde des résultats
    import pickle
    with open('resultats_cristallisation.pkl', 'wb') as f:
        pickle.dump({'comparaison': comparaison, 'dimensionnement': dimensionnement}, f)
    print("\nRésultats sauvegardés dans 'resultats_cristallisation.pkl'")