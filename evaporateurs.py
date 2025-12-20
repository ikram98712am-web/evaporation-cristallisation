"""
evaporateurs.py
----------------
Modèles d'évaporation améliorés avec :
- 1 effet (classe EvaporateurEffet)
- système multi-effets (fonction systeme_multi_effets)
- Optimisation et analyse de sensibilité
"""

import numpy as np
from scipy.optimize import fsolve, minimize
import matplotlib.pyplot as plt
from thermodynamique import (
    temperature_ebullition_solution,
    enthalpie_solution_sucree,
    chaleur_latente_vaporisation,
    temperature_ebullition_eau,
    epe_saccharose_duhring,  # À implémenter selon corrélation de Dühring
    enthalpie_vapeur_surchauffee,
)


class EvaporateurEffet:
    """
    Représente un évaporateur simple (un effet).
    Bilans matière + énergie complets avec EPE.
    """

    def __init__(self, numero, P_effet, P_vapeur_chaude, U, 
                 A=None, pertes_frac=0.03, Rf=0.0002, 
                 surchauffe_vapeur=10.0):
        """
        numero           : index de l'effet (1, 2, 3, ...)
        P_effet          : pression interne [Pa]
        P_vapeur_chaude  : pression vapeur de chauffe [Pa]
        U                : coef. global d'échange [kW/m2/K]
        A                : surface (si déjà dimensionnée), sinon None
        pertes_frac      : fraction de pertes de chaleur (ex: 0.03 = 3%)
        Rf               : résistance d'encrassement [m²·K/W]
        surchauffe_vapeur: surchauffe de la vapeur de chauffe [K]
        """
        self.numero = numero
        self.P_effet = P_effet
        self.P_vapeur_chaude = P_vapeur_chaude
        self.U_nom = U
        self.U_eff = 1/(1/U + Rf) if U > 0 else 0  # Correction encrassement
        self.A = A
        self.pertes_frac = pertes_frac
        self.Rf = Rf
        self.surchauffe = surchauffe_vapeur
        
    def calculer_bilans(self, F_in, X_in, T_in, X_out_cible, 
                       vapeur_secondaire=None):
        """
        Calcule les bilans pour cet effet avec EPE.
        
        Entrées :
        - F_in  : débit d'alimentation [kg/h]
        - X_in  : fraction massique saccharose en entrée
        - T_in  : température d'entrée [K]
        - X_out_cible : fraction massique saccharose voulue en sortie
        - vapeur_secondaire : tuple (débit, pression) si vapeur de l'effet précédent
        
        Sortie : dict avec résultats complets
        """
        # ------------ Bilans matière ------------
        if X_out_cible <= 0 or X_out_cible <= X_in:
            raise ValueError("X_out_cible doit être > X_in et > 0")
        
        L_out = F_in * X_in / X_out_cible        # kg/h liquide en sortie
        V_out = F_in - L_out                     # kg/h d'eau évaporée
        
        # ------------ Températures avec EPE ------------
        # Température d'ébullition de l'eau pure à P_effet
        T_eb_eau = temperature_ebullition_eau(self.P_effet)
        
        # Élévation du point d'ébullition (corrélation de Dühring)
        EPE = epe_saccharose_duhring(X_out_cible, T_eb_eau)
        
        # Température d'ébullition de la solution
        T_eb_solution = T_eb_eau + EPE
        
        # Température vapeur de chauffe (surchauffée)
        T_sat_steam = temperature_ebullition_eau(self.P_vapeur_chaude)
        T_vapeur_chaude = T_sat_steam + self.surchauffe
        
        # ------------ Bilans énergétiques complets ------------
        # Enthalpies solution (kJ/kg)
        h_in = enthalpie_solution_sucree(T_in, X_in)
        h_out = enthalpie_solution_sucree(T_eb_solution, X_out_cible)
        
        # Chauffage solution jusqu'à T_eb
        Q_sensible = F_in * (h_out - h_in)   # kJ/h
        
        # Enthalpie de la vapeur produite (vapeur saturée à P_effet)
        h_vapeur_saturee = enthalpie_vapeur_surchauffee(self.P_effet, T_eb_eau)
        
        # Chaleur latente pour évaporer V_out
        Q_latent = V_out * (h_vapeur_saturee - h_out)  # kJ/h
        
        # Pertes de chaleur
        Q_total_sans_pertes = Q_sensible + Q_latent
        Q_total_avec_pertes = Q_total_sans_pertes / (1.0 - self.pertes_frac)
        Q_total_kW = Q_total_avec_pertes / 3600.0  # kW
        
        # ------------ Débit vapeur de chauffe requis ------------
        if vapeur_secondaire is None:
            # Vapeur vive (premier effet)
            h_vapeur_chaude = enthalpie_vapeur_surchauffee(
                self.P_vapeur_chaude, T_vapeur_chaude
            )
            h_condensat = enthalpie_solution_sucree(T_sat_steam, 0)
            delta_h_steam = h_vapeur_chaude - h_condensat
        else:
            # Vapeur secondaire (effets suivants)
            V_sec, P_sec = vapeur_secondaire
            T_sat_sec = temperature_ebullition_eau(P_sec)
            h_vapeur_sec = enthalpie_vapeur_surchauffee(P_sec, T_sat_sec)
            h_condensat_sec = enthalpie_solution_sucree(T_sat_sec, 0)
            delta_h_steam = h_vapeur_sec - h_condensat_sec
        
        m_vapeur_chaude = Q_total_avec_pertes / delta_h_steam  # kg/h
        
        # ------------ Surface d'échange ------------
        # ΔT utile (vapeur chaude -> solution bouillante)
        deltaT = max(T_vapeur_chaude - T_eb_solution, 1.0)  # K, min 1 K
        A_calc = Q_total_kW / (self.U_eff * deltaT)         # m²
        
        return {
            "numero": self.numero,
            "F_in": F_in,
            "X_in": X_in,
            "T_in": T_in,
            "L_out": L_out,
            "X_out": X_out_cible,
            "V_out": V_out,
            "T_eb_eau": T_eb_eau,
            "EPE": EPE,
            "T_eb_solution": T_eb_solution,
            "T_vapeur_chaude": T_vapeur_chaude,
            "Q_sensible_kJh": Q_sensible,
            "Q_latent_kJh": Q_latent,
            "Q_total_kW": Q_total_kW,
            "m_vapeur_chaude": m_vapeur_chaude,
            "deltaT": deltaT,
            "A_calc": A_calc,
            "U_nom": self.U_nom,
            "U_eff": self.U_eff,
            "P_effet": self.P_effet,
            "P_vapeur_chaude": self.P_vapeur_chaude,
        }


def systeme_multi_effets(
    F_alim,
    X_alim,
    T_alim,
    pressions_effets,
    P_vapeur_vive,
    U_list,
    X_sortie_cible,
    pertes_frac=0.03,
    Rf=0.0002,
    surchauffe=10.0,
    methode='sequentielle',
    tol=1e-6,
    max_iter=100
):
    """
    Système multi-effets en série, co-courant.
    
    Paramètres :
    ------------
    F_alim, X_alim, T_alim : alim globale [kg/h, -, K]
    pressions_effets : liste [Pa] pour chaque effet
    P_vapeur_vive : pression vapeur vive pour le 1er effet [Pa]
    U_list : liste de U [kW/m2/K] pour chaque effet
    X_sortie_cible : X visée en sortie du DERNIER effet
    methode : 'sequentielle' ou 'globale' (résolution simultanée)
    
    Retour : liste de dict (un par effet) et résultats globaux
    """
    n = len(pressions_effets)
    if len(U_list) != n:
        raise ValueError("U_list doit avoir la même longueur que pressions_effets")
    
    if methode == 'sequentielle':
        return _solve_sequentiel(
            F_alim, X_alim, T_alim, pressions_effets, P_vapeur_vive,
            U_list, X_sortie_cible, pertes_frac, Rf, surchauffe
        )
    else:
        return _solve_global(
            F_alim, X_alim, T_alim, pressions_effets, P_vapeur_vive,
            U_list, X_sortie_cible, pertes_frac, Rf, surchauffe,
            tol, max_iter
        )


def _solve_sequentiel(F_alim, X_alim, T_alim, pressions_effets, 
                     P_vapeur_vive, U_list, X_sortie_cible,
                     pertes_frac, Rf, surchauffe):
    """Résolution séquentielle (plus simple mais approximative)."""
    n = len(pressions_effets)
    
    # Distribution des concentrations (linéaire en fraction de solides)
    solides_in = F_alim * X_alim
    solides_out = solides_in / X_sortie_cible * X_sortie_cible
    
    # Calcul des débits intermédiaires (conservation solides)
    L_vals = [F_alim]
    X_vals = [X_alim]
    
    for i in range(n):
        # Estimation de la concentration (progression géométrique)
        if i < n-1:
            X_next = X_vals[i] * (X_sortie_cible / X_alim) ** (1/n)
        else:
            X_next = X_sortie_cible
        X_vals.append(X_next)
        
        L_next = solides_in / X_next
        L_vals.append(L_next)
    
    effets_res = []
    F = F_alim
    X = X_alim
    T = T_alim
    P_steam = P_vapeur_vive
    vapeur_precedente = None
    
    for i in range(n):
        effet = EvaporateurEffet(
            numero=i + 1,
            P_effet=pressions_effets[i],
            P_vapeur_chaude=P_steam,
            U=U_list[i],
            pertes_frac=pertes_frac,
            Rf=Rf,
            surchauffe_vapeur=surchauffe if i == 0 else 0.0  # Surchauffe seulement 1er effet
        )
        
        res = effet.calculer_bilans(
            F_in=F,
            X_in=X,
            T_in=T,
            X_out_cible=X_vals[i+1],
            vapeur_secondaire=vapeur_precedente
        )
        effets_res.append(res)
        
        # Mise à jour pour l'effet suivant
        F = res["L_out"]
        X = res["X_out"]
        T = res["T_eb_solution"]
        
        # La vapeur produite devient vapeur de chauffe pour l'effet suivant
        vapeur_precedente = (res["V_out"], pressions_effets[i])
        P_steam = pressions_effets[i]
    
    return effets_res


def _solve_global(F_alim, X_alim, T_alim, pressions_effets,
                 P_vapeur_vive, U_list, X_sortie_cible,
                 pertes_frac, Rf, surchauffe, tol, max_iter):
    """Résolution globale simultanée (plus précise)."""
    n = len(pressions_effets)
    
    # Variables : V1, V2, ..., Vn, T1, T2, ..., Tn
    # Initial guess
    V_total_est = F_alim * (1 - X_alim / X_sortie_cible)
    V_guess = [V_total_est / n] * n
    
    T_sat_steam = temperature_ebullition_eau(P_vapeur_vive)
    T_cond = temperature_ebullition_eau(pressions_effets[-1])
    T_guess = np.linspace(T_sat_steam - 10, T_cond + 10, n)
    
    x0 = np.concatenate([V_guess, T_guess])
    
    def equations(x):
        V = x[:n]
        T = x[n:]
        
        # Bilans matière
        L = [F_alim]
        X = [X_alim]
        for i in range(n):
            L_next = L[i] - V[i]
            X_next = L[i] * X[i] / L_next if L_next > 0 else X[i]
            L.append(L_next)
            X.append(X_next)
        
        # Vérification concentration finale
        eqs = []
        
        # Équation concentration finale
        eqs.append(X[-1] - X_sortie_cible)
        
        # Bilans énergie
        for i in range(n):
            if i == 0:
                # Premier effet : vapeur vive
                T_sat_steam = temperature_ebullition_eau(P_vapeur_vive)
                T_vap = T_sat_steam + surchauffe
                h_vap = enthalpie_vapeur_surchauffee(P_vapeur_vive, T_vap)
                h_cond = enthalpie_solution_sucree(T_sat_steam, 0)
                delta_h = h_vap - h_cond
                Q_vap = V[i] * chaleur_latente_vaporisation(pressions_effets[i])
                
                # Bilan approximatif
                eqs.append(Q_vap - delta_h * 0.9)  # Facteur d'efficacité
            else:
                # Effets suivants : vapeur de l'effet précédent
                Q_vap = V[i] * chaleur_latente_vaporisation(pressions_effets[i])
                Q_dispo = V[i-1] * chaleur_latente_vaporisation(pressions_effets[i-1])
                eqs.append(Q_vap - Q_dispo * 0.95)  # Pertes entre effets
        
        return eqs
    
    # Résolution
    sol = fsolve(equations, x0, xtol=tol, maxfev=max_iter)
    
    # Reconstruction des résultats
    V = sol[:n]
    T = sol[n:]
    
    effets_res = []
    L = F_alim
    X = X_alim
    Temp = T_alim
    
    for i in range(n):
        L_next = L - V[i]
        X_next = L * X / L_next
        
        effet = EvaporateurEffet(
            numero=i + 1,
            P_effet=pressions_effets[i],
            P_vapeur_chaude=P_vapeur_vive if i == 0 else pressions_effets[i-1],
            U=U_list[i],
            pertes_frac=pertes_frac,
            Rf=Rf,
            surchauffe_vapeur=surchauffe if i == 0 else 0.0
        )
        
        # Calcul complet avec les valeurs trouvées
        res = effet.calculer_bilans(
            F_in=L,
            X_in=X,
            T_in=Temp,
            X_out_cible=X_next,
            vapeur_secondaire=None if i == 0 else (V[i-1], pressions_effets[i-1])
        )
        
        # Ajustement avec la température calculée
        res["T_eb_solution"] = T[i]
        effets_res.append(res)
        
        # Mise à jour pour l'effet suivant
        L = L_next
        X = X_next
        Temp = T[i]
    
    return effets_res


def calculer_economy_globale(effets):
    """
    Calcule l'économie de vapeur globale.
    
    Parameters
    ----------
    effets : list
        Liste des résultats par effet
        
    Returns
    -------
    float
        Économie de vapeur (eau évaporée / vapeur vive)
    """
    if not effets:
        return 0.0
    
    vapeur_vive = effets[0].get("m_vapeur_chaude", 0)
    if vapeur_vive <= 0:
        return 0.0
    
    eau_totale_evaporee = sum(e.get("V_out", 0) for e in effets)
    return eau_totale_evaporee / vapeur_vive


def optimiser_nombre_effets(
    F_alim,
    X_alim,
    T_alim,
    P_vapeur_vive,
    P_condenseur,
    X_sortie_cible,
    n_min=2,
    n_max=5,
    cout_investissement=1500,  # €/m²
    cout_energie=30,           # €/tonne vapeur
    heures_an=8000
):
    """
    Optimisation du nombre d'effets.
    
    Returns
    -------
    dict
        Résultats pour chaque configuration
    """
    resultats = {}
    
    for n in range(n_min, n_max + 1):
        # Distribution linéaire des pressions
        pressions = np.linspace(P_vapeur_vive, P_condenseur, n)
        
        # Coefficients U typiques (décroissants avec n)
        U_base = [2.5, 2.2, 1.8, 1.5, 1.3]  # kW/m²/K
        U_list = U_base[:n]
        
        # Simulation
        effets = systeme_multi_effets(
            F_alim=F_alim,
            X_alim=X_alim,
            T_alim=T_alim,
            pressions_effets=pressions,
            P_vapeur_vive=P_vapeur_vive,
            U_list=U_list,
            X_sortie_cible=X_sortie_cible,
            methode='sequentielle'
        )
        
        # Calculs économiques
        surface_totale = sum(e.get("A_calc", 0) for e in effets)
        vapeur_vive = effets[0].get("m_vapeur_chaude", 0)
        economy = calculer_economy_globale(effets)
        
        # Coûts
        investissement = surface_totale * cout_investissement
        energie_annuelle = (vapeur_vive * heures_an / 1000) * cout_energie
        cout_total_annuel = (investissement / 10) + energie_annuelle  # Amortissement 10 ans
        
        resultats[n] = {
            'n_effets': n,
            'surface_totale_m2': surface_totale,
            'vapeur_vive_kg_h': vapeur_vive,
            'economy': economy,
            'investissement_euros': investissement,
            'energie_annuelle_euros': energie_annuelle,
            'cout_total_annuel_euros': cout_total_annuel,
            'effets': effets
        }
    
    return resultats


def analyse_sensibilite(
    F_alim_base,
    X_alim_base,
    T_alim_base,
    P_vapeur_base,
    P_condenseur_base,
    X_sortie_base,
    n_effets=3,
    variations=None
):
    """
    Analyse de sensibilité paramétrique.
    
    variations : dict avec paramètres et plages de variation
    """
    if variations is None:
        variations = {
            'P_vapeur': np.linspace(2.5e5, 4.5e5, 5),      # 2.5 à 4.5 bar
            'X_sortie': np.linspace(0.60, 0.70, 5),       # 60% à 70%
            'F_alim': np.linspace(0.8*F_alim_base, 1.2*F_alim_base, 5),
            'T_alim': np.linspace(75+273, 95+273, 5),     # 75 à 95°C
        }
    
    resultats_sens = {}
    
    for param, valeurs in variations.items():
        resultats_param = []
        
        for val in valeurs:
            # Préparation des paramètres
            params = {
                'F_alim': F_alim_base,
                'X_alim': X_alim_base,
                'T_alim': T_alim_base,
                'P_vapeur_vive': P_vapeur_base,
                'P_condenseur': P_condenseur_base,
                'X_sortie_cible': X_sortie_base,
                'n_effets': n_effets
            }
            
            # Modification du paramètre étudié
            if param == 'P_vapeur':
                params['P_vapeur_vive'] = val
            elif param == 'X_sortie':
                params['X_sortie_cible'] = val
            elif param == 'F_alim':
                params['F_alim'] = val
            elif param == 'T_alim':
                params['T_alim'] = val
            
            # Distribution des pressions
            pressions = np.linspace(
                params['P_vapeur_vive'], 
                params['P_condenseur'], 
                n_effets
            )
            
            # Coefficients U
            U_base = [2.5, 2.2, 1.8, 1.5, 1.3]
            U_list = U_base[:n_effets]
            
            # Simulation
            try:
                effets = systeme_multi_effets(
                    F_alim=params['F_alim'],
                    X_alim=params['X_alim'],
                    T_alim=params['T_alim'],
                    pressions_effets=pressions,
                    P_vapeur_vive=params['P_vapeur_vive'],
                    U_list=U_list,
                    X_sortie_cible=params['X_sortie_cible']
                )
                
                # Extraction des indicateurs
                vapeur_vive = effets[0].get("m_vapeur_chaude", 0)
                surface_totale = sum(e.get("A_calc", 0) for e in effets)
                economy = calculer_economy_globale(effets)
                temperatures = [e.get("T_eb_solution", 0) for e in effets]
                
                resultats_param.append({
                    'valeur_param': val,
                    'vapeur_vive_kg_h': vapeur_vive,
                    'surface_totale_m2': surface_totale,
                    'economy': economy,
                    'temperatures_K': temperatures
                })
            except Exception as e:
                print(f"Erreur pour {param}={val}: {e}")
                resultats_param.append(None)
        
        resultats_sens[param] = resultats_param
    
    return resultats_sens


def visualiser_resultats(effets, title="Résultats évaporation"):
    """Visualisation graphique des résultats."""
    if not effets:
        return
    
    n = len(effets)
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(title, fontsize=16)
    
    # 1. Températures et pressions
    ax1 = axes[0, 0]
    numeros = [e['numero'] for e in effets]
    T_eb = [e['T_eb_solution'] - 273.15 for e in effets]  # Conversion en °C
    P_eff = [e['P_effet'] / 1e5 for e in effets]  # Conversion en bar
    
    ax1.plot(numeros, T_eb, 'bo-', label='T ébullition (°C)')
    ax1.set_xlabel('Numéro effet')
    ax1.set_ylabel('Température (°C)', color='b')
    ax1.tick_params(axis='y', labelcolor='b')
    ax1.grid(True, alpha=0.3)
    
    ax1b = ax1.twinx()
    ax1b.plot(numeros, P_eff, 'rs--', label='Pression (bar)')
    ax1b.set_ylabel('Pression (bar)', color='r')
    ax1b.tick_params(axis='y', labelcolor='r')
    
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax1b.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')
    
    # 2. Concentrations et débits
    ax2 = axes[0, 1]
    X_vals = [e['X_out'] * 100 for e in effets]  # Conversion en %
    L_vals = [e['L_out'] for e in effets]
    
    ax2.plot(numeros, X_vals, 'go-', label='Concentration (%)')
    ax2.set_xlabel('Numéro effet')
    ax2.set_ylabel('Concentration (%)', color='g')
    ax2.tick_params(axis='y', labelcolor='g')
    ax2.grid(True, alpha=0.3)
    
    ax2b = ax2.twinx()
    ax2b.plot(numeros, L_vals, 'ms--', label='Débit liquide (kg/h)')
    ax2b.set_ylabel('Débit liquide (kg/h)', color='m')
    ax2b.tick_params(axis='y', labelcolor='m')
    
    lines1, labels1 = ax2.get_legend_handles_labels()
    lines2, labels2 = ax2b.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, loc='upper right')
    
    # 3. Surfaces et puissances
    ax3 = axes[1, 0]
    A_vals = [e['A_calc'] for e in effets]
    Q_vals = [e['Q_total_kW'] for e in effets]
    
    bars = ax3.bar(numeros, A_vals, alpha=0.7, color='c', label='Surface (m²)')
    ax3.set_xlabel('Numéro effet')
    ax3.set_ylabel('Surface (m²)', color='c')
    ax3.tick_params(axis='y', labelcolor='c')
    
    # Ajout des valeurs sur les barres
    for bar in bars:
        height = bar.get_height()
        ax3.text(bar.get_x() + bar.get_width()/2., height,
                f'{height:.1f}', ha='center', va='bottom')
    
    ax3b = ax3.twinx()
    ax3b.plot(numeros, Q_vals, 'y^-', linewidth=2, label='Puissance (kW)')
    ax3b.set_ylabel('Puissance (kW)', color='y')
    ax3b.tick_params(axis='y', labelcolor='y')
    
    lines1, labels1 = ax3.get_legend_handles_labels()
    lines2, labels2 = ax3b.get_legend_handles_labels()
    ax3.legend(lines1 + lines2, labels1 + labels2, loc='upper left')
    ax3.grid(True, alpha=0.3)
    
    # 4. Tableau récapitulatif
    ax4 = axes[1, 1]
    ax4.axis('tight')
    ax4.axis('off')
    
    # Données du tableau
    tableau_data = []
    headers = ['Effet', 'L (kg/h)', 'X (%)', 'T (°C)', 'V (kg/h)', 'A (m²)']
    
    for e in effets:
        tableau_data.append([
            e['numero'],
            f"{e['L_out']:.0f}",
            f"{e['X_out']*100:.1f}",
            f"{e['T_eb_solution']-273.15:.1f}",
            f"{e['V_out']:.0f}",
            f"{e['A_calc']:.1f}"
        ])
    
    # Ajout des totaux
    tableau_data.append([
        'Total',
        f"{effets[-1]['L_out']:.0f}",
        f"{effets[-1]['X_out']*100:.1f}",
        '-',
        f"{sum(e['V_out'] for e in effets):.0f}",
        f"{sum(e['A_calc'] for e in effets):.1f}"
    ])
    
    # Création du tableau
    table = ax4.table(
        cellText=tableau_data,
        colLabels=headers,
        cellLoc='center',
        loc='center',
        colWidths=[0.15, 0.15, 0.15, 0.15, 0.15, 0.15]
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.5)
    
    plt.tight_layout()
    plt.show()


def example_utilisation():
    """Exemple d'utilisation des fonctions."""
    # Données de base
    F_alim = 20000.0  # kg/h
    X_alim = 0.15     # 15%
    T_alim = 85 + 273.15  # K
    P_vapeur_vive = 3.5e5  # Pa (3.5 bar)
    P_condenseur = 0.15e5  # Pa (0.15 bar)
    X_sortie_cible = 0.65  # 65%
    
    # Configuration triple effet
    n_effets = 3
    pressions = np.linspace(P_vapeur_vive, P_condenseur, n_effets)
    U_list = [2.5, 2.2, 1.8]  # kW/m²/K
    
    print("=== Simulation évaporation triple effet ===\n")
    
    # Simulation
    effets = systeme_multi_effets(
        F_alim=F_alim,
        X_alim=X_alim,
        T_alim=T_alim,
        pressions_effets=pressions,
        P_vapeur_vive=P_vapeur_vive,
        U_list=U_list,
        X_sortie_cible=X_sortie_cible,
        methode='sequentielle'
    )
    
    # Affichage résultats
    print(f"Configuration : {n_effets} effets")
    print(f"Alimentation : {F_alim:.1f} kg/h, {X_alim*100:.1f}%")
    print(f"Concentrat final : {effets[-1]['L_out']:.1f} kg/h, {effets[-1]['X_out']*100:.1f}%")
    print(f"Eau évaporée totale : {sum(e['V_out'] for e in effets):.1f} kg/h")
    print(f"Vapeur vive consommée : {effets[0]['m_vapeur_chaude']:.1f} kg/h")
    print(f"Économie de vapeur : {calculer_economy_globale(effets):.3f}")
    print(f"Surface totale : {sum(e['A_calc'] for e in effets):.1f} m²")
    
    # Visualisation
    visualiser_resultats(effets, "Évaporation triple effet - Résultats")
    
    # Optimisation nombre d'effets
    print("\n=== Optimisation nombre d'effets ===")
    resultats_opt = optimiser_nombre_effets(
        F_alim=F_alim,
        X_alim=X_alim,
        T_alim=T_alim,
        P_vapeur_vive=P_vapeur_vive,
        P_condenseur=P_condenseur,
        X_sortie_cible=X_sortie_cible,
        n_min=2,
        n_max=5
    )
    
    for n, res in resultats_opt.items():
        print(f"\n{n} effets :")
        print(f"  Surface totale : {res['surface_totale_m2']:.1f} m²")
        print(f"  Vapeur vive : {res['vapeur_vive_kg_h']:.1f} kg/h")
        print(f"  Économie : {res['economy']:.3f}")
        print(f"  Coût total annuel : {res['cout_total_annuel_euros']:.0f} €")
    
    # Analyse de sensibilité
    print("\n=== Analyse de sensibilité ===")
    resultats_sens = analyse_sensibilite(
        F_alim_base=F_alim,
        X_alim_base=X_alim,
        T_alim_base=T_alim,
        P_vapeur_base=P_vapeur_vive,
        P_condenseur_base=P_condenseur,
        X_sortie_base=X_sortie_cible,
        n_effets=3
    )
    
    return effets, resultats_opt, resultats_sens


if __name__ == "__main__":
    # Exécution de l'exemple
    effets, resultats_opt, resultats_sens = example_utilisation()