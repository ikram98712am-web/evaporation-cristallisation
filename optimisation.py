"""
optimisation_complete.py
------------------------
Analyse technico-économique complète pour le projet d'évaporation-cristallisation.
- CAPEX détaillé avec méthodes paramétriques
- OPEX complet (énergie, personnel, maintenance)
- VAN, TRI, ROI, période de récupération
- Optimisation nombre d'effets
- Intégration énergétique
- Analyse de sensibilité économique
"""

import numpy as np
from scipy.optimize import fsolve, minimize
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from evaporateurs import systeme_multi_effets, calculer_economy_globale
from cristallisation import dimensionner_cristalliseur_complet

# ============================================
# 1. PARAMÈTRES ÉCONOMIQUES (MAROC)
# ============================================

# Temps de fonctionnement
HEURES_FONCTIONNEMENT_AN = 8000.0  # h/an
JOURS_FONCTIONNEMENT_AN = 330.0    # jours/an

# Coûts énergétiques (MAD)
PRIX_GAZ_NATUREL_MAD_M3 = 5.50     # MAD/m³ gaz naturel
PCI_GAZ_NATUREL_MJ_M3 = 35.0       # MJ/m³
PRIX_ELECTRICITE_MAD_KWH = 1.20    # MAD/kWh
PRIX_EAU_INDUSTRIELLE_MAD_M3 = 15.0  # MAD/m³
PRIX_EAU_REFROIDISSEMENT_MAD_M3 = 1.5  # MAD/m³

# Rendements
RENDEMENT_CHAUDIERE = 0.85
RENDEMENT_POMPE = 0.75
RENDEMENT_MOTEUR = 0.92

# Coûts main d'œuvre
SALAIRE_OPERATEUR_MAD_H = 45.0     # MAD/h
SALAIRE_TECHNICIEN_MAD_H = 70.0    # MAD/h
SALAIRE_INGENIEUR_MAD_H = 120.0    # MAD/h
NOMBRE_OPERATEURS_POSTE = 2
NOMBRE_TECHNICIENS = 1
NOMBRE_INGENIEURS = 0.5  # mi-temps

# Coûts maintenance
MAINTENANCE_PREVENTIVE_FRACTION = 0.02  # 2% du CAPEX
MAINTENANCE_CORRECTIVE_FRACTION = 0.01  # 1% du CAPEX
CONSOMMABLES_FRACTION = 0.005           # 0.5% du CAPEX

# Inflation et taux
TAUX_INFLATION = 0.03    # 3% par an
TAUX_ACTUALISATION = 0.10  # 10%
TAU_IMPOSITION = 0.30    # 30%

# Durées
DUREE_AMORTISSEMENT = 10  # ans
DUREE_VIE_PROJET = 15     # ans
PERIODE_CONSTRUCTION = 1   # an

# Taux de change
TAUX_CHANGE_EUR_MAD = 10.5  # 1 EUR = 10.5 MAD
TAUX_CHANGE_USD_MAD = 9.5   # 1 USD = 9.5 MAD

# ============================================
# 2. FONCTIONS ÉCONOMIQUES DE BASE
# ============================================

def facteur_recuperation_capital(i, n):
    """
    Capital Recovery Factor (CRF).
    
    CRF = i(1+i)^n / [(1+i)^n - 1]
    """
    if i == 0:
        return 1.0 / n
    return i * (1 + i) ** n / ((1 + i) ** n - 1)


def facteur_valeur_actuelle(i, n):
    """
    Present Worth Factor (PWF).
    
    PWF = [(1+i)^n - 1] / [i(1+i)^n]
    """
    if i == 0:
        return n
    return ((1 + i) ** n - 1) / (i * (1 + i) ** n)


def cout_vapeur_thermique(P_vapeur, T_surchauffe=10.0):
    """
    Calcule le coût de production de vapeur.
    
    Parameters
    ----------
    P_vapeur : float
        Pression de vapeur [Pa]
    T_surchauffe : float
        Surchauffe [K]
        
    Returns
    -------
    float
        Coût [MAD/tonne vapeur]
    """
    # Approximation énergétique
    # Chaleur latente moyenne ~ 2200 kJ/kg
    # + chaleur sensible pour la surchauffe
    
    h_latent = 2200.0  # kJ/kg
    Cp_vapeur = 2.0    # kJ/kg·K
    delta_h = h_latent + Cp_vapeur * T_surchauffe
    
    # Énergie nécessaire
    energie_kWh_kg = delta_h / 3600.0
    
    # Consommation gaz
    rendement_total = RENDEMENT_CHAUDIERE * 0.9  # Pertes distribution
    energie_gaz_kWh_kg = energie_kWh_kg / rendement_total
    
    # Conversion en gaz naturel
    PCI_kWh_m3 = PCI_GAZ_NATUREL_MJ_M3 / 3.6  # MJ/m³ -> kWh/m³
    consommation_gaz_m3_kg = energie_gaz_kWh_kg / PCI_kWh_m3
    
    # Coût
    cout_par_kg = consommation_gaz_m3_kg * PRIX_GAZ_NATUREL_MAD_M3
    cout_par_tonne = cout_par_kg * 1000.0
    
    return cout_par_tonne


# ============================================
# 3. CAPEX - COÛTS D'INVESTISSEMENT
# ============================================

def cout_evaporateur_modulaire(A_totale, n_effets, P_max=3.5e5):
    """
    Coût des évaporateurs selon formule modulaire.
    
    C_evap = C_base × (A_totale)^α × f_n_effets × f_pression
    
    Parameters
    ----------
    A_totale : float
        Surface totale d'échange [m²]
    n_effets : int
        Nombre d'effets
    P_max : float
        Pression maximale [Pa]
        
    Returns
    -------
    float
        Coût en MAD
    """
    # Coût de base (référence 100 m², 3 effets, 3.5 bar)
    C_base_EUR = 500000.0
    C_base_MAD = C_base_EUR * TAUX_CHANGE_EUR_MAD
    
    # Facteur surface (économie d'échelle)
    alpha = 0.65
    facteur_surface = (A_totale / 100.0) ** alpha
    
    # Facteur nombre d'effets (coût supplémentaire par effet)
    facteur_n_effets = 1.0 + (n_effets - 3) * 0.15
    
    # Facteur pression (coût plus élevé pour pression élevée)
    P_bar = P_max / 1e5
    facteur_pression = 1.0 + 0.1 * (P_bar - 3.5)
    
    cout_MAD = C_base_MAD * facteur_surface * facteur_n_effets * facteur_pression
    
    return cout_MAD


def cout_cristalliseur_avec_agitation(V_batch, puissance_agitation_kW):
    """
    Coût du cristalliseur avec système d'agitation.
    
    Parameters
    ----------
    V_batch : float
        Volume utile [m³]
    puissance_agitation_kW : float
        Puissance d'agitation [kW]
        
    Returns
    -------
    float
        Coût en MAD
    """
    # Coût cuve
    C_cuve_base_EUR = 80000.0
    C_cuve_MAD = C_cuve_base_EUR * TAUX_CHANGE_EUR_MAD
    facteur_volume = (V_batch / 10.0) ** 0.6
    
    C_cuve = C_cuve_MAD * facteur_volume
    
    # Coût agitation
    C_agitation_base_EUR = 25000.0
    C_agitation_MAD = C_agitation_base_EUR * TAUX_CHANGE_EUR_MAD
    facteur_puissance = (puissance_agitation_kW / 10.0) ** 0.7
    
    C_agitation = C_agitation_MAD * facteur_puissance
    
    # Coût système de refroidissement
    C_refroidissement = 0.2 * C_cuve  # Estimation
    
    return C_cuve + C_agitation + C_refroidissement


def cout_echangeurs_chaleur(A_totale, type_echangeur='tube-calandre'):
    """
    Coût des échangeurs de chaleur.
    
    Parameters
    ----------
    A_totale : float
        Surface totale [m²]
    type_echangeur : str
        Type d'échangeur
        
    Returns
    -------
    float
        Coût en MAD
    """
    # Coûts selon type
    if type_echangeur == 'tube-calandre':
        C_base_EUR = 3000.0
        alpha = 0.7
    elif type_echangeur == 'plaque':
        C_base_EUR = 2500.0
        alpha = 0.65
    else:  # air-cooled
        C_base_EUR = 4000.0
        alpha = 0.75
    
    C_base_MAD = C_base_EUR * TAUX_CHANGE_EUR_MAD
    
    return C_base_MAD * (A_totale ** alpha)


def cout_instrumentation_controle(n_effets, A_totale):
    """
    Coût du système d'instrumentation et contrôle.
    """
    # Coût de base pour une unité standard
    C_base_EUR = 150000.0
    C_base_MAD = C_base_EUR * TAUX_CHANGE_EUR_MAD
    
    # Facteurs
    facteur_complexite = 1.0 + (n_effets - 2) * 0.1
    facteur_taille = (A_totale / 500.0) ** 0.3
    
    return C_base_MAD * facteur_complexite * facteur_taille


def cout_tuyauterie_isolation(n_effets, P_max=3.5e5):
    """
    Coût de la tuyauterie et isolation.
    """
    C_base_EUR = 80000.0
    C_base_MAD = C_base_EUR * TAUX_CHANGE_EUR_MAD
    
    P_bar = P_max / 1e5
    facteur_pression = 1.0 + 0.05 * (P_bar - 3.5)
    facteur_n_effets = 1.0 + (n_effets - 2) * 0.15
    
    return C_base_MAD * facteur_pression * facteur_n_effets


def calculer_capex_complet(A_evap, n_effets, V_crist, P_agitation, 
                          A_echangeurs=100.0, P_max=3.5e5):
    """
    Calcule le CAPEX complet (Total Capital Investment).
    
    Parameters
    ----------
    A_evap : float
        Surface évaporateurs [m²]
    n_effets : int
        Nombre d'effets
    V_crist : float
        Volume cristalliseur [m³]
    P_agitation : float
        Puissance agitation [kW]
    A_echangeurs : float
        Surface échangeurs [m²]
    P_max : float
        Pression maximale [Pa]
        
    Returns
    -------
    dict
        Détails du CAPEX
    """
    # Coûts équipements principaux (Battery Limits)
    C_evaporateurs = cout_evaporateur_modulaire(A_evap, n_effets, P_max)
    C_cristalliseur = cout_cristalliseur_avec_agitation(V_crist, P_agitation)
    C_echangeurs = cout_echangeurs_chaleur(A_echangeurs)
    
    C_equipement = C_evaporateurs + C_cristalliseur + C_echangeurs
    
    # Coûts auxiliaires
    C_instrumentation = cout_instrumentation_controle(n_effets, A_evap)
    C_tuyauterie = cout_tuyauterie_isolation(n_effets, P_max)
    C_pompes = 0.15 * C_equipement  # Estimation
    C_reservoirs = 0.10 * C_equipement  # Estimation
    C_electrique = 0.12 * C_equipement  # Estimation
    C_civil = 0.08 * C_equipement  # Estimation
    
    # Coûts indirects
    C_ingenierie = 0.20 * C_equipement  # 20%
    C_construction = 0.15 * C_equipement  # 15%
    C_gerance = 0.08 * C_equipement  # 8%
    C_essais = 0.05 * C_equipement  # 5%
    
    # Calculs totaux
    C_direct = (C_equipement + C_instrumentation + C_tuyauterie + 
                C_pompes + C_reservoirs + C_electrique + C_civil)
    
    C_indirect = C_ingenierie + C_construction + C_gerance + C_essais
    
    C_total_installe = C_direct + C_indirect
    
    # Capital de démarrage
    C_demarrage = 0.10 * C_total_installe
    
    # Fonds de roulement
    fonds_roulement = 0.15 * C_total_installe
    
    # Total Capital Investment (TCI)
    TCI = C_total_installe + C_demarrage + fonds_roulement
    
    return {
        'equipements_principaux': {
            'evaporateurs': C_evaporateurs,
            'cristalliseur': C_cristalliseur,
            'echangeurs': C_echangeurs,
            'total': C_equipement
        },
        'equipements_auxiliaires': {
            'instrumentation': C_instrumentation,
            'tuyauterie': C_tuyauterie,
            'pompes': C_pompes,
            'reservoirs': C_reservoirs,
            'electrique': C_electrique,
            'civil': C_civil,
            'total': C_direct - C_equipement
        },
        'cout_direct': C_direct,
        'couts_indirects': {
            'ingenierie': C_ingenierie,
            'construction': C_construction,
            'gerance': C_gerance,
            'essais': C_essais,
            'total': C_indirect
        },
        'total_installe': C_total_installe,
        'cout_demarrage': C_demarrage,
        'fonds_roulement': fonds_roulement,
        'TCI': TCI
    }


# ============================================
# 4. OPEX - COÛTS D'EXPLOITATION
# ============================================

def calculer_opex_annuel(vapeur_kg_h, electricite_kW, eau_m3_h, 
                        personnel_details=None, capex_details=None):
    """
    Calcule les coûts d'exploitation annuels.
    
    Parameters
    ----------
    vapeur_kg_h : float
        Consommation vapeur [kg/h]
    electricite_kW : float
        Puissance électrique [kW]
    eau_m3_h : float
        Consommation eau [m³/h]
    personnel_details : dict, optional
        Détails du personnel
    capex_details : dict, optional
        Détails du CAPEX pour calculer maintenance
        
    Returns
    -------
    dict
        Détails OPEX annuel
    """
    if personnel_details is None:
        personnel_details = {
            'operateurs': NOMBRE_OPERATEURS_POSTE,
            'techniciens': NOMBRE_TECHNICIENS,
            'ingenieurs': NOMBRE_INGENIEURS
        }
    
    # 1. Énergie - Vapeur
    vapeur_tonne_an = (vapeur_kg_h / 1000.0) * HEURES_FONCTIONNEMENT_AN
    cout_vapeur_par_tonne = cout_vapeur_thermique(3.5e5)  # Coût moyen
    cout_vapeur_an = vapeur_tonne_an * cout_vapeur_par_tonne
    
    # 2. Énergie - Électricité
    electricite_kWh_an = electricite_kW * HEURES_FONCTIONNEMENT_AN
    cout_elec_an = electricite_kWh_an * PRIX_ELECTRICITE_MAD_KWH
    
    # 3. Eau
    eau_m3_an = eau_m3_h * HEURES_FONCTIONNEMENT_AN
    cout_eau_proc_an = eau_m3_an * 0.5 * PRIX_EAU_INDUSTRIELLE_MAD_M3  # 50% eau procédé
    cout_eau_refroid_an = eau_m3_an * 0.5 * PRIX_EAU_REFROIDISSEMENT_MAD_M3  # 50% refroidissement
    
    cout_eau_total_an = cout_eau_proc_an + cout_eau_refroid_an
    
    # 4. Personnel
    heures_an = HEURES_FONCTIONNEMENT_AN
    
    cout_operateurs = (personnel_details['operateurs'] * 3 * heures_an * 
                       SALAIRE_OPERATEUR_MAD_H)
    cout_techniciens = (personnel_details['techniciens'] * heures_an * 
                        SALAIRE_TECHNICIEN_MAD_H)
    cout_ingenieurs = (personnel_details['ingenieurs'] * heures_an * 
                       SALAIRE_INGENIEUR_MAD_H)
    
    cout_personnel_an = cout_operateurs + cout_techniciens + cout_ingenieurs
    
    # 5. Maintenance
    cout_maintenance_an = 0.0
    if capex_details:
        TCI = capex_details.get('TCI', 0)
        cout_maintenance_an = (MAINTENANCE_PREVENTIVE_FRACTION + 
                              MAINTENANCE_CORRECTIVE_FRACTION + 
                              CONSOMMABLES_FRACTION) * TCI
    
    # 6. Autres
    cout_assurance_an = 0.005 * (capex_details['TCI'] if capex_details else 1e6)
    cout_analyses_an = 50000.0  # MAD/an (fixe)
    cout_divers_an = 30000.0    # MAD/an (fixe)
    
    # Total OPEX
    OPEX_total = (cout_vapeur_an + cout_elec_an + cout_eau_total_an + 
                  cout_personnel_an + cout_maintenance_an + cout_assurance_an + 
                  cout_analyses_an + cout_divers_an)
    
    return {
        'energie': {
            'vapeur': cout_vapeur_an,
            'electricite': cout_elec_an,
            'eau': cout_eau_total_an,
            'total': cout_vapeur_an + cout_elec_an + cout_eau_total_an
        },
        'personnel': {
            'operateurs': cout_operateurs,
            'techniciens': cout_techniciens,
            'ingenieurs': cout_ingenieurs,
            'total': cout_personnel_an
        },
        'maintenance': cout_maintenance_an,
        'autres': {
            'assurance': cout_assurance_an,
            'analyses': cout_analyses_an,
            'divers': cout_divers_an,
            'total': cout_assurance_an + cout_analyses_an + cout_divers_an
        },
        'OPEX_total': OPEX_total,
        'OPEX_par_heure': OPEX_total / HEURES_FONCTIONNEMENT_AN
    }


# ============================================
# 5. REVENUS ET RENTABILITÉ
# ============================================

def calculer_revenus(sucre_kg_h, prix_sucre_MAD_kg=8.0,  # 8 MAD/kg = 8000 MAD/tonne
                     sous_produits_MAD_h=0.0):
    """
    Calcule les revenus annuels.
    
    Parameters
    ----------
    sucre_kg_h : float
        Production de sucre [kg/h]
    prix_sucre_MAD_kg : float
        Prix de vente du sucre [MAD/kg]
    sous_produits_MAD_h : float
        Valeur des sous-produits [MAD/h]
        
    Returns
    -------
    dict
        Détails des revenus
    """
    # Production annuelle
    sucre_tonne_an = (sucre_kg_h / 1000.0) * HEURES_FONCTIONNEMENT_AN
    
    # Revenus sucre
    revenus_sucre_an = sucre_tonne_an * prix_sucre_MAD_kg * 1000  # Conversion kg->tonne
    
    # Revenus sous-produits
    revenus_sous_produits_an = sous_produits_MAD_h * HEURES_FONCTIONNEMENT_AN
    
    # Total revenus
    revenus_totaux_an = revenus_sucre_an + revenus_sous_produits_an
    
    return {
        'production_sucre_tonne_an': sucre_tonne_an,
        'prix_sucre_MAD_tonne': prix_sucre_MAD_kg * 1000,
        'revenus_sucre_an': revenus_sucre_an,
        'revenus_sous_produits_an': revenus_sous_produits_an,
        'revenus_totaux_an': revenus_totaux_an
    }


def calculer_flux_tresorerie(capex_details, opex_annuel, revenus_annuels, 
                            duree_projet=DUREE_VIE_PROJET):
    """
    Calcule les flux de trésorerie sur la durée du projet.
    
    Parameters
    ----------
    capex_details : dict
        Détails du CAPEX
    opex_annuel : dict
        OPEX annuel
    revenus_annuels : dict
        Revenus annuels
    duree_projet : int
        Durée du projet [ans]
        
    Returns
    -------
    dict
        Flux de trésorerie
    """
    TCI = capex_details['TCI']
    OPEX_total = opex_annuel['OPEX_total']
    revenus_totaux = revenus_annuels['revenus_totaux_an']
    
    # Flux annuels
    flux_bruts = []
    flux_nets = []
    flux_actualises = []
    
    # Année 0 : Investissement
    flux_bruts.append(-TCI)
    flux_nets.append(-TCI)
    flux_actualises.append(-TCI)
    
    # Années d'exploitation
    for annee in range(1, duree_projet + 1):
        # Revenus bruts
        revenus_bruts = revenus_totaux * (1 + TAUX_INFLATION) ** (annee - 1)
        
        # OPEX avec inflation
        opex_annee = OPEX_total * (1 + TAUX_INFLATION) ** (annee - 1)
        
        # EBITDA
        EBITDA = revenus_bruts - opex_annee
        
        # Amortissement
        amortissement = TCI / DUREE_AMORTISSEMENT
        
        # Résultat avant impôt
        resultat_avant_impot = EBITDA - amortissement
        
        # Impôt
        impot = max(0, resultat_avant_impot) * TAU_IMPOSITION
        
        # Résultat net
        resultat_net = resultat_avant_impot - impot
        
        # Flux de trésorerie
        flux_tresorerie = resultat_net + amortissement
        
        flux_bruts.append(flux_tresorerie)
        flux_nets.append(resultat_net)
        flux_actualises.append(flux_tresorerie / (1 + TAUX_ACTUALISATION) ** annee)
    
    # Calcul des cumuls
    cumul_brut = np.cumsum(flux_bruts)
    cumul_net = np.cumsum(flux_nets)
    cumul_actualise = np.cumsum(flux_actualises)
    
    return {
        'annees': list(range(duree_projet + 1)),
        'flux_bruts': flux_bruts,
        'flux_nets': flux_nets,
        'flux_actualises': flux_actualises,
        'cumul_brut': cumul_brut,
        'cumul_net': cumul_net,
        'cumul_actualise': cumul_actualise,
        'EBITDA': revenus_totaux - OPEX_total,
        'amortissement_annuel': TCI / DUREE_AMORTISSEMENT
    }


# ============================================
# 6. INDICATEURS DE RENTABILITÉ
# ============================================

def calculer_VAN(flux_actualises):
    """
    Calcule la Valeur Actuelle Nette (VAN).
    """
    return np.sum(flux_actualises)


def calculer_TRI(flux_bruts, precision=1e-6, max_iter=1000):
    """
    Calcule le Taux de Rendement Interne (TRI) par méthode itérative.
    """
    def van_taux(taux):
        van = 0
        for t, flux in enumerate(flux_bruts):
            van += flux / (1 + taux) ** t
        return van
    
    # Recherche par dichotomie
    taux_min, taux_max = -0.5, 2.0
    van_min = van_taux(taux_min)
    van_max = van_taux(taux_max)
    
    for _ in range(max_iter):
        taux_moyen = (taux_min + taux_max) / 2
        van_moyen = van_taux(taux_moyen)
        
        if abs(van_moyen) < precision:
            return taux_moyen
        
        if van_min * van_moyen < 0:
            taux_max = taux_moyen
            van_max = van_moyen
        else:
            taux_min = taux_moyen
            van_min = van_moyen
    
    return (taux_min + taux_max) / 2


def calculer_periode_recuperation(flux_tresorerie, investissement_initial):
    """
    Calcule la période de récupération du capital.
    
    Returns
    -------
    float
        Période de récupération [ans]
    """
    cumul = 0
    for annee, flux in enumerate(flux_tresorerie):
        if annee == 0:
            continue
        cumul += flux
        if cumul >= investissement_initial:
            # Interpolation pour plus de précision
            flux_annee = flux
            manquant = investissement_initial - (cumul - flux_annee)
            fraction_annee = manquant / flux_annee if flux_annee != 0 else 1
            return annee - 1 + fraction_annee
    return np.inf


def calculer_indice_rentabilite(VAN, investissement_initial):
    """
    Calcule l'Indice de Rentabilité (IR).
    
    IR = VAN / Investissement initial
    """
    if investissement_initial <= 0:
        return np.inf
    return VAN / investissement_initial


def calculer_BEP(cout_fixe_annuel, marge_contrib_unitaire, prix_vente_unitaire):
    """
    Calcule le Break-Even Point (BEP).
    
    BEP (quantité) = Coûts fixes / (Prix - Coût variable unitaire)
    BEP (CA) = Coûts fixes / Marge sur coûts variables
    """
    if marge_contrib_unitaire <= 0:
        return np.inf
    
    bep_quantite = cout_fixe_annuel / marge_contrib_unitaire
    bep_ca = bep_quantite * prix_vente_unitaire
    
    return {
        'bep_quantite_tonne_an': bep_quantite,
        'bep_chiffre_affaires_MAD_an': bep_ca
    }


def analyse_rentabilite_complete(capex_details, opex_annuel, revenus_annuels, 
                                flux_tresorerie):
    """
    Analyse complète de la rentabilité.
    """
    TCI = capex_details['TCI']
    
    # VAN
    VAN = calculer_VAN(flux_tresorerie['flux_actualises'])
    
    # TRI
    TRI = calculer_TRI(flux_tresorerie['flux_bruts'])
    
    # Période de récupération
    periode_recup = calculer_periode_recuperation(flux_tresorerie['flux_bruts'][1:], TCI)
    
    # Indice de rentabilité
    IR = calculer_indice_rentabilite(VAN, TCI)
    
    # BEP
    cout_fixe_annuel = opex_annuel['OPEX_total']
    production_annuelle = revenus_annuels['production_sucre_tonne_an']
    prix_tonne = revenus_annuels['prix_sucre_MAD_tonne']
    
    # Estimation coût variable unitaire (60% des coûts variables totaux)
    cout_variable_tonne = (opex_annuel['energie']['total'] * 0.6 + 
                          opex_annuel['autres']['total']) / production_annuelle
    marge_unitaire = prix_tonne - cout_variable_tonne
    
    BEP = calculer_BEP(cout_fixe_annuel, marge_unitaire, prix_tonne)
    
    # Taux de rendement comptable
    profit_annuel_moyen = np.mean(flux_tresorerie['flux_nets'][1:])
    TRC = profit_annuel_moyen / TCI if TCI > 0 else 0
    
    return {
        'VAN_MAD': VAN,
        'TRI': TRI,
        'periode_recuperation_ans': periode_recup,
        'indice_rentabilite': IR,
        'TRC': TRC,
        'BEP': BEP,
        'cout_variable_unitaire_MAD_tonne': cout_variable_tonne,
        'marge_unitaire_MAD_tonne': marge_unitaire
    }


# ============================================
# 7. OPTIMISATION DU NOMBRE D'EFFETS
# ============================================

def optimiser_nombre_effets_avec_economique(
    F_alim=20000.0,
    X_alim=0.15,
    T_alim=85.0 + 273.15,
    P_vapeur_vive=3.5e5,
    P_condenseur=0.15e5,
    X_sortie_cible=0.65,
    V_cristalliseur=10.0,
    puissance_agitation=50.0
):
    """
    Optimise le nombre d'effets avec analyse économique complète.
    
    Returns
    -------
    list
        Résultats pour 2, 3, 4, 5 effets
    """
    resultats = []
    nombres_effets = [2, 3, 4, 5]
    
    for n in nombres_effets:
        print(f"Analyse pour {n} effets...")
        
        # Distribution des pressions et coefficients U
        pressions = np.linspace(P_vapeur_vive, P_condenseur, n)
        U_base = np.linspace(2.5, 1.8, n)  # kW/m²·K
        
        # Simulation évaporation
        effets = systeme_multi_effets(
            F_alim=F_alim,
            X_alim=X_alim,
            T_alim=T_alim,
            pressions_effets=pressions,
            P_vapeur_vive=P_vapeur_vive,
            U_list=U_base,
            X_sortie_cible=X_sortie_cible,
            methode='sequentielle'
        )
        
        # Calcul des indicateurs techniques
        surface_totale = sum(e.get('A_calc', 0) for e in effets)
        vapeur_vive = effets[0].get('m_vapeur_chaude', 0) if effets else 0
        eau_evaporee = sum(e.get('V_out', 0) for e in effets)
        economy = calculer_economy_globale(effets) if effets else 0
        
        # Production sucre (simplifiée)
        dernier_effet = effets[-1] if effets else {}
        production_sucre_kg_h = dernier_effet.get('L_out', 0) * X_sortie_cible
        
        # CAPEX
        capex_details = calculer_capex_complet(
            A_evap=surface_totale,
            n_effets=n,
            V_crist=V_cristalliseur,
            P_agitation=puissance_agitation,
            A_echangeurs=surface_totale * 0.5,  # Estimation
            P_max=P_vapeur_vive
        )
        
        # OPEX
        # Estimation consommation eau (condenseur dernier effet)
        eau_condenseur_m3_h = dernier_effet.get('V_out', 0) / 1000.0 * 30
        
        # Estimation puissance électrique
        puissance_elec_kW = 50.0 + n * 10.0  # Augmente avec le nombre d'effets
        
        opex_details = calculer_opex_annuel(
            vapeur_kg_h=vapeur_vive,
            electricite_kW=puissance_elec_kW,
            eau_m3_h=eau_condenseur_m3_h,
            capex_details=capex_details
        )
        
        # Revenus
        revenus_details = calculer_revenus(
            sucre_kg_h=production_sucre_kg_h,
            prix_sucre_MAD_kg=8.0
        )
        
        # Flux de trésorerie
        flux_details = calculer_flux_tresorerie(
            capex_details=capex_details,
            opex_annuel=opex_details,
            revenus_annuels=revenus_details
        )
        
        # Analyse de rentabilité
        rentabilite = analyse_rentabilite_complete(
            capex_details=capex_details,
            opex_annuel=opex_details,
            revenus_annuels=revenus_details,
            flux_tresorerie=flux_details
        )
        
        # Coût de l'eau évaporée
        cout_eau_evaporee_MAD_tonne = (opex_details['OPEX_total'] / 
                                      (eau_evaporee * HEURES_FONCTIONNEMENT_AN / 1000))
        
        resultats.append({
            'n_effets': n,
            'technique': {
                'surface_totale_m2': surface_totale,
                'vapeur_vive_kg_h': vapeur_vive,
                'eau_evaporee_kg_h': eau_evaporee,
                'economy': economy,
                'production_sucre_kg_h': production_sucre_kg_h
            },
            'capex': capex_details,
            'opex': opex_details,
            'revenus': revenus_details,
            'flux_tresorerie': flux_details,
            'rentabilite': rentabilite,
            'cout_eau_evaporee_MAD_tonne': cout_eau_evaporee_MAD_tonne,
            'cout_total_annuel_MAD': (capex_details['TCI'] * 
                                     facteur_recuperation_capital(TAUX_ACTUALISATION, DUREE_AMORTISSEMENT) + 
                                     opex_details['OPEX_total'])
        })
    
    return resultats


# ============================================
# 8. VISUALISATION DES RÉSULTATS
# ============================================

def visualiser_comparaison_n_effets(resultats_optimisation):
    """
    Visualise la comparaison du nombre d'effets.
    """
    n_vals = [r['n_effets'] for r in resultats_optimisation]
    
    fig = plt.figure(figsize=(16, 12))
    gs = GridSpec(3, 3, figure=fig, hspace=0.3, wspace=0.3)
    
    # 1. CAPEX vs OPEX
    ax1 = fig.add_subplot(gs[0, 0])
    capex_vals = [r['capex']['TCI'] / 1e6 for r in resultats_optimisation]  # Millions MAD
    opex_vals = [r['opex']['OPEX_total'] / 1e6 for r in resultats_optimisation]
    
    ax1.plot(n_vals, capex_vals, 'bo-', markersize=8, linewidth=2, label='CAPEX')
    ax1.plot(n_vals, opex_vals, 'rs--', markersize=8, linewidth=2, label='OPEX')
    ax1.set_xlabel('Nombre d\'effets')
    ax1.set_ylabel('Coût [Million MAD]')
    ax1.set_title('CAPEX vs OPEX')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. Consommation vapeur et économie
    ax2 = fig.add_subplot(gs[0, 1])
    vapeur_vals = [r['technique']['vapeur_vive_kg_h'] / 1000 for r in resultats_optimisation]  # tonne/h
    economy_vals = [r['technique']['economy'] for r in resultats_optimisation]
    
    ax2_1 = ax2.twinx()
    bars = ax2.bar(n_vals, vapeur_vals, alpha=0.7, color='lightblue')
    ax2.set_xlabel('Nombre d\'effets')
    ax2.set_ylabel('Vapeur vive [tonne/h]', color='blue')
    ax2.tick_params(axis='y', labelcolor='blue')
    
    ax2_1.plot(n_vals, economy_vals, 'g^-', linewidth=2, markersize=8)
    ax2_1.set_ylabel('Économie de vapeur [-]', color='green')
    ax2_1.tick_params(axis='y', labelcolor='green')
    
    ax2.set_title('Consommation et économie de vapeur')
    ax2.grid(True, alpha=0.3, axis='x')
    
    # 3. Indicateurs de rentabilité
    ax3 = fig.add_subplot(gs[0, 2])
    van_vals = [r['rentabilite']['VAN_MAD'] / 1e6 for r in resultats_optimisation]
    tri_vals = [r['rentabilite']['TRI'] * 100 for r in resultats_optimisation]  # %
    
    ax3_1 = ax3.twinx()
    bars_van = ax3.bar(n_vals, van_vals, alpha=0.7, color='orange')
    ax3.set_xlabel('Nombre d\'effets')
    ax3.set_ylabel('VAN [Million MAD]', color='orange')
    ax3.tick_params(axis='y', labelcolor='orange')
    
    ax3_1.plot(n_vals, tri_vals, 'm^-', linewidth=2, markersize=8)
    ax3_1.set_ylabel('TRI [%]', color='purple')
    ax3_1.tick_params(axis='y', labelcolor='purple')
    
    ax3.set_title('VAN et TRI')
    ax3.grid(True, alpha=0.3, axis='x')
    
    # 4. Coût de l'eau évaporée
    ax4 = fig.add_subplot(gs[1, 0])
    cout_eau_vals = [r['cout_eau_evaporee_MAD_tonne'] for r in resultats_optimisation]
    
    ax4.plot(n_vals, cout_eau_vals, 'ko-', linewidth=2, markersize=8)
    ax4.set_xlabel('Nombre d\'effets')
    ax4.set_ylabel('Coût eau évaporée [MAD/tonne]')
    ax4.set_title('Coût spécifique de l\'eau évaporée')
    ax4.grid(True, alpha=0.3)
    
    # 5. Production et revenus
    ax5 = fig.add_subplot(gs[1, 1])
    production_vals = [r['revenus']['production_sucre_tonne_an'] for r in resultats_optimisation]
    revenus_vals = [r['revenus']['revenus_totaux_an'] / 1e6 for r in resultats_optimisation]  # Millions
    
    ax5_1 = ax5.twinx()
    ax5.plot(n_vals, production_vals, 'bo-', linewidth=2, markersize=8, label='Production')
    ax5.set_xlabel('Nombre d\'effets')
    ax5.set_ylabel('Production sucre [tonne/an]', color='blue')
    ax5.tick_params(axis='y', labelcolor='blue')
    
    ax5_1.plot(n_vals, revenus_vals, 'rs--', linewidth=2, markersize=8, label='Revenus')
    ax5_1.set_ylabel('Revenus [Million MAD/an]', color='red')
    ax5_1.tick_params(axis='y', labelcolor='red')
    
    ax5.set_title('Production et revenus')
    ax5.grid(True, alpha=0.3)
    
    # 6. Tableau comparatif
    ax6 = fig.add_subplot(gs[1, 2])
    ax6.axis('tight')
    ax6.axis('off')
    
    tableau_data = [['Effets', 'VAN (M MAD)', 'TRI (%)', 'Période Récup. (ans)', 'Coût eau (MAD/t)']]
    
    for r in resultats_optimisation:
        tableau_data.append([
            str(r['n_effets']),
            f"{r['rentabilite']['VAN_MAD']/1e6:.1f}",
            f"{r['rentabilite']['TRI']*100:.1f}",
            f"{r['rentabilite']['periode_recuperation_ans']:.1f}",
            f"{r['cout_eau_evaporee_MAD_tonne']:.0f}"
        ])
    
    tableau = ax6.table(cellText=tableau_data, loc='center', cellLoc='center')
    tableau.auto_set_font_size(False)
    tableau.set_fontsize(9)
    tableau.scale(1, 1.5)
    
    # 7. Flux de trésorerie cumulés
    ax7 = fig.add_subplot(gs[2, :])
    
    for i, r in enumerate(resultats_optimisation):
        n = r['n_effets']
        cumul = r['flux_tresorerie']['cumul_actualise'] / 1e6
        annees = r['flux_tresorerie']['annees']
        
        ax7.plot(annees, cumul, 'o-', linewidth=2, markersize=4, label=f'{n} effets')
    
    ax7.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    ax7.set_xlabel('Années')
    ax7.set_ylabel('VAN cumulée [Million MAD]')
    ax7.set_title('Flux de trésorerie actualisés cumulés')
    ax7.legend()
    ax7.grid(True, alpha=0.3)
    
    plt.suptitle('Optimisation du nombre d\'effets - Analyse économique', fontsize=14, y=0.98)
    plt.tight_layout()
    plt.show()


# ============================================
# 9. FONCTION PRINCIPALE
# ============================================

def demo_optimisation_complete():
    """
    Démonstration complète de l'analyse d'optimisation.
    """
    print("=" * 80)
    print("ANALYSE D'OPTIMISATION ÉCONOMIQUE COMPLÈTE")
    print("=" * 80)
    
    print("\n1. OPTIMISATION DU NOMBRE D'EFFETS")
    print("-" * 40)
    
    resultats = optimiser_nombre_effets_avec_economique(
        F_alim=20000.0,
        X_alim=0.15,
        T_alim=85.0 + 273.15,
        P_vapeur_vive=3.5e5,
        P_condenseur=0.15e5,
        X_sortie_cible=0.65,
        V_crystalliseur=10.0,
        puissance_agitation=50.0
    )
    
    print("\nRésultats de l'optimisation :")
    print("-" * 60)
    print(f"{'Effets':<8} {'VAN (M MAD)':<12} {'TRI (%)':<10} {'Récup. (ans)':<12} {'Coût eau (MAD/t)':<15}")
    print("-" * 60)
    
    for r in resultats:
        print(f"{r['n_effets']:<8} {r['rentabilite']['VAN_MAD']/1e6:<12.1f} "
              f"{r['rentabilite']['TRI']*100:<10.1f} "
              f"{r['rentabilite']['periode_recuperation_ans']:<12.1f} "
              f"{r['cout_eau_evaporee_MAD_tonne']:<15.0f}")
    
    print("-" * 60)
    
    # Détermination de la configuration optimale
    van_vals = [r['rentabilite']['VAN_MAD'] for r in resultats]
    idx_optimal = np.argmax(van_vals)
    optimal = resultats[idx_optimal]
    
    print(f"\n✓ CONFIGURATION OPTIMALE : {optimal['n_effets']} effets")
    print(f"  VAN : {optimal['rentabilite']['VAN_MAD']/1e6:.1f} M MAD")
    print(f"  TRI : {optimal['rentabilite']['TRI']*100:.1f} %")
    print(f"  Période de récupération : {optimal['rentabilite']['periode_recuperation_ans']:.1f} ans")
    print(f"  Coût de l'eau évaporée : {optimal['cout_eau_evaporee_MAD_tonne']:.0f} MAD/tonne")
    
    # Visualisation
    print("\n2. GÉNÉRATION DES VISUALISATIONS")
    print("-" * 40)
    
    visualiser_comparaison_n_effets(resultats)
    
    # Analyse de sensibilité supplémentaire
    print("\n3. ANALYSE DE SENSIBILITÉ (exemple)")
    print("-" * 40)
    
    # Impact du prix du sucre
    prix_sucre_scenarios = [6.0, 8.0, 10.0, 12.0]  # MAD/kg
    van_sensibilite = []
    
    for prix in prix_sucre_scenarios:
        # Réutilisation des résultats avec nouveau prix
        r_optimal = optimal.copy()
        # Recalcul simplifié
        revenus_nouveau = r_optimal['revenus']['production_sucre_tonne_an'] * prix * 1000
        flux_modifies = r_optimal['flux_tresorerie']['flux_bruts'].copy()
        # Modification approximative (simplifiée)
        flux_modifies[1:] = [f * (prix/8.0) for f in flux_modifies[1:]]
        van_nouveau = calculer_VAN(flux_modifies) / (1 + TAUX_ACTUALISATION)
        van_sensibilite.append(van_nouveau)
    
    print(f"Impact du prix du sucre sur la VAN (configuration {optimal['n_effets']} effets):")
    for prix, van in zip(prix_sucre_scenarios, van_sensibilite):
        print(f"  Prix {prix:.1f} MAD/kg : VAN = {van/1e6:.1f} M MAD")
    
    print("\n" + "=" * 80)
    print("ANALYSE TERMINÉE AVEC SUCCÈS")
    print("=" * 80)
    
    return resultats, optimal


# ============================================
# 10. EXÉCUTION
# ============================================

if __name__ == "__main__":
    # Exécute l'analyse d'optimisation
    resultats, optimal = demo_optimisation_complete()
    
    # Sauvegarde des résultats
    import pickle
    with open('resultats_optimisation.pkl', 'wb') as f:
        pickle.dump({'resultats': resultats, 'optimal': optimal}, f)
    
    print("\nRésultats sauvegardés dans 'resultats_optimisation.pkl'")
    
    # Génération d'un rapport texte
    with open('rapport_optimisation.txt', 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("RAPPORT D'OPTIMISATION ÉCONOMIQUE\n")
        f.write("=" * 80 + "\n\n")
        
        f.write("CONFIGURATION OPTIMALE :\n")
        f.write("-" * 40 + "\n")
        f.write(f"Nombre d'effets : {optimal['n_effets']}\n")
        f.write(f"VAN : {optimal['rentabilite']['VAN_MAD']/1e6:.1f} M MAD\n")
        f.write(f"TRI : {optimal['rentabilite']['TRI']*100:.1f} %\n")
        f.write(f"Période de récupération : {optimal['rentabilite']['periode_recuperation_ans']:.1f} ans\n")
        f.write(f"Coût eau évaporée : {optimal['cout_eau_evaporee_MAD_tonne']:.0f} MAD/tonne\n\n")
        
        f.write("DÉTAILS TECHNIQUES :\n")
        f.write("-" * 40 + "\n")
        tech = optimal['technique']
        f.write(f"Surface totale : {tech['surface_totale_m2']:.1f} m²\n")
        f.write(f"Vapeur vive : {tech['vapeur_vive_kg_h']:.0f} kg/h\n")
        f.write(f"Eau évaporée : {tech['eau_evaporee_kg_h']:.0f} kg/h\n")
        f.write(f"Économie de vapeur : {tech['economy']:.2f}\n")
        f.write(f"Production sucre : {tech['production_sucre_kg_h']:.0f} kg/h\n\n")
        
        f.write("INVESTISSEMENT (CAPEX) :\n")
        f.write("-" * 40 + "\n")
        capex = optimal['capex']
        f.write(f"Total Capital Investment : {capex['TCI']/1e6:.1f} M MAD\n")
        f.write(f"Équipements principaux : {capex['equipements_principaux']['total']/1e6:.1f} M MAD\n")
        f.write(f"Équipements auxiliaires : {(capex['cout_direct'] - capex['equipements_principaux']['total'])/1e6:.1f} M MAD\n")
        f.write(f"Coûts indirects : {capex['couts_indirects']['total']/1e6:.1f} M MAD\n\n")
        
        f.write("EXPLOITATION (OPEX) :\n")
        f.write("-" * 40 + "\n")
        opex = optimal['opex']
        f.write(f"Total OPEX annuel : {opex['OPEX_total']/1e6:.1f} M MAD/an\n")
        f.write(f"  Énergie : {opex['energie']['total']/1e6:.2f} M MAD\n")
        f.write(f"  Personnel : {opex['personnel']['total']/1e6:.2f} M MAD\n")
        f.write(f"  Maintenance : {opex['maintenance']/1e6:.2f} M MAD\n")
        f.write(f"  Autres : {opex['autres']['total']/1e6:.2f} M MAD\n\n")
        
        f.write("REVENUS :\n")
        f.write("-" * 40 + "\n")
        revenus = optimal['revenus']
        f.write(f"Production annuelle : {revenus['production_sucre_tonne_an']:.0f} tonnes/an\n")
        f.write(f"Prix de vente : {revenus['prix_sucre_MAD_tonne']:.0f} MAD/tonne\n")
        f.write(f"Revenus annuels : {revenus['revenus_totaux_an']/1e6:.1f} M MAD/an\n")
    
    print("Rapport généré dans 'rapport_optimisation.txt'")