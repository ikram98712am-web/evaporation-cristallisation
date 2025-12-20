"""
thermodynamique.py
------------------
Module thermodynamique complet avec toutes les corrections selon l'énoncé.
- Corrélation de Dühring correcte
- Coefficient d'encrassement
- Propriétés solutions sucrées
- Validation complète
"""

import numpy as np
from CoolProp.CoolProp import PropsSI

R = 8.314  # J/mol/K

# ========= Eau pure =========

def enthalpie_eau_liquide(T):
    """
    Enthalpie eau liquide saturée à T [K].
    
    Parameters
    ----------
    T : float
        Température [K]
        
    Returns
    -------
    float
        Enthalpie [kJ/kg]
    """
    try:
        h = PropsSI("H", "T", T, "Q", 0, "Water")  # J/kg
        return h / 1000.0
    except Exception as e:
        raise ValueError(f"Erreur calcul enthalpie eau liquide à T={T}K: {e}")


def enthalpie_vapeur_saturee(P):
    """
    Enthalpie vapeur saturée à P [Pa].
    
    Parameters
    ----------
    P : float
        Pression [Pa]
        
    Returns
    -------
    float
        Enthalpie [kJ/kg]
    """
    try:
        h = PropsSI("H", "P", P, "Q", 1, "Water")  # J/kg
        return h / 1000.0
    except Exception as e:
        raise ValueError(f"Erreur calcul enthalpie vapeur à P={P}Pa: {e}")


def chaleur_latente_vaporisation(P):
    """
    Chaleur latente L à P [Pa].
    
    Parameters
    ----------
    P : float
        Pression [Pa]
        
    Returns
    -------
    float
        Chaleur latente [kJ/kg]
    """
    try:
        h_vap = PropsSI("H", "P", P, "Q", 1, "Water")
        h_liq = PropsSI("H", "P", P, "Q", 0, "Water")
        return (h_vap - h_liq) / 1000.0
    except Exception as e:
        raise ValueError(f"Erreur calcul chaleur latente à P={P}Pa: {e}")


def temperature_ebullition_eau(P):
    """
    T_eb de l'eau à P [Pa].
    
    Parameters
    ----------
    P : float
        Pression [Pa]
        
    Returns
    -------
    float
        Température d'ébullition [K]
    """
    try:
        return PropsSI("T", "P", P, "Q", 0.5, "Water")
    except Exception as e:
        raise ValueError(f"Erreur calcul T_eb eau à P={P}Pa: {e}")


# ========= EPE solution sucrée - CORRÉLATION DE DÜHRING CORRECTE =========

def elevation_point_ebullition(concentration_massique, P):
    """
    EPE(X) [K] pour une solution de saccharose.
    
    Corrélation de Dühring selon l'énoncé :
    - Si x < 50% : A = 0.03, B = 0.00015
    - Si x ≥ 50% : A = 0.045, B = 0.0003
    
    EPE = A*x + B*x² (avec x en %)
    
    Parameters
    ----------
    concentration_massique : float
        Fraction massique saccharose (0-1)
    P : float
        Pression [Pa] (non utilisée dans cette corrélation simplifiée)
        
    Returns
    -------
    float
        Élévation point d'ébullition [K]
        
    References
    ----------
    Guide aide_PIC_evapo.pdf, Section 2.2, équations (4) et (5)
    """
    if not 0 <= concentration_massique <= 1:
        raise ValueError(f"Concentration doit être entre 0 et 1, reçu: {concentration_massique}")
    
    x_pct = concentration_massique * 100.0  # Conversion en %
    
    # Coefficients selon la concentration
    if x_pct < 50.0:
        A = 0.03
        B = 0.00015
    else:
        A = 0.045
        B = 0.0003
    
    EPE = A * x_pct + B * x_pct**2
    return EPE  # en K (ou °C, c'est une différence)


def temperature_ebullition_solution(concentration_massique, P):
    """
    T_eb solution sucrée à P [Pa] : T_eau + EPE.
    
    Parameters
    ----------
    concentration_massique : float
        Fraction massique saccharose (0-1)
    P : float
        Pression [Pa]
        
    Returns
    -------
    float
        Température d'ébullition [K]
    """
    T_eau = temperature_ebullition_eau(P)
    EPE = elevation_point_ebullition(concentration_massique, P)
    return T_eau + EPE


# ========= Propriétés solution sucrée =========

def capacite_calorifique_solution(concentration_massique, T=None):
    """
    Capacité calorifique Cp de la solution [J/(kg·K)].
    
    Formule : Cp = (1-x)*Cp_eau + x*Cp_saccharose
    
    Parameters
    ----------
    concentration_massique : float
        Fraction massique saccharose (0-1)
    T : float, optional
        Température [K] (non utilisée dans version simplifiée)
        
    Returns
    -------
    float
        Capacité calorifique [J/(kg·K)]
    """
    Cp_eau = 4180.0  # J/(kg·K)
    Cp_saccharose = 1250.0  # J/(kg·K)
    
    x = concentration_massique
    return (1.0 - x) * Cp_eau + x * Cp_saccharose


def masse_volumique_solution(concentration_massique, T):
    """
    Masse volumique de la solution [kg/m³].
    
    Formule empirique : ρ = 1000 + 400*x - 0.3*T
    
    Parameters
    ----------
    concentration_massique : float
        Fraction massique saccharose (0-1)
    T : float
        Température [K]
        
    Returns
    -------
    float
        Masse volumique [kg/m³]
    """
    T_C = T - 273.15  # Conversion K -> °C
    x = concentration_massique
    return 1000.0 + 400.0 * x - 0.3 * T_C

def epe_saccharose_duhring(X, T_ref=None):
    """
    Élévation du point d'ébullition (EPE) pour solutions de saccharose.
    Corrélation de Dühring basée sur des données expérimentales.
    """
    if X <= 0:
        return 0.0
    
    # Limiter la concentration
    X = min(X, 0.75)
    
    # Coefficients pour la corrélation
    a = 0.42  # K
    b = 1.85  # K
    c = 2.30  # K
    
    # Calcul de l'EPE
    EPE = a * X + b * X**2 + c * X**3
    
    # Correction pour la température
    if T_ref is not None:
        T_ref_C = T_ref - 273.15
        f_temp = 1.0 - 0.002 * max(0, T_ref_C - 100.0)
        EPE *= f_temp
    
    return min(EPE, 20.0)
def enthalpie_solution_sucree(T, concentration_massique):
    """
    Enthalpie massique de la solution [kJ/kg].
    
    Approximation : mélange linéaire eau + sucre
    
    Parameters
    ----------
    T : float
        Température [K]
    concentration_massique : float
        Fraction massique saccharose (0-1)
        
    Returns
    -------
    float
        Enthalpie [kJ/kg]
    """
    X = concentration_massique
    h_eau = enthalpie_eau_liquide(T)  # kJ/kg d'eau
    
    # Enthalpie du saccharose (approximation simplifiée)
    T_C = T - 273.15
    h_sucre = 1.25 * T_C / 1000.0  # kJ/kg (Cp_saccharose ~ 1.25 kJ/kg/K)
    
    return X * h_sucre + (1.0 - X) * h_eau


# ========= Coefficient de transfert avec encrassement =========

def coefficient_transfert_effectif(U_propre, Rf=0.0002):
    """
    Coefficient global de transfert avec encrassement.
    
    Formule : 1/U_eff = 1/U_propre + Rf
    
    Parameters
    ----------
    U_propre : float
        Coefficient de transfert propre [W/(m²·K)]
    Rf : float, optional
        Résistance d'encrassement [m²·K/W], défaut = 0.0002
        
    Returns
    -------
    float
        Coefficient effectif [W/(m²·K)]
        
    References
    ----------
    Énoncé Section 4.1.1 : Rf = 0.0002 m²·K/W
    """
    if U_propre <= 0:
        raise ValueError(f"U_propre doit être > 0, reçu: {U_propre}")
    if Rf < 0:
        raise ValueError(f"Rf doit être >= 0, reçu: {Rf}")
    
    return 1.0 / (1.0 / U_propre + Rf)


# ========= Fonctions de validation =========

def valider_bilans_matiere(F_in, L_out, V_out, tolerance=0.01):
    """
    Vérifie que le bilan matière se ferme.
    
    Parameters
    ----------
    F_in : float
        Débit d'entrée [kg/h]
    L_out : float
        Débit liquide sortie [kg/h]
    V_out : float
        Débit vapeur sortie [kg/h]
    tolerance : float
        Tolérance relative (défaut 1%)
        
    Returns
    -------
    bool
        True si bilan fermé
        
    Raises
    ------
    ValueError
        Si erreur > tolérance
    """
    erreur_relative = abs(F_in - (L_out + V_out)) / F_in
    
    if erreur_relative > tolerance:
        raise ValueError(
            f"Bilan matière non fermé: F_in={F_in:.1f}, "
            f"L_out+V_out={L_out+V_out:.1f}, "
            f"erreur={erreur_relative*100:.2f}%"
        )
    
    return True


def valider_bilan_saccharose(F_in, X_in, L_out, X_out, tolerance=0.01):
    """
    Vérifie que le bilan saccharose se ferme.
    
    Parameters
    ----------
    F_in : float
        Débit d'entrée [kg/h]
    X_in : float
        Fraction massique entrée
    L_out : float
        Débit liquide sortie [kg/h]
    X_out : float
        Fraction massique sortie
    tolerance : float
        Tolérance relative (défaut 1%)
        
    Returns
    -------
    bool
        True si bilan fermé
    """
    saccharose_in = F_in * X_in
    saccharose_out = L_out * X_out
    erreur_relative = abs(saccharose_in - saccharose_out) / saccharose_in
    
    if erreur_relative > tolerance:
        raise ValueError(
            f"Bilan saccharose non fermé: "
            f"in={saccharose_in:.1f}, out={saccharose_out:.1f}, "
            f"erreur={erreur_relative*100:.2f}%"
        )
    
    return True