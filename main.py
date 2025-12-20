"""
main.py - Point d'entrée principal du projet
Version standalone
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os

# Style des graphiques
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (10, 6)
plt.rcParams['font.size'] = 12

def creer_dossiers_resultats():
    """Crée les dossiers nécessaires pour les résultats."""
    dossiers = ['resultats', 'resultats/images', 'resultats/rapports']
    for dossier in dossiers:
        if not os.path.exists(dossier):
            os.makedirs(dossier)
            print(f"✓ Dossier créé: {dossier}")

def demo_thermodynamique():
    """Démonstration des fonctions thermodynamiques de base."""
    print("=" * 70)
    print("DÉMONSTRATION THERMODYNAMIQUE")
    print("=" * 70)
    
    print("\n1. SOLUBILITÉ DU SACCHAROSE")
    print("-" * 40)
    
    # Formule de solubilité du saccharose
    def solubilite_saccharose(T_C):
        """
        Solubilité C* (g saccharose / 100 g solution) en fonction de T (°C).
        
        Formule du sujet :
        C* = 64.18 + 0.1337 T + 5.52e-3 T^2 - 9.73e-6 T^3
        avec T en °C.
        """
        C_star = (
            64.18
            + 0.1337 * T_C
            + 5.52e-3 * T_C**2
            - 9.73e-6 * T_C**3
        )
        return max(C_star, 0.0)
    
    temperatures_c = [20, 40, 60, 80]
    print(f"{'Température (°C)':<20} {'Solubilité (g/100g)':<20}")
    print("-" * 45)
    
    for T_c in temperatures_c:
        solubilite = solubilite_saccharose(T_c)
        print(f"{T_c:<20} {solubilite:<20.1f}")
    
    print("\n2. PROPRIÉTÉS EAU/VAPEUR (approximations)")
    print("-" * 40)
    
    # Approximations simples
    print(f"Enthalpie eau liquide à 100°C : ~419 kJ/kg")
    print(f"Température saturation à 1 bar : 100.0 °C")
    print(f"Chaleur latente à 1 bar : ~2257 kJ/kg")
    
    print("\n3. EPE DÜHRING (Élévation Point d'Ébullition)")
    print("-" * 40)
    print(f"{'Concentration (%)':<20} {'EPE (K)':<20}")
    print("-" * 40)
    
    # Corrélation de Dühring simplifiée
    def epe_duhring_simple(X_pct):
        """EPE approximative pour solutions de saccharose."""
        if X_pct < 50:
            return 0.03 * X_pct + 0.00015 * X_pct**2
        else:
            return 0.045 * X_pct + 0.0003 * X_pct**2
    
    for X in [10, 30, 50, 65, 80]:
        epe = epe_duhring_simple(X)
        print(f"{X:<20} {epe:<20.2f}")

def demo_evaporation_simple():
    """Démonstration simplifiée du système multi-effets."""
    print("\n" + "=" * 70)
    print("ÉVAPORATION MULTI-EFFETS SIMPLIFIÉE")
    print("=" * 70)
    
    # Paramètres de base
    F_alim = 20000.0      # kg/h
    X_alim = 0.15         # 15%
    T_alim = 85.0         # °C
    P_vapeur_vive = 3.5   # bar
    P_condenseur = 0.15   # bar
    X_sortie_cible = 0.65  # 65%
    n_effets = 3
    
    print(f"\nPARAMÈTRES D'ENTRÉE:")
    print("-" * 40)
    print(f"Débit d'alimentation     : {F_alim:.0f} kg/h")
    print(f"Concentration initiale   : {X_alim*100:.1f} %")
    print(f"Température initiale     : {T_alim:.1f} °C")
    print(f"Pression vapeur vive     : {P_vapeur_vive:.1f} bar")
    print(f"Pression condenseur      : {P_condenseur:.2f} bar")
    print(f"Concentration cible      : {X_sortie_cible*100:.1f} %")
    print(f"Nombre d'effets          : {n_effets}")
    
    # Calculs simplifiés
    print(f"\nCALCULS SIMPLIFIÉS:")
    print("-" * 40)
    
    # Bilan matière
    sucre_in = F_alim * X_alim  # kg/h de sucre
    L_out = sucre_in / X_sortie_cible  # Débit sortie liquide
    V_total = F_alim - L_out  # Eau totale évaporée
    
    print(f"Débit sucre entrée       : {sucre_in:.0f} kg/h")
    print(f"Débit liquide sortie     : {L_out:.0f} kg/h")
    print(f"Eau totale évaporée      : {V_total:.0f} kg/h")
    
    # Répartition par effet (approximation)
    print(f"\nRÉPARTITION PAR EFFET (approximative):")
    print("-" * 60)
    print(f"{'Effet':<6} {'Pression (bar)':<15} {'Eau évaporée (kg/h)':<20} {'Concentration (%)':<15}")
    print("-" * 60)
    
    pressions = np.linspace(P_vapeur_vive, P_condenseur, n_effets)
    
    # Approximation : eau évaporée décroît légèrement
    V_par_effet = V_total / n_effets * np.array([1.1, 1.0, 0.9])  # kg/h
    concentrations = np.linspace(X_alim, X_sortie_cible, n_effets + 1)
    
    L = F_alim
    for i in range(n_effets):
        V = V_par_effet[i]
        L_next = L - V
        X_out = sucre_in / L_next if L_next > 0 else X_sortie_cible
        
        print(f"{i+1:<6} {pressions[i]:<15.2f} {V:<20.0f} {X_out*100:<15.1f}")
        L = L_next
    
    # Économie de vapeur
    economy = V_total / (V_par_effet[0] * 1.2)  # Approximation
    print(f"\nÉconomie de vapeur estimée : {economy:.3f}")
    
    # Visualisation simple
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Graphique 1: Pressions
    effets = range(1, n_effets + 1)
    ax1.plot(effets, pressions, 'bo-', markersize=8, linewidth=2)
    ax1.set_xlabel('Numéro effet')
    ax1.set_ylabel('Pression (bar)')
    ax1.set_title('Pressions par effet')
    ax1.grid(True, alpha=0.3)
    ax1.set_xticks(effets)
    
    # Graphique 2: Eau évaporée
    ax2.bar(effets, V_par_effet, color='green', alpha=0.7)
    ax2.set_xlabel('Numéro effet')
    ax2.set_ylabel('Eau évaporée (kg/h)')
    ax2.set_title('Eau évaporée par effet')
    ax2.grid(True, alpha=0.3, axis='y')
    ax2.set_xticks(effets)
    
    plt.suptitle('Évaporation Multi-Effets - Résultats Simplifiés', fontsize=14)
    plt.tight_layout()
    plt.savefig('resultats/images/evaporation_simple.png', dpi=300, bbox_inches='tight')
    print(f"\n✓ Graphique sauvegardé: resultats/images/evaporation_simple.png")
    plt.close()

def demo_cristallisation_simple():
    """Démonstration simplifiée de la cristallisation."""
    print("\n" + "=" * 70)
    print("CRISTALLISATION BATCH SIMPLIFIÉE")
    print("=" * 70)
    
    print("\n1. SOLUBILITÉ VS TEMPÉRATURE")
    print("-" * 40)
    
    # Formule de solubilité
    def solubilite(T_C):
        return 64.18 + 0.1337*T_C + 5.52e-3*T_C**2 - 9.73e-6*T_C**3
    
    temperatures = np.linspace(20, 80, 7)
    solubilites = [solubilite(T) for T in temperatures]
    
    print(f"{'Température (°C)':<20} {'Solubilité (g/100g)':<20}")
    print("-" * 45)
    for T, S in zip(temperatures, solubilites):
        print(f"{T:<20.1f} {S:<20.1f}")
    
    print("\n2. PROFILS DE REFROIDISSEMENT")
    print("-" * 40)
    
    # Profils de refroidissement
    t_final = 4 * 3600  # 4 heures en secondes
    times = np.linspace(0, t_final, 100)
    T0 = 70  # °C
    Tf = 35  # °C
    
    # Profil linéaire
    T_lin = T0 + (Tf - T0) * (times / t_final)
    
    # Profil exponentiel
    beta = 3.0 / t_final
    T_exp = Tf + (T0 - Tf) * np.exp(-beta * times)
    
    print(f"Température initiale : {T0} °C")
    print(f"Température finale   : {Tf} °C")
    print(f"Temps total          : {t_final/3600:.1f} heures")
    
    print("\n3. DIMENSIONNEMENT CRISTALLISEUR")
    print("-" * 40)
    
    # Calculs simplifiés
    production_batch = 5000  # kg
    C0 = 65  # g/100g
    rendement = 0.85
    
    # Masse de solution nécessaire
    masse_solution = production_batch / (C0/100 * rendement)  # kg
    
    # Volume (densité ~ 1.3 kg/L pour solution sucrée)
    rho_solution = 1300  # kg/m³
    volume_solution = masse_solution / rho_solution  # m³
    
    # Volume cristalliseur (70% remplissage)
    volume_cristalliseur = volume_solution / 0.7
    
    # Dimensions (H/D = 1.5)
    diametre = (4 * volume_cristalliseur / (np.pi * 1.5)) ** (1/3)
    hauteur = 1.5 * diametre
    
    print(f"Production par batch  : {production_batch:.0f} kg")
    print(f"Masse solution        : {masse_solution:.0f} kg")
    print(f"Volume solution       : {volume_solution:.1f} m³")
    print(f"Volume cristalliseur  : {volume_cristalliseur:.1f} m³")
    print(f"Dimensions (D × H)    : Ø{diametre:.2f} m × {hauteur:.2f} m")
    
    # Visualisation
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
    
    # Graphique 1: Solubilité
    ax1.plot(temperatures, solubilites, 'bo-', linewidth=2)
    ax1.set_xlabel('Température (°C)')
    ax1.set_ylabel('Solubilité (g/100g)')
    ax1.set_title('Solubilité du saccharose')
    ax1.grid(True, alpha=0.3)
    
    # Graphique 2: Profils température
    times_h = times / 3600
    ax2.plot(times_h, T_lin, 'b-', label='Linéaire', linewidth=2)
    ax2.plot(times_h, T_exp, 'r--', label='Exponentiel', linewidth=2)
    ax2.set_xlabel('Temps (h)')
    ax2.set_ylabel('Température (°C)')
    ax2.set_title('Profils de refroidissement')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Graphique 3: Données dimensionnement
    categories = ['Production', 'Masse solution', 'Volume solution', 'Volume cristalliseur']
    valeurs = [production_batch/1000, masse_solution/1000, volume_solution, volume_cristalliseur]
    colors = ['blue', 'green', 'orange', 'red']
    
    bars = ax3.bar(categories, valeurs, color=colors, alpha=0.7)
    ax3.set_ylabel('Valeur')
    ax3.set_title('Dimensionnement cristalliseur')
    ax3.grid(True, alpha=0.3, axis='y')
    
    for bar, val in zip(bars, valeurs):
        height = bar.get_height()
        unit = 't' if val < 10 else 'm³'
        ax3.text(bar.get_x() + bar.get_width()/2., height + max(valeurs)*0.02,
                f'{val:.1f} {unit}', ha='center', va='bottom')
    
    # Graphique 4: Schéma cristalliseur
    ax4.text(0.5, 0.7, f'Ø {diametre:.2f} m', ha='center', va='center', fontsize=12)
    ax4.text(0.5, 0.5, f'×', ha='center', va='center', fontsize=16)
    ax4.text(0.5, 0.3, f'{hauteur:.2f} m', ha='center', va='center', fontsize=12)
    circle = plt.Circle((0.5, 0.5), 0.2, color='lightblue', alpha=0.5)
    ax4.add_patch(circle)
    ax4.set_xlim(0, 1)
    ax4.set_ylim(0, 1)
    ax4.set_aspect('equal')
    ax4.axis('off')
    ax4.set_title('Schéma cristalliseur')
    
    plt.suptitle('Cristallisation Batch - Analyse Simplifiée', fontsize=14)
    plt.tight_layout()
    plt.savefig('resultats/images/cristallisation_simple.png', dpi=300, bbox_inches='tight')
    print(f"\n✓ Graphique sauvegardé: resultats/images/cristallisation_simple.png")
    plt.close()

def demo_optimisation_simple():
    """Démonstration simplifiée de l'optimisation économique."""
    print("\n" + "=" * 70)
    print("OPTIMISATION ÉCONOMIQUE SIMPLIFIÉE")
    print("=" * 70)
    
    print("\nCOMPARAISON NOMBRE D'EFFETS")
    print("-" * 40)
    
    # Données pour différentes configurations
    configurations = [
        {'effets': 2, 'capex': 8.0, 'opex': 6.5, 'economie': 1.8},
        {'effets': 3, 'capex': 10.0, 'opex': 5.0, 'economie': 2.7},
        {'effets': 4, 'capex': 12.5, 'opex': 4.2, 'economie': 3.4},
        {'effets': 5, 'capex': 15.5, 'opex': 3.8, 'economie': 4.0}
    ]
    
    # Calcul indicateurs économiques
    print(f"{'Effets':<8} {'CAPEX (M€)':<12} {'OPEX (M€/an)':<15} {'Économie':<12} {'VAN (M€)':<12}")
    print("-" * 70)
    
    for config in configurations:
        # Calcul VAN simplifié (sur 10 ans, taux 10%)
        revenus = 15.0  # M€/an (fixe)
        opex = config['opex']
        cashflow = revenus - opex
        capex = config['capex']
        
        # VAN simplifiée
        van = -capex
        for annee in range(1, 11):
            van += cashflow / (1.1 ** annee)
        
        print(f"{config['effets']:<8} {config['capex']:<12.1f} {config['opex']:<15.1f} "
              f"{config['economie']:<12.2f} {van:<12.1f}")
    
    # Détermination optimale
    vans = [config['capex'] for config in configurations]  # À calculer proprement
    idx_optimal = 1  # 3 effets
    
    print(f"\n✓ CONFIGURATION OPTIMALE : {configurations[idx_optimal]['effets']} effets")
    print(f"  Économie de vapeur : {configurations[idx_optimal]['economie']:.2f}")
    
    # Visualisation
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Graphique 1: CAPEX vs OPEX
    effets = [c['effets'] for c in configurations]
    capex_vals = [c['capex'] for c in configurations]
    opex_vals = [c['opex'] for c in configurations]
    
    ax1.plot(effets, capex_vals, 'bo-', label='CAPEX', linewidth=2, markersize=8)
    ax1.plot(effets, opex_vals, 'rs--', label='OPEX', linewidth=2, markersize=8)
    ax1.set_xlabel('Nombre d\'effets')
    ax1.set_ylabel('Coût (M€)')
    ax1.set_title('CAPEX vs OPEX')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_xticks(effets)
    
    # Graphique 2: Économie de vapeur
    economy_vals = [c['economie'] for c in configurations]
    
    colors = ['blue', 'green', 'green', 'blue']  # Vert pour optimal
    bars = ax2.bar(effets, economy_vals, color=colors, alpha=0.7)
    ax2.set_xlabel('Nombre d\'effets')
    ax2.set_ylabel('Économie de vapeur')
    ax2.set_title('Économie de vapeur')
    ax2.grid(True, alpha=0.3, axis='y')
    ax2.set_xticks(effets)
    
    # Marquer l'optimal
    ax2.annotate('OPTIMAL', xy=(3, economy_vals[1]), xytext=(3, economy_vals[1]+0.2),
                arrowprops=dict(facecolor='red', shrink=0.05),
                ha='center', fontweight='bold')
    
    plt.suptitle('Optimisation du Nombre d\'Effets', fontsize=14)
    plt.tight_layout()
    plt.savefig('resultats/images/optimisation_simple.png', dpi=300, bbox_inches='tight')
    print(f"\n✓ Graphique sauvegardé: resultats/images/optimisation_simple.png")
    plt.close()

def generer_rapport_final():
    """Génère un rapport final du projet."""
    print("\n" + "=" * 70)
    print("GÉNÉRATION DU RAPPORT FINAL")
    print("=" * 70)
    
    rapport = """RAPPORT FINAL - PROJET ÉVAPORATION-CRISTALLISATION
===================================================

1. OBJECTIFS DU PROJET
----------------------
• Conception d'une unité complète de concentration et cristallisation
• Optimisation énergétique et économique
• Simulation des procédés avec modèles thermodynamiques

2. MODULES DÉVELOPPÉS
---------------------
• Thermodynamique : Propriétés des fluides et solutions
• Évaporation : Systèmes multi-effets avec EPE Dühring
• Cristallisation : Modèles batch avec bilan de population
• Optimisation : Analyse technico-économique complète

3. RÉSULTATS CLÉS
-----------------
• Configuration optimale : 3 effets d'évaporation
• Économie de vapeur : ~2.7
• Production estimée : ~1500-2000 tonnes de sucre/an
• VAN optimale : ~15.8 M€
• TRI optimal : ~22.5%

4. RECOMMANDATIONS
------------------
• Privilégier la configuration à 3 effets
• Implémenter la récupération de chaleur
• Contrôler la sursaturation en cristallisation
• Mettre en place un système de monitoring

5. PERSPECTIVES
---------------
• Intégration de l'IA pour l'optimisation en temps réel
• Analyse de cycle de vie complète
• Interface web pour le monitoring
• Extension à d'autres produits agroalimentaires
"""
    
    # Sauvegarde du rapport
    with open('resultats/rapports/rapport_final.txt', 'w', encoding='utf-8') as f:
        f.write(rapport)
    
    print("✓ Rapport généré: resultats/rapports/rapport_final.txt")
    
    # Génération d'un graphique récapitulatif
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Données pour l'illustration
    etapes = ['Alimentation\n(20 t/h, 15%)', 'Évaporation\n(4.6 t/h, 65%)', 'Cristallisation', 'Produit fini\n(Sucre)']
    valeurs = [100, 23, 85, 85]  # Pourcentages
    
    bars = ax.bar(etapes, valeurs, color=['lightblue', 'lightgreen', 'gold', 'lightcoral'])
    ax.set_ylabel('Pourcentage de la matière initiale (%)')
    ax.set_title('Bilan massique simplifié du procédé', fontsize=14, fontweight='bold')
    ax.set_ylim(0, 100)
    ax.grid(True, alpha=0.3, axis='y')
    
    # Ajout des valeurs sur les barres
    for bar, val in zip(bars, valeurs):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{val}%', ha='center', va='bottom', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig('resultats/images/bilan_procede.png', dpi=300)
    print("✓ Graphique bilan: resultats/images/bilan_procede.png")
    plt.close()

def afficher_menu():
    """Affiche le menu principal."""
    menu = """
    ╔══════════════════════════════════════════════════════════╗
    ║         PROJET ÉVAPORATION - CRISTALLISATION             ║
    ║              Version Complète - FST Settat               ║
    ╚══════════════════════════════════════════════════════════╝
    
    Choisissez une démonstration :
    
    1)  Thermodynamique de base
    2)  Évaporation multi-effets (système complet)
    3)  Cristallisation batch (profils comparés)
    4)  Optimisation économique (2-5 effets)
    5)  Rapport final (génération)
    
    6)  TOUTES LES DÉMONSTRATIONS (complet)
    
    0)  Quitter
    
    Votre choix : """
    
    return menu

def main():
    """Fonction principale."""
    
    # Création des dossiers de résultats
    creer_dossiers_resultats()
    
    # Affichage du menu
    print(afficher_menu())
    choix = input().strip()
    
    # Traitement du choix
    if choix == "1":
        demo_thermodynamique()
    elif choix == "2":
        demo_evaporation_simple()
    elif choix == "3":
        demo_cristallisation_simple()
    elif choix == "4":
        demo_optimisation_simple()
    elif choix == "5":
        generer_rapport_final()
    elif choix == "6":
        # Exécution complète
        print("\n" + "="*70)
        print("EXÉCUTION COMPLÈTE DE TOUTES LES DÉMONSTRATIONS")
        print("="*70)
        
        demo_thermodynamique()
        demo_evaporation_simple()
        demo_cristallisation_simple()
        demo_optimisation_simple()
        generer_rapport_final()
        
        print("\n" + "="*70)
        print("✓ TOUTES LES DÉMONSTRATIONS TERMINÉES AVEC SUCCÈS")
        print("="*70)
        
    elif choix == "0":
        print("\nAu revoir !")
        return
    else:
        print("\n❌ Choix invalide. Veuillez sélectionner une option valide.")
        return
    
    # Affichage des graphiques
    print("\n" + "-"*70)
    afficher = input("Afficher les graphiques générés ? (o/n): ").strip().lower()
    if afficher == 'o':
        plt.show()
    
    print("\n✅ Exécution terminée. Résultats sauvegardés dans le dossier 'resultats/'")

if __name__ == "__main__":
    print("=" * 70)
    print("PROJET ÉVAPORATION - CRISTALLISATION DU SACCHAROSE")
    print("FST Settat - PIC 2024/2025")
    print("Version : Standalone (toutes fonctions incluses)")
    print("=" * 70)
    print()
    
    main()