# Projet Évaporation - Cristallisation du Saccharose

**Conception et Simulation d'une Unité Intégrée d'Évaporation à Multiples Effets et de Cristallisation**

---

## 📋 Informations du Projet

- **Filière**: Procédés et Ingénierie Chimique (PIC)
- **Niveau**: 3ème Année - Cycle d'Ingénieur
- **Université**: Hassan 1 - FST Settat
- **Année Universitaire**: 2025-2026

---

## 🎯 Objectifs

Ce projet vise à :

1. Dimensionner un système d'évaporation à multiples effets
2. Modéliser la cinétique de cristallisation batch
3. Résoudre le bilan de population des cristaux
4. Optimiser le procédé sous contraintes technico-économiques
5. Réaliser une analyse de sensibilité paramétrique

---

## 📁 Structure du Projet

```
projet_evaporation_cristallisation/
│
├── thermodynamique.py          # Propriétés physiques (EPE Dühring, CoolProp)
├── evaporateurs.py             # Modèles évaporation multi-effets
├── cristallisation.py          # Cinétique + bilan de population
├── optimisation.py             # Analyse technico-économique (VAN, TRI)
├── main.py                     # Point d'entrée - démos
├── app_flask.py                # Interface web Flask
├                               # Tests unitaires pytest
├── requirements.txt            # Dépendances Python
├── README.md                   # Ce fichier
│
├── resultats/                  # Résultats de simulation
│   ├── graphiques/
│  
│
└── rapport/                    # Rapport LaTeX
    ├── rapport.tex
   
```

---

## 🚀 Installation

### 1. Créer un environnement virtuel (recommandé)

```bash
python -m venv venv

# Windows
venv\Scripts\activate

# Linux/Mac
source venv/bin/activate
```

### 2. Installer les dépendances

```bash
pip install -r requirements.txt
```

### 3. Vérifier l'installation

```bash
python -c "from CoolProp.CoolProp import PropsSI; print('CoolProp OK')"
pytest test_projet.py -v
```

---

## 💻 Utilisation

### Mode Console Interactif

```bash
python main.py
```

**Menu disponible** :
1. Test thermodynamique (EPE Dühring corrigé)
2. Évaporation multi-effets (cas de base du sujet)
3. Cristallisation batch avec bilan de population
4. Comparaison des 3 profils de refroidissement
5. Optimisation nombre d'effets (2-5) avec VAN/TRI
6. Procédé complet (évaporation + cristallisation)
7. **Tout exécuter** (démo complète)

### Interface Web Flask

```bash
python app_flask.py
```

Puis ouvrir : http://127.0.0.1:5000/

---

## 📊 Fonctionnalités Implémentées

### ✅ Partie 1 : Évaporation 

- [x] **Modélisation thermodynamique** 
  - Corrélation de Dühring correcte (coefficients selon concentration)
  - Bilans matière et énergie complets
  - Coefficient d'encrassement Rf = 0.0002 m²·K/W
  - Validation bilans (erreur < 1%)

- [x] **Optimisation énergétique**
  - Économie de vapeur calculée
  - Étude 2, 3, 4, 5 effets
  - Configuration optimale (VAN maximale)

- [x] **Analyse de sensibilité** 
  - Pression vapeur de chauffe
  - Concentration finale
  - Débit alimentation
  - Température alimentation

### ✅ Partie 2 : Cristallisation

- [x] **Modélisation cinétique** 
  - Solubilité (équation 1 de l'énoncé)
  - Nucléation (équation 3)
  - Croissance (équation 4)
  - **Bilan de population** résolu (méthode des classes)
  - Distribution de taille n(L,t)
  - Moments m0, m1, m3

- [x] **Stratégie de refroidissement** 
  - Profil linéaire
  - Profil exponentiel
  - Profil optimal (S constant)
  - Comparaison DTG, L50, CV

- [x] **Dimensionnement cristalliseur** 
  - Volume requis
  - Puissance agitation
  - Surface serpentin
  - Temps de résidence

### ✅ Partie 3 : Intégration 

- [x] **Intégration énergétique** 
  - Récupération chaleur condensats
  - Préchauffage alimentation
  - Économie énergétique calculée

- [x] **Analyse technico-économique** 
  - CAPEX détaillé (évaporateurs, cristalliseur, échangeurs)
  - OPEX complet (vapeur, électricité, eau, main d'œuvre, maintenance)
  - **VAN** (Valeur Actuelle Nette)
  - **TRI** (Taux de Rendement Interne)
  - **ROI** (temps de retour)
  - Coût de production par tonne

---

## 🧪 Tests Unitaires

Exécuter les tests avec **pytest** :

```bash
# Tous les tests
pytest test_projet.py -v

# Tests spécifiques
pytest test_projet.py::TestThermodynamique -v

# Avec couverture
pytest test_projet.py --cov=. --cov-report=html
```

**Couverture actuelle** : > 80%

---

## 📈 Résultats Clés

### Cas de Base (3 effets)

| Paramètre | Valeur |
|-----------|--------|
| Alimentation | 20 000 kg/h à 15% |
| Sortie évaporation | 4 615 kg/h à 65% |
| Vapeur vive | ~5 750 kg/h |
| Économie de vapeur | ~2.67 |
| Surface totale | ~180 m² |

### Optimisation Économique

| Nombre d'effets | VAN (MMAD) | TRI (%) | ROI (ans) |
|-----------------|------------|---------|-----------|
| 2 | 12.5 | 18.2 | 4.8 |
| **3** | **15.8** | **22.5** | **3.9** ✓ |
| 4 | 14.2 | 20.1 | 4.3 |
| 5 | 12.9 | 17.8 | 4.9 |

**Recommandation** : **3 effets** (VAN maximale)

---

## 📚 Références Techniques

### Équations Clés

#### EPE (Dühring)
```
Si X < 50% : EPE = 0.03·X + 0.00015·X²
Si X ≥ 50% : EPE = 0.045·X + 0.0003·X²
```

#### Solubilité Saccharose
```
C* = 64.18 + 0.1337·T + 5.52×10⁻³·T² - 9.73×10⁻⁶·T³
```

#### Croissance Cristaux
```
G = kg·S^g·exp(-Eg/RT)
avec kg = 2.8×10⁻⁷ m/s, g = 1.5, Eg = 45 kJ/mol
```

### Documentation Externe

- [CoolProp Documentation](http://www.coolprop.org)
- [Perry's Chemical Engineers' Handbook](https://www.accessengineeringlibrary.com)
- Mullin, J.W. "Crystallization" (4th ed.)

---

## ⚠️ Corrections Importantes

### Par rapport au code initial

1. **EPE Dühring** : Coefficients corrigés selon l'énoncé
2. **Coefficient d'encrassement** : Rf = 0.0002 appliqué
3. **Bilan de population** : Résolution complète de l'EDP
4. **3 profils refroidissement** : Tous implémentés
5. **VAN/TRI** : Calculs complets ajoutés
6. **Intégration énergétique** : Récupération chaleur condensats
7. **Tests unitaires** : Validation complète

---

## 👥 Contribution

### Répartition des Tâches (à adapter)

- **Étudiant 1** : Évaporation + Thermodynamique
- **Étudiant 2** : Cristallisation + Optimisation

### Travail en Équipe

- Git/GitHub pour versioning
- Réunions hebdomadaires
- Documentation continue

---

## 📝 Livrables

### Code Python (40%)
- [x] Modules complets et documentés
- [x] Tests unitaires (pytest)
- [x] Gestion d'erreurs
- [x] Docstrings style NumPy

### Rapport Technique (40%)
- [ ] LaTeX (template fourni)
- [ ] 10 pages max
- [ ] Figures haute qualité
- [ ] Bibliographie

### Présentation Orale (20%)
- [ ] PowerPoint/Beamer
- [ ] 10 minutes
- [ ] Démo code

---

## 📅 Planning

- **Semaine 1** : Évaporation + validation
- **Semaine 2 (J1-3)** : Cristallisation + bilan population
- **Semaine 2 (J4-5)** : Intégration + optimisation
- **Semaine 3 (J1-3)** : Rapport + présentation

**Date limite** : 15/12/2025

---

## 🐛 Problèmes Connus

- Bilan de population : calcul intensif (optimiser avec Numba)
- Profil optimal : résolution inverse simplifiée

---

## 📧 Contact

Pour questions techniques :
- Email : amiri.ikr.fst@uhp.ac.ma / smajouk@gmail.com
- Responsable module : Prof.BAKHER

---

## 📄 Licence

Ce projet est réalisé dans le cadre académique de la FST Settat.

---
