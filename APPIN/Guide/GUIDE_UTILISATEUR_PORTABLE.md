# APPIN - Guide Utilisateur (Version Portable)

## 📦 Contenu du dossier

```
APPIN_Portable/
├── Launch_APPIN.bat    ← Double-cliquez ici pour lancer
├── R-Portable/         ← R (ne pas modifier)
├── Python/             ← Python + TensorFlow (ne pas modifier)
├── app/                ← Application (ne pas modifier)
└── data/               ← Vos données (optionnel)
```

---

## 🚀 Démarrage rapide

### 1. Lancer l'application

**Double-cliquez sur `Launch_APPIN.bat`**

Une fenêtre noire apparaît — **ne la fermez pas**, c'est normal !

Votre navigateur s'ouvrira automatiquement sur l'application.

> Si le navigateur ne s'ouvre pas automatiquement, ouvrez-le manuellement et allez à : **http://localhost:3838**

### 2. Charger vos données

1. Dans le panneau de gauche, cliquez sur **"📂 1. Load Data"**
2. Cliquez sur **"Browse..."** ou **"Select Folder"**
3. Naviguez vers votre dossier contenant les spectres Bruker
4. Sélectionnez le dossier et validez

### 3. Configurer et visualiser

1. **Sélectionnez le type de spectre** (très important !) :
   - TOCSY
   - HSQC
   - COSY
   - UFCOSY

2. Cliquez sur **"Generate Plot"** pour afficher le spectre

3. Ajustez le seuil (Threshold) si nécessaire

### 4. Détecter les pics

Deux méthodes disponibles :

| Méthode | Description | Vitesse |
|---------|-------------|---------|
| **Local Max** | Détection par maxima locaux + clustering DBSCAN | Rapide |
| **CNN** | Détection par réseau de neurones | Plus lent mais plus précis |

### 5. Fermer l'application

- Fermez simplement la **fenêtre noire**
- Ou appuyez sur **Ctrl+C** dans la fenêtre noire

---

## ⚠️ Résolution de problèmes

### L'application ne se lance pas

- Vérifiez que vous avez double-cliqué sur `Launch_APPIN.bat`
- Essayez : **clic droit** → **Exécuter en tant qu'administrateur**
- Vérifiez qu'aucun antivirus ne bloque l'exécution

### Le navigateur ne s'ouvre pas

1. Ouvrez manuellement votre navigateur (Chrome, Firefox, Edge...)
2. Allez à l'adresse : `http://localhost:3838`

### Erreur "Port 3838 déjà utilisé"

Une autre instance d'APPIN tourne peut-être déjà :

1. Fermez toutes les fenêtres noires APPIN
2. Ou redémarrez votre ordinateur
3. Réessayez

### Erreur TensorFlow / CNN

- Le CNN nécessite beaucoup de mémoire RAM
- Fermez d'autres applications (Chrome, Word, etc.)
- Utilisez **Local Max** comme alternative

### L'application est très lente

- Les spectres TOCSY volumineux sont longs à charger
- Un spinner s'affiche pendant le traitement — patientez
- La première exécution est plus lente (chargement des modèles)

### Erreur "R not found" ou "Python not found"

Les dossiers `R-Portable` ou `Python` sont manquants ou corrompus :

1. Re-téléchargez APPIN_Portable
2. Extrayez à nouveau l'archive
3. Ne déplacez pas les sous-dossiers individuellement

---

## 💡 Conseils d'utilisation

### Performance

- **Fermez** les applications gourmandes en mémoire avant de lancer APPIN
- Les spectres **HSQC** sont généralement plus rapides à traiter que les **TOCSY**
- Utilisez **Local Max** pour un premier aperçu rapide, puis **CNN** pour affiner

### Organisation des données

Placez vos dossiers Bruker dans le dossier `data/` pour y accéder facilement :

```
APPIN_Portable/
└── data/
    ├── Experiment_001/
    │   └── pdata/1/
    ├── Experiment_002/
    │   └── pdata/1/
    └── ...
```

### Sauvegarde des résultats

Les résultats exportés (CSV, etc.) sont sauvegardés dans le dossier que vous choisissez lors de l'export.

---

## 📧 Support

Pour toute question ou problème :

**Julien Guibert**  
INRAe Toxalim / MetaboHUB

---

## 📋 Informations techniques

| Composant | Version |
|-----------|---------|
| APPIN | 1.0.0 |
| R | 4.2.2 |
| Python | 3.10.11 |
| TensorFlow | 2.16.1 |

**Types de spectres supportés** : TOCSY, HSQC, COSY, UFCOSY

**Format d'entrée** : Dossiers Bruker (avec sous-dossier pdata/1/)
