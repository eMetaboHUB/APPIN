# APPIN - Guide de création du bundle portable

Ce guide explique comment créer la version portable d'APPIN étape par étape.

## Prérequis

Sur ta machine de développement :
- Windows 10/11
- R 4.2.2 installé
- Python 3.10+ installé
- APPIN fonctionnel
- ~10 Go d'espace disque libre

---

## Méthode 1 : Script automatique (recommandé)

### Étape 1 : Configurer le script

1. Ouvre `build_portable.ps1` dans un éditeur de texte
2. Modifie les chemins si nécessaire :

```powershell
$SOURCE_APP_DIR = "C:\Users\juguibert\Documents\APPIN_1.0.0"  # Ton dossier APPIN
$OUTPUT_DIR = "C:\Users\juguibert\Documents\APPIN_Portable"   # Dossier de sortie
```

### Étape 2 : Exécuter le script

1. Clic droit sur `build_portable.ps1`
2. Sélectionne **"Exécuter avec PowerShell"**
3. Si demandé, autorise l'exécution (tape `Y` puis Entrée)
4. Attends la fin du processus (~15-30 minutes)

### Étape 3 : Tester

1. Va dans le dossier `APPIN_Portable` créé
2. Double-clique sur `Launch_APPIN.bat`
3. Vérifie que tout fonctionne

### Étape 4 : Distribuer

1. Compresse le dossier `APPIN_Portable` en `.zip`
2. Partage le fichier zip avec tes utilisateurs

---

## Méthode 2 : Création manuelle (Utile surtout pour R normalement)

Si le script automatique ne fonctionne pas, voici les étapes manuelles :

### Étape 1 : Créer la structure de dossiers

```
APPIN_Portable/
├── R-Portable/
├── Python/
├── app/
└── data/
```

### Étape 2 : Installer R Portable

**Option A** - Copier R existant :
```
Copier tout le contenu de C:\Program Files\R\R-4.2.2\ vers APPIN_Portable\R-Portable\
```

**Option B** - Télécharger R Portable :
1. Va sur https://sourceforge.net/projects/rportable/
2. Télécharge la version 4.2.2
3. Extrait dans `APPIN_Portable\R-Portable\`

### Étape 3 : Installer Python Embedded

1. Télécharge Python embedded :
   https://www.python.org/ftp/python/3.10.11/python-3.10.11-embed-amd64.zip

2. Extrait dans `APPIN_Portable\Python\`

3. Active pip - édite le fichier `python310._pth` :
   - Trouve la ligne `#import site`
   - Enlève le `#` pour avoir `import site`

4. Installe pip :
   ```cmd
   cd APPIN_Portable\Python
   curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
   python.exe get-pip.py
   del get-pip.py
   ```

5. Installe TensorFlow :
   ```cmd
   python.exe -m pip install tensorflow==2.16.1 tf_keras numpy
   ```

6. Si il manque rpy Tools aller le chercher :

C:\Users\juguibert\Documents\APPIN_Portable\R-Portable\library\reticulate\python\rpytools


### Étape 4 : Copier l'application

Copie ces fichiers/dossiers de ton APPIN vers `APPIN_Portable\app\` :

- `Shine.R` → renommer en `app.R`
- `renv.lock`
- `Function/` (dossier entier)
- `R/` (dossier entier)
- `www/` (dossier entier)
- `models/` ou `Model/` (dossier des modèles CNN)

### Étape 5 : Installer les packages R

Ouvre une invite de commandes et exécute :

```cmd
cd APPIN_Portable\app
..\R-Portable\bin\Rscript.exe -e "install.packages('renv', repos='https://cran.rstudio.com/')"
..\R-Portable\bin\Rscript.exe -e "renv::restore(prompt=FALSE)"
```

⚠️ Cette étape peut prendre 20-30 minutes !

### Étape 6 : Ajouter le launcher

Copie `Launch_APPIN.bat` à la racine de `APPIN_Portable\`

### Étape 7 : Tester

Double-clique sur `Launch_APPIN.bat` et vérifie que tout fonctionne.

---

## Résolution de problèmes

### "renv::restore() échoue"

Certains packages peuvent ne pas s'installer. Solutions :

1. Installe les packages problématiques manuellement :
   ```r
   install.packages("nom_du_package")
   ```

2. Ou utilise la bibliothèque de packages de ton R principal :
   ```cmd
   xcopy "C:\Users\juguibert\AppData\Local\R\win-library\4.2" "APPIN_Portable\R-Portable\library" /E /I
   ```

### "TensorFlow ne fonctionne pas"

Vérifie l'installation Python :
```cmd
APPIN_Portable\Python\python.exe -c "import tensorflow; print(tensorflow.__version__)"
```

Si ça échoue, réinstalle TensorFlow :
```cmd
APPIN_Portable\Python\python.exe -m pip uninstall tensorflow
APPIN_Portable\Python\python.exe -m pip install tensorflow==2.16.1
```

### "Le bundle est trop gros"

Le bundle complet fait ~4-6 Go. Pour réduire :

1. Supprime les fichiers de cache :
   ```
   APPIN_Portable\R-Portable\library\*\doc\
   APPIN_Portable\R-Portable\library\*\help\
   APPIN_Portable\Python\Lib\site-packages\tensorflow\include\
   ```

2. Utilise un outil de compression efficace (7-Zip avec compression Ultra)

---

## Checklist finale avant distribution

- [ ] `Launch_APPIN.bat` lance l'application
- [ ] Le navigateur s'ouvre automatiquement
- [ ] Les données peuvent être chargées
- [ ] "Generate Plot" fonctionne
- [ ] "Local Max" détecte des pics
- [ ] "CNN" fonctionne (si activé)
- [ ] L'export fonctionne
- [ ] Le guide utilisateur est inclus

---

## Taille attendue

| Composant | Taille approximative |
|-----------|---------------------|
| R-Portable | ~500 Mo |
| R packages (library) | ~1.5 Go |
| Python embedded | ~100 Mo |
| TensorFlow + deps | ~2 Go |
| Application | ~100 Mo |
| **Total** | **~4-5 Go** |

Compressé en .zip : ~1.5-2 Go
