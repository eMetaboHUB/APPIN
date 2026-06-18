# APPIN - Guide d'installation Docker

Ce guide explique comment installer et utiliser APPIN via Docker.

## Prérequis

- **Docker Desktop** installé sur votre machine
  - [Télécharger Docker Desktop pour Windows](https://www.docker.com/products/docker-desktop/)
  - [Télécharger Docker Desktop pour Mac](https://www.docker.com/products/docker-desktop/)
  - Pour Linux : [Instructions d'installation](https://docs.docker.com/engine/install/)

## Installation

### 1. Télécharger APPIN

Téléchargez et extrayez l'archive APPIN dans un dossier de votre choix, par exemple :
- Windows : `C:\Users\VotreNom\Documents\APPIN`
- Mac/Linux : `~/Documents/APPIN`

### 2. Construire l'image Docker

Ouvrez un terminal (CMD, PowerShell, ou Terminal) et naviguez vers le dossier APPIN :

```bash
cd chemin/vers/APPIN
```

Construisez l'image Docker :

```bash
docker build -t appin .
```

> ⏱️ Cette étape peut prendre plusieurs minutes lors de la première exécution.

## Utilisation

### Lancer APPIN

Pour lancer l'application avec accès à vos données spectrales (Vous devez bien mettre le chemin qui mène jusqu'à vos données) :

```bash
docker run -p 80:3838 -v /chemin/vers/vos/donnees:/data appin
```

**Exemples :**

Windows (CMD/PowerShell) :
```cmd
docker run -p 80:3838 -v C:\Users\VotreNom\Documents\Spectres:/data appin
```

Mac/Linux :
```bash
docker run -p 80:3838 -v ~/Documents/Spectres:/data appin
```

### Accéder à l'application

Ouvrez votre navigateur et allez sur : **http://localhost**

Vos fichiers spectraux seront accessibles dans le dossier `/data` de l'application.

## Options avancées

### Allouer plus de ressources

Pour des analyses intensives, vous pouvez allouer plus de mémoire et de CPU :

```bash
docker run -p 80:3838 -v /chemin/vers/vos/donnees:/data --memory=8g --cpus=4 appin
```

- `--memory=8g` : Alloue 8 Go de RAM au container
- `--cpus=4` : Alloue 4 cœurs CPU

### Lancer en arrière-plan

Pour lancer APPIN en arrière-plan :

```bash
docker run -d -p 80:3838 -v /chemin/vers/vos/donnees:/data --name appin appin
```

Voir les logs :
```bash
docker logs appin
```

Arrêter l'application :
```bash
docker stop appin
```

Relancer l'application :
```bash
docker start appin
```

### Utiliser un port différent

Si le port 80 est déjà utilisé, vous pouvez en choisir un autre :

```bash
docker run -p 8080:3838 -v /chemin/vers/vos/donnees:/data appin
```

Accédez ensuite à : **http://localhost:8080**

## Configuration WSL2 (Windows uniquement)

Si vous utilisez Docker avec WSL2 et souhaitez augmenter les ressources disponibles :

1. Créez/modifiez le fichier `%USERPROFILE%\.wslconfig` :

Pour ouvrir le fichier vous devez ecrire dans le cmd : notepad %USERPROFILE%\.wslconfig

```ini
[wsl2]
memory=16GB
processors=4
swap=4GB
```

2. Redémarrez WSL :
```cmd
wsl --shutdown
```

3. Relancez Docker Desktop

## Dépannage

### Le port est déjà utilisé

```
Error: port is already allocated
```

Arrêtez le container existant :
```bash
docker ps
docker stop CONTAINER_ID
```

Ou utilisez un autre port (voir section ci-dessus).

### L'application ne se lance pas

Vérifiez les logs :
```bash
docker logs CONTAINER_ID
```

### Reconstruire l'image après une mise à jour

```bash
docker build --no-cache -t appin .
```

## Commandes utiles

| Commande | Description |
|----------|-------------|
| `docker ps` | Lister les containers en cours |
| `docker ps -a` | Lister tous les containers |
| `docker stop CONTAINER_ID` | Arrêter un container |
| `docker rm CONTAINER_ID` | Supprimer un container |
| `docker images` | Lister les images |
| `docker stats` | Voir l'utilisation des ressources en temps réel |
| `docker logs -f CONTAINER_ID` | Suivre les logs en temps réel |
