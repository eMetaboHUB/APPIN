# APPIN - Guide Docker avec shiny2docker

## Etape 1 : Generer le Dockerfile avec shiny2docker

```r
# Installer shiny2docker
install.packages("shiny2docker")
library(shiny2docker)

# Se placer dans le dossier de l'application
setwd("C:/Users/juguibert/Documents/APPIN_1.0.0")

# Generer le Dockerfile automatiquement
shiny2docker(path = ".")
```

Cela cree deux fichiers :
- `Dockerfile`
- `.dockerignore`

---

## Etape 2 : Adapter le Dockerfile pour Python/TensorFlow

Le Dockerfile genere par shiny2docker ne gere pas Python/TensorFlow.
Remplacer le contenu du Dockerfile par la version ci-dessous qui inclut :
- Installation de Python3 et creation d'un environnement virtuel
- Installation de TensorFlow 2.16.1
- Configuration de reticulate pour utiliser le bon Python

### Dockerfile complet pour APPIN :

```dockerfile
FROM rocker/geospatial:4.2.2

RUN apt-get update -y && apt-get install -y \
    make pandoc libx11-dev libcurl4-openssl-dev libssl-dev \
    zlib1g-dev libglpk-dev libxml2-dev libfftw3-dev libicu-dev \
    libjpeg-dev libpng-dev libtiff-dev python3 python3-pip python3-venv git \
    && rm -rf /var/lib/apt/lists/*

# Create virtualenv and install TensorFlow
RUN python3 -m venv /opt/venv
RUN /opt/venv/bin/pip install --upgrade pip
RUN /opt/venv/bin/pip install --no-cache-dir tensorflow==2.16.1 tf_keras numpy

# Verify TensorFlow installation
RUN /opt/venv/bin/python -c "import tensorflow as tf; print('TF version:', tf.__version__)"

# Create symlink so reticulate can find python3.10
RUN ln -sf /opt/venv/bin/python /usr/bin/python3.10 || true

RUN mkdir -p /usr/local/lib/R/etc/ /usr/lib/R/etc/

# Force RETICULATE_PYTHON to use venv
RUN echo 'Sys.setenv(RETICULATE_PYTHON = "/opt/venv/bin/python")' | tee /usr/local/lib/R/etc/Rprofile.site | tee /usr/lib/R/etc/Rprofile.site
RUN echo 'options(renv.config.pak.enabled = FALSE, repos = c(CRAN = "https://cran.rstudio.com/"), download.file.method = "libcurl", Ncpus = 4)' | tee -a /usr/local/lib/R/etc/Rprofile.site | tee -a /usr/lib/R/etc/Rprofile.site

RUN echo 'RETICULATE_PYTHON=/opt/venv/bin/python' | tee /usr/local/lib/R/etc/Renviron.site | tee /usr/lib/R/etc/Renviron.site

ENV RETICULATE_PYTHON=/opt/venv/bin/python
ENV PATH="/opt/venv/bin:$PATH"

RUN R -e 'install.packages("remotes")'
RUN R -e 'remotes::install_version("renv", version = "1.0.3")'

COPY renv.lock renv.lock
RUN --mount=type=cache,id=renv-cache,target=/root/.cache/R/renv R -e 'renv::restore()'

# Clear any reticulate Python config cache
RUN rm -rf /root/.local/share/r-reticulate || true

# Verify reticulate sees TensorFlow
RUN R -e 'library(reticulate); py_config(); tf <- import("tensorflow"); print(paste("R sees TensorFlow:", tf$"__version__"))'

WORKDIR /srv/shiny-server/
COPY . /srv/shiny-server/
RUN mv /srv/shiny-server/Shine.R /srv/shiny-server/app.R

EXPOSE 3838
CMD R -e 'shiny::runApp("/srv/shiny-server",host="0.0.0.0",port=3838)'
```

---

## Etape 3 : Construire l'image Docker

Dans un terminal (PowerShell ou cmd), se placer dans le dossier et lancer :

```bash
# Se placer dans le dossier de l'application
cd C:\Users\juguibert\Documents\APPIN_1.0.0

# Construire l'image
docker build -t appin .
```

Duree : 15-45 minutes selon la connexion internet.

---

## Etape 4 : Lancer le conteneur

```bash
docker run -d -p 3838:3838 --name appin_container appin
```

Ouvrir dans le navigateur : http://localhost:3838

---

## Commandes utiles

```bash
# Voir les conteneurs en cours
docker ps

# Voir les logs
docker logs appin_container

# Suivre les logs en temps reel
docker logs -f appin_container

# Arreter
docker stop appin_container

# Redemarrer
docker start appin_container

# Supprimer le conteneur
docker rm appin_container

# Supprimer l'image
docker rmi appin
```

---

## Explication des ajouts pour Python/TensorFlow

### 1. Installation de Python dans apt-get

```dockerfile
python3 python3-pip python3-venv
```

### 2. Creation de l'environnement virtuel

```dockerfile
RUN python3 -m venv /opt/venv
RUN /opt/venv/bin/pip install --upgrade pip
RUN /opt/venv/bin/pip install --no-cache-dir tensorflow==2.16.1 tf_keras numpy
```

### 3. Configuration de reticulate

Plusieurs methodes pour s'assurer que R/reticulate trouve le bon Python :

```dockerfile
# Variable d'environnement
ENV RETICULATE_PYTHON=/opt/venv/bin/python

# Dans Rprofile.site (au demarrage de R)
RUN echo 'Sys.setenv(RETICULATE_PYTHON = "/opt/venv/bin/python")' | tee /usr/local/lib/R/etc/Rprofile.site

# Dans Renviron.site
RUN echo 'RETICULATE_PYTHON=/opt/venv/bin/python' | tee /usr/local/lib/R/etc/Renviron.site
```

### 4. Verification

```dockerfile
# Verifier que TensorFlow est installe
RUN /opt/venv/bin/python -c "import tensorflow as tf; print('TF version:', tf.__version__)"

# Verifier que R voit TensorFlow
RUN R -e 'library(reticulate); tf <- import("tensorflow"); print(tf$"__version__")'
```

### 5. Renommage du fichier principal

```dockerfile
RUN mv /srv/shiny-server/Shine.R /srv/shiny-server/app.R
```

---

## Resume

| Etape | Commande |
|-------|----------|
| 1. Generer Dockerfile de base | `shiny2docker(path = ".")` |
| 2. Remplacer par version avec Python | Copier le Dockerfile ci-dessus |
| 3. Construire l'image | `docker build -t appin .` |
| 4. Lancer le conteneur | `docker run -d -p 3838:3838 appin` |
| 5. Acceder a l'application | http://localhost:3838 |
