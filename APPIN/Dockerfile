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