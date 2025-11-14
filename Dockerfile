# Use a base R image from the Rocker project
FROM rocker/tidyverse:4.5.1

WORKDIR /pipeline
# Install any required system dependencies (example for a common R package dependency)
# You may need to add more `RUN apt-get install -y` commands for specific packages
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libc6 \
    # Add other system dependencies as needed
    && rm -rf /var/lib/apt/lists/*

RUN pip install --break-system-packages crowelab_pyir

RUN Rscript -e "\
    if (!require('BiocManager', quietly = TRUE)) {install.packages('BiocManager')}; \
    BiocManager::install('msa'); \
    BiocManager::install('ggtree'); \
    install.packages('tidytree'); \
    install.packages('seqinr')"

# Set workdir
COPY . /pipeline

# Default to bash so you can override the command easily in Nextflow
CMD ["bash"]

