# Use R as base image to have both R and Python environments
FROM rocker/tidyverse:4.2.0

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3-pip \
    python3-dev \
    libpython3-dev \
    libbz2-dev \
    liblzma-dev \
    libpcre2-dev \
    zlib1g-dev \
    libicu-dev \
    build-essential \
    libnetcdf-dev \
    libhdf5-dev \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
COPY requirements.txt /tmp/
RUN pip3 install --no-cache-dir -r /tmp/requirements.txt

# Install Bioconductor and required R packages
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); \
    BiocManager::install(c('DEP', 'limma', 'mzR', 'MSnbase', 'vsn'), update = FALSE); \
    install.packages(c('PerformanceAnalytics', 'ncdf4'), repos='https://cloud.r-project.org/')"

# Set working directory
WORKDIR /app

# Copy project files
COPY . /app/

# Create necessary directories if they don't exist
RUN mkdir -p data/raw data/processed notebooks tests config

# Set environment variables
ENV PYTHONPATH=/app
ENV R_LIBS=/usr/local/lib/R/site-library

# Default command
CMD ["bash"] 
