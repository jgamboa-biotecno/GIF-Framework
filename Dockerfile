# GIF Framework — Genomic Intelligence Framework
# Copyright 2026 Javier Gamboa
# Licensed under the Apache License, Version 2.0
FROM python:3.10-slim

LABEL maintainer="Javier Gamboa <jgamboa@biotecno.org>"
LABEL description="GIF Framework — Genomic Intelligence Framework"
LABEL license="Apache-2.0"

# System dependencies (blast+ for AMR detection, prodigal for gene prediction,
# git required by chewBBACA for schema download)
RUN apt-get update && apt-get install -y --no-install-recommends \
    ncbi-blast+ \
    prodigal \
    wget \
    git \
    && rm -rf /var/lib/apt/lists/*

# Bioinformatics Python tools
RUN pip install --no-cache-dir chewbbaca

COPY . /app
WORKDIR /app
RUN pip install --no-cache-dir .

ENTRYPOINT ["gif"]
CMD ["--help"]
