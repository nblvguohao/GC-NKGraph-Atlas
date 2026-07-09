# GC-NKGraph-Atlas — Reproducible Pipeline Container
#
# This Dockerfile provides a minimal, self-contained environment for running
# the mechanism-card-driven target-discovery pipeline.  It ships the full
# codebase plus the interactive web playground, and defaults to SYNTHETIC
# DATA MODE so reviewers can verify the pipeline executes end-to-end without
# downloading any real data.
#
# Quick start (synthetic data — no network after build):
#   docker build -t gc-nkgraph-atlas .
#   docker run --rm gc-nkgraph-atlas
#
# With real data (mount your data/ directory):
#   docker run --rm -v /path/to/data:/workspace/data gc-nkgraph-atlas \
#       python src/pipeline.py --config configs/experiment_config.yaml
#
# Web playground:
#   docker run --rm -p 8080:8080 gc-nkgraph-atlas web

# ---------------------------------------------------------------------------
# Stage 1: Builder — conda environment from environment.yml
# ---------------------------------------------------------------------------
FROM continuumio/miniconda3:24.1.2 AS builder

WORKDIR /workspace

# Copy dependency files first for layer caching
COPY environment.yml requirements.txt pyproject.toml ./

# Create conda environment
RUN conda env create -f environment.yml -n gc-nkgraph && \
    conda clean -afy

# ---------------------------------------------------------------------------
# Stage 2: Runtime — slim, production-ready
# ---------------------------------------------------------------------------
FROM continuumio/miniconda3:24.1.2 AS runtime

WORKDIR /workspace

# Copy conda environment from builder
COPY --from=builder /opt/conda/envs/gc-nkgraph /opt/conda/envs/gc-nkgraph

# Copy the full codebase (excludes data/ which is .gitignored)
COPY . .

# Make conda environment available by default
ENV PATH="/opt/conda/envs/gc-nkgraph/bin:$PATH"
ENV CONDA_DEFAULT_ENV=gc-nkgraph
ENV PYTHONUNBUFFERED=1

# Expose port for web playground
EXPOSE 8080

# ---------------------------------------------------------------------------
# Entrypoint: run the full pipeline in synthetic mode as a smoke test.
# This DOES NOT require any real data and completes in < 2 minutes on CPU.
# ---------------------------------------------------------------------------
COPY docker-entrypoint.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/docker-entrypoint.sh

ENTRYPOINT ["docker-entrypoint.sh"]

# Default: run synthetic pipeline + tests to verify everything works
CMD ["verify"]
