# Builder base
FROM mambaorg/micromamba:latest as builder
USER root
COPY --chmod=755 install_tools.sh /usr/local/bin/
USER $MAMBA_USER

# Create separate builders to facilitate developing this Docker image
# Should run in parallel with modern BuildKit.
# Enable multiple processes with Podman: podman build --tag autocycler --jobs 0
FROM builder as builder_flye
RUN install_tools.sh flye
FROM builder as builder_canu
RUN install_tools.sh canu bioconda::canu
FROM builder as builder_metamdbg
RUN install_tools.sh metamdbg
FROM builder as builder_necat
RUN install_tools.sh necat
FROM builder as builder_minipolish
RUN install_tools.sh minipolish miniasm bioconda::minipolish racon
FROM builder as builder_verkko
RUN install_tools.sh Verkko
FROM builder as builder_nextdenovo
RUN install_tools.sh nextdenovo nextdenovo nextpolish
FROM builder as builder_raven
RUN install_tools.sh raven raven-assembler
FROM builder as builder_wtdbg
RUN install_tools.sh wtdbg
FROM builder as builder_plassembler
RUN install_tools.sh plassembler

# Actual Autocycler build starts here
FROM mambaorg/micromamba:latest

# Install dependencies (libgomp1 for lja)
USER root
RUN apt-get update && \
    apt-get install -y vim curl libgomp1 && \
    rm -rf /var/lib/apt/lists/*

# Switch back to default user
USER $MAMBA_USER

# Create a default autocycler env
RUN micromamba create -y -n autocycler python && micromamba clean --all --yes
ENV PATH="/opt/conda/envs/autocycler/bin:$PATH"

# Download and install Autocycler (binary and scripts) from the latest release.
RUN curl -L -o autocycler.tar.gz $(curl -s https://api.github.com/repos/rrwick/Autocycler/releases/latest | grep "browser_download_url.*linux-x86_64-musl.*tar.gz" | cut -d '"' -f 4) && \
    tar -xzf autocycler.tar.gz && \
    mv autocycler /opt/conda/envs/autocycler/bin/ && \
    rm environment.yml autocycler.tar.gz

# Copy certain tools from public containers
COPY --from=docker.io/staphb/minimap2 /usr/local/bin/minimap2 /usr/local/bin/minimap2
COPY --from=docker.io/troder/lja /opt/LJA/bin/lja /usr/local/bin/lja

# Copy the installed tools from each builder stage
COPY --chown=$MAMBA_USER_ID:$MAMBA_USER_GID --chmod=755 --from=builder_raven /opt/conda/envs/raven /opt/conda/envs/raven
COPY --chown=$MAMBA_USER_ID:$MAMBA_USER_GID --chmod=755 --from=builder_flye /opt/conda/envs/flye /opt/conda/envs/flye
COPY --chown=$MAMBA_USER_ID:$MAMBA_USER_GID --chmod=755 --from=builder_canu /opt/conda/envs/canu /opt/conda/envs/canu
COPY --chown=$MAMBA_USER_ID:$MAMBA_USER_GID --chmod=755 --from=builder_metamdbg /opt/conda/envs/metamdbg /opt/conda/envs/metamdbg
COPY --chown=$MAMBA_USER_ID:$MAMBA_USER_GID --chmod=755 --from=builder_necat /opt/conda/envs/necat /opt/conda/envs/necat
COPY --chown=$MAMBA_USER_ID:$MAMBA_USER_GID --chmod=755 --from=builder_minipolish /opt/conda/envs/minipolish /opt/conda/envs/minipolish
COPY --chown=$MAMBA_USER_ID:$MAMBA_USER_GID --chmod=755 --from=builder_verkko /opt/conda/envs/Verkko /opt/conda/envs/Verkko
COPY --chown=$MAMBA_USER_ID:$MAMBA_USER_GID --chmod=755 --from=builder_nextdenovo /opt/conda/envs/nextdenovo /opt/conda/envs/nextdenovo
COPY --chown=$MAMBA_USER_ID:$MAMBA_USER_GID --chmod=755 --from=builder_wtdbg /opt/conda/envs/wtdbg /opt/conda/envs/wtdbg
# pbassembler: huge install -.-
COPY --chown=$MAMBA_USER_ID:$MAMBA_USER_GID --chmod=755 --from=builder_plassembler /opt/conda/envs/plassembler /opt/conda/envs/plassembler

# Set default command to Autocycler help and show missing/installed tools.
CMD ["/bin/sh", "-c", "autocycler --help"]
