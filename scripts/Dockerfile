# -----------------------------------------------------------------------------
# Autocycler Dockerfile
#
# This Dockerfile builds an image for Autocycler. It uses micromamba to create
# a dedicated Conda environment (named "autocycler") for Autocycler and various
# long-read assemblers. Installation is done one package at a time, tolerating
# failures, since not all assemblers may be available on all platforms.
#
# Usage:
#   Build the image (from the same directory as this file):
#     docker build -t autocycler .
#
#   Run the container (to show Autocycler help and which tools are installed):
#     docker run --rm autocycler
#
#   Run an interactive session with host files mounted (e.g. mounting
#   local directory /path/to/data to /data in the container):
#     docker run --rm -it -v /path/to/data:/data autocycler bash
#
# -----------------------------------------------------------------------------

FROM mambaorg/micromamba:latest

# Create an "autocycler" environment with long-read assemblers.
RUN micromamba create -y -n autocycler -c conda-forge -c bioconda curl && \
    for tool in any2fasta canu flye lja metamdbg miniasm minimap2 minipolish3 necat nextdenovo nextpolish racon raven-assembler seqtk wtdbg; do \
        if micromamba install -y -n autocycler -c conda-forge -c bioconda "$tool"; then \
            printf " $tool" >> /tmp/installed_tools.txt; \
        else \
            printf " $tool" >> /tmp/missing_tools.txt; \
        fi \
    done && \
    printf "\n" >> /tmp/installed_tools.txt && printf "\n" >> /tmp/missing_tools.txt && \
    micromamba clean --all --yes

# Add the autocycler environment to PATH.
ENV PATH="/opt/conda/envs/autocycler/bin:$PATH"

# Download and install Autocycler (binary and scripts) from the latest release.
RUN curl -L -o autocycler.tar.gz $(curl -s https://api.github.com/repos/rrwick/Autocycler/releases/latest | grep "browser_download_url.*linux-x86_64-musl.*tar.gz" | cut -d '"' -f 4) && \
    curl -L -o autocycler_scripts.tar.gz $(curl -s https://api.github.com/repos/rrwick/Autocycler/releases/latest | grep "browser_download_url.*autocycler-helper-scripts.*tar.gz" | cut -d '"' -f 4) && \
    tar -xzf autocycler.tar.gz && \
    tar -xzf autocycler_scripts.tar.gz && \
    mv autocycler *.sh *.py /opt/conda/envs/autocycler/bin/ && \
    rm environment.yml autocycler.tar.gz autocycler_scripts.tar.gz

# Set default command to Autocycler help and show missing/installed tools.
CMD ["/bin/sh", "-c", "autocycler --help; echo \"\nInstalled tools:\"; cat /tmp/installed_tools.txt; echo \"\nMissing tools:\"; cat /tmp/missing_tools.txt; echo"]
