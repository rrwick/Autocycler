FROM mambaorg/micromamba:latest

# Create an "autocycler" environment with long-read assemblers.
RUN micromamba create -y -n autocycler python=3.12 -c conda-forge -c bioconda curl && \
    for tool in canu flye lja metamdbg miniasm minimap2 minipolish necat nextdenovo nextpolish racon raven-assembler wtdbg; do \
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
    tar -xzf autocycler.tar.gz && \
    mv autocycler /opt/conda/envs/autocycler/bin/ && \
    rm environment.yml autocycler.tar.gz

# Set default command to Autocycler help and show missing/installed tools.
CMD ["/bin/sh", "-c", "autocycler --help; echo \"\nInstalled tools:\"; cat /tmp/installed_tools.txt; echo \"\nMissing tools:\"; cat /tmp/missing_tools.txt; echo"]
