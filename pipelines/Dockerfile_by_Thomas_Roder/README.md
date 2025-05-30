# **Advanced Autocycler Docker Image**

## **Overview**

An alternative Docker image for Autocycler.

**Key differences:**

* Assemblers installed into separate Micromamba environments
* Successfully installs all targeted assemblers
* Reduced Build Speed: Docker's multi-stage builds are parallelizable
* Disadvantages: increased image size and Dockerfile complexity

## **Building the Image**

Ensure Dockerfile and `install_tools.sh` are in the same directory.

```bash
# Using Docker  
docker build --tag autocycler .
# Using Podman (with parallel jobs)  
podman build --tag autocycler --jobs 0 .
```

## **Running the Container**

```bash
# Interactive Session
docker run --rm -it autocycler bash
# Mounting Data
docker run --rm -it -v /path/to/your/data:/data autocycler bash
```

## **Installed Software & Environment Activation**

Autocycler itself and its helper scripts are installed in the autocycler Conda environment, which is added to the system
PATH. Thus, autocycler can be run directly.  
The following table lists the installed assemblers/tools and the commands required to activate their respective
environments if needed for manual use in an interactive session:

| Assembler                                    | Activation Command              |
|:---------------------------------------------|:--------------------------------|
| Autocycler                                   | no activation needed            |
| Flye                                         | micromamba activate flye        |
| Canu                                         | micromamba activate canu        |
| MetaMDBG                                     | micromamba activate metamdbg    |
| NECAT                                        | micromamba activate necat       |
| Minipolish (incl. Miniasm, Racon, Any2fasta) | micromamba activate minipolish  |
| Verkko                                       | micromamba activate Verkko      |
| NextDenovo (incl. NextPolish)                | micromamba activate nextdenovo  |
| Raven                                        | micromamba activate raven       |
| WTDBG2 (wtdbg)                               | micromamba activate wtdbg       |
| Plassembler                                  | micromamba activate plassembler |
| LJA                                          | no activation needed            |
| Seqtk (not an assembler)                     | no activation needed            |
| Minimap2 (not an assembler)                  | no activation needed            |
