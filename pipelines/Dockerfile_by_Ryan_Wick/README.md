# Autocycler Dockerfile (by Ryan Wick)

This directory contains a Dockerfile for building a container with Autocycler and a suite of long-read assemblers. It uses micromamba to create a Conda environment and installs each tool individually, allowing the build to succeed even if some tools are unavailable on your platform.



## Build the image

From this directory, run:
```bash
docker build -t autocycler .
```



## Run the image

To run Autocycler and display which tools were successfully installed:
```bash
docker run --rm autocycler
```

To run an interactive session with a mounted directory (e.g. your data):
```bash
docker run --rm -it -v /path/to/data:/data autocycler bash
```



## Notes

* All tools are installed into the same Conda environment named `autocycler`.
* Successfully installed tools are listed in `/tmp/installed_tools.txt`. Tools that failed to install are listed in `/tmp/missing_tools.txt`.
* If needed, you can modify the Dockerfile to add/remove specific tools (see the `for tool in` loop).
