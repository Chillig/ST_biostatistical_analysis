# Spatial Transcriptomics Density Clustering and Correlation Analysis

This repository provides tools to analyse spatial transcriptomics data by identifying 
gene-of-interest (GOI) clusters, quantifying responder gene interactions, and 
computing weighted correlations. It runs the [analysis](https://doi.org/10.1038/s41467-022-35319-w) 
shown in Figure 5 D-H.  <br>

For detailed background and workflow explanation, 
see [Information_analysis.md](./Information_analysis.md).  <br>

Here we show how to clone ST_biostatistical_analysis → build the docker image → run Jupyter. 
See [Detailed Workflow](#detailed-workflow) for more information on how to get started. <br>
```bash
# Clone repo
git clone git@github.com:Chillig/ST_biostatistical_analysis.git
cd ST_biostatistical_analysis

# Build image
docker build -t spatial_clustering:latest .

# Run container with Jupyter
docker run -m 24g -it --rm \
  -p 8888:8888 \
  -v "$(pwd)":/app \
  -v /path/to/data:/data \
  -v "$(pwd)/results":/app/results \
  spatial_clustering
```


## Detailed Workflow
Get GitHub repository
```bash
git clone git@github.com:Chillig/ST_biostatistical_analysis.git
cd ST_biostatistical_analysis
```

## Quick overview:
1. Install Docker → via package manager (Linux) or [Docker Desktop](https://www.docker.com/products/docker-desktop) (MacOS/Windows)
2. Build the image → using the provided [Dockerfile](Dockerfile)
3. Run the analysis → either via Jupyter Notebook or directly with Python

Start Docker Desktop and verify Docker installation in terminal:
```bash
docker --version
```

### Build Docker Image
This repository contains a [Dockerfile](Dockerfile) and environment 
specification [python38_conda_env.yml](python38_conda_env.yml) to build the Docker Image. <br>
To build the Docker Image, run from inside _ST_Gene_Responders_Clustering_ directory:
```bash
docker build -t spatial_clustering:latest .
```
Here, _spatial_clustering:latest_ is the _<repository>:<tag>_ name,
but you can use any name you prefer. <br>
Verify that image exists:
```bash
docker images
```

### Run Docker Container
A Docker Container is a running instance of a Docker image.
You can run it with Jupyter Notebook support and mounted directories:
```bash
docker run -m 24g -it --rm \ 
  -p 8888:8888 \
  -v "$(pwd)":/app \
  -v /path/to/data:/data \
  -v /path/to/scripts:/scripts \
  -v /path/to/output:/output \
  spatial_clustering
```
Open [Run_clustering.ipynb](./spatial_correlation/Run_clustering.ipynb) to follow the analysis pipeline step by step. <br>

Alternatively, you can run the analysis directly:
```bash
docker run -m 24g -it --rm \ 
  -v "$(pwd)":/app \
  -v /path/to/data:/data \
  -v /path/to/scripts:/scripts \
  -v /path/to/output:/output \
  spatial_clustering python /scripts/main.py
```
Explanation:
- -m 24g → limit container memory to 24 GB.
- -it → interactive terminal.
- -p 8888:8888 → maps host port 8888 → container port 8888 (for Jupyter or Flask).
- -v "$(pwd)":/app → mounts your current dir into /app in container.
- --rm → remove container after exit (optional).
- -v host_dir:container_dir → mounts directories.
- spatial_clustering → the image you built earlier.