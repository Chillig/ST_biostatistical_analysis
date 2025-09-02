FROM continuumio/miniconda3

WORKDIR /app

# Copy only environment file
COPY python38_conda_env.yml .

# Install mamba for faster env creation
RUN conda install -n base -c conda-forge mamba \
    && mamba env create -f python38_conda_env.yml \
    && conda clean -afy

# Activate environment by default
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "analysis_env"]

# Expose Jupyter port
EXPOSE 8888

# Default: start JupyterLab
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root"]

