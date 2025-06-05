FROM jupyter/datascience-notebook:python-3.11

# Set working directory inside container
WORKDIR /home/database_comparator

# Copy project folders into the container
COPY notebooks/ notebooks/
COPY tests/ tests/

# System dependencies (e.g. BLAST for your library)
USER root
RUN apt-get update && apt-get install -y ncbi-blast+ && apt-get clean

# Make target folders writable by notebook user
RUN chmod -R a+w /home/database_comparator/notebooks \
    && chmod -R a+w /home/database_comparator/tests

# Revert to Jupyter's non-root user
USER $NB_UID

# Install your Python package from PyPI
RUN pip install --no-cache-dir database_comparator openpyxl

# Install spreadsheet editor for CSV/TSV
RUN pip install jupyterlab-spreadsheet-editor
RUN pip install jupyterlab-spreadsheet


# Start JupyterLab with no token for convenience
CMD ["start.sh", "jupyter", "lab", "--NotebookApp.token=''", "--NotebookApp.password=''"]
