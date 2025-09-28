FROM jupyter/datascience-notebook:python-3.11

# Nastav pracovní adresář v rámci kontejneru
WORKDIR /home/database_comparator

# Zkopíruj projektové složky a README
COPY notebooks/ notebooks/
COPY README.md README.md

# Přepni na root kvůli instalacím a úpravám
USER root

# Nainstaluj systémové závislosti (např. BLAST)
RUN apt-get update && apt-get install -y ncbi-blast+ && apt-get clean

# Uprav práva k zápisu pro celou pracovní složku (ne jen notebooks/)
RUN chmod -R a+w /home/database_comparator

# Nainstaluj Python balíčky (stále jako root)
RUN pip install --no-cache-dir \
    database_comparator \
    openpyxl \
    jupytext \
    jupyterlab-spreadsheet-editor \
    jupyterlab-spreadsheet

# Převod README.md na interaktivní notebook (ještě jako root)
RUN jupytext --to notebook README.md --output README.ipynb && rm README.md

# Přepni zpět na JupyterLab usera
USER $NB_UID

# Spusť JupyterLab s otevřeným README.ipynb a bez tokenu
CMD ["start.sh", "jupyter", "lab", "--NotebookApp.token=''", "--NotebookApp.password=''", "--NotebookApp.default_url=/home/database_comparator/README.ipynb"]

EXPOSE 8888
