FROM rocker/rstudio:latest

ARG path="/home/rstudio"
WORKDIR ${path}

RUN mkdir database_comparator
RUN chown -R rstudio /home/rstudio/database_comparator

COPY DEFAULT_config_file.txt database_comparator/
#ncbi-blast
RUN apt-get update && apt-get -y install ncbi-blast+
#R
RUN install2.r --error \
    reticulate \
    png \
    jsonlite


# Python3
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-dev \
    python3-venv \
    && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*


# Install Python packages using reticulate
RUN pip install database_comparator


WORKDIR /home/rstudio/database_search
EXPOSE 8787
