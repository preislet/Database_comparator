FROM rocker/rstudio:latest

ARG path="/home/rstudio"
WORKDIR ${path}

RUN mkdir database_search
RUN chown -R rstudio /home/rstudio/database_search

#ncbi-blast
RUN apt-get update && apt-get -y install ncbi-blast+
#R
RUN install2.r --error \
    reticulate \
    png \
    jsonlite \

# Install Python packages using reticulate
RUN R -e "reticulate::py_install('Database_comparator')"


WORKDIR /home/rstudio/database_search
EXPOSE 8787