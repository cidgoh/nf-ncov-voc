FROM continuumio/miniconda3:latest
LABEL authors="Muhammad Zohaib Anwar, Ivan Gill, Madeline Iseminger, Anoosha Sehar & William Hsiao" \
      description="Docker image containing all requirements for the nf-ncov-voc"

COPY environments/extras.yml /extras.yml
COPY environments/environment.yml /environment.yml
RUN apt-get update && apt-get install -y curl g++ git make procps && apt-get clean -y
RUN /opt/conda/bin/conda env create -f /environment.yml
RUN /opt/conda/bin/conda env create -f /extras.yml && /opt/conda/bin/conda clean -a
ENV PATH=/opt/conda/envs/nf-ncov-voc/bin:$PATH
ENV PATH=/opt/conda/envs/nf-ncov-voc-extra/bin:$PATH
