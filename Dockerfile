FROM ubuntu:24.10
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get -y update && apt-get -y upgrade && apt-get -y install r-base r-base-dev
RUN apt-get -y install gdebi-core
RUN apt-get -y install wget
RUN wget https://download2.rstudio.org/server/jammy/amd64/rstudio-server-2023.12.1-402-amd64.deb
RUN gdebi -n rstudio-server-2023.12.1-402-amd64.deb
RUN useradd rstudio -p "\$y\$j9T\$/.6YKeUOB4ifaPjuG/xaC1\$0162SW98NtTo5c6I7uXbwlNlKGuu9LTcUanCzz6DF/C" -d /home/rstudio -m
RUN apt-get update
RUN apt-get install -y \
    apt-transport-https \
    ca-certificates \
    curl \
    gnupg \
    lsb-release
RUN apt-get update && apt install -y libudunits2-dev libgdal-dev
RUN apt-get update
RUN apt-get -y install gfortran
RUN apt-get -y install build-essential
RUN apt-get -y install fort77
RUN apt-get -y install xorg-dev
RUN apt-get -y install liblzma-dev  libblas-dev gfortran
RUN apt-get -y install gobjc++
RUN apt-get -y install aptitude
RUN apt-get -y install libbz2-dev
RUN apt-get -y install libpcre3-dev
RUN aptitude -y install libreadline-dev
RUN apt-get -y install libcurl4-openssl-dev
RUN apt install -y build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev
RUN apt-get install -y libcurl4-openssl-dev libssl-dev
RUN apt-get install -y libgit2-dev
RUN apt-get install -y libharfbuzz-dev
RUN apt-get install -y libfribidi-dev
RUN apt-get install -y cmake
RUN apt-get install -y libcairo2-dev
RUN  Rscript -e 'install.packages(c("dplyr", "reshape2","BiocManager", "readxl"), dependencies = TRUE)'
RUN  Rscript -e 'BiocManager::install("crisprScore")'
RUN  Rscript -e 'install.packages("devtools")'
RUN  Rscript -e 'library(devtools); install_github("crisprVerse/crisprScoreData"); install_github("crisprVerse/crisprScore")'
EXPOSE 8787
ENTRYPOINT ["tail"]
CMD ["-f","/dev/null"]