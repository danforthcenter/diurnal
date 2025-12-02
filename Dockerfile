FROM rocker/shiny-verse

RUN apt update && apt install -y default-jre && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('shiny', 'shinyBS', 'ggplot2', 'plotly', 'jsonlite', 'data.table', 'DT', 'stringr', 'remotes')); remotes::install_github('daqana/dqshiny')"

RUN rm -fr /srv/shiny-server

ADD diurnal /srv/shiny-server
