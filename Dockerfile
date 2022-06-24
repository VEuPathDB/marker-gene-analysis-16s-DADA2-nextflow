FROM golob/dada2:1.14.1.ub.1804

MAINTAINER rdemko2332@gmail.com

WORKDIR /usr/bin/

RUN Rscript -e 'install.packages("data.table")' \
  && Rscript -e 'install.packages("optparse")'

COPY /bin/* /usr/bin/

RUN chmod +x *

WORKDIR /work
