FROM golob/dada2:1.14.1.ub.1804

MAINTAINER rdemko2332@gmail.com

WORKDIR /usr/bin/

RUN Rscript -e 'install.packages("data.table")'

COPY /bin/* /usr/bin/

Run chmod +x *

WORKDIR /work

ENV PERL5LIB=/usr/bin/
