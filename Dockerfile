FROM sagemath/sagemath:latest

# Ignore APT warnings about not having a TTY
ENV DEBIAN_FRONTEND noninteractive

RUN pwd


# RUN cd /home/sage
RUN sage -pip install .
# RUN cd /home


COPY . ${HOME}
