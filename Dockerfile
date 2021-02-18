FROM sagemath/sagemath:latest

# Ignore APT warnings about not having a TTY
ENV DEBIAN_FRONTEND noninteractive

RUN pwd

# ImageMagick and Graphviz (for graph examples)
RUN sudo apt-get -q update && sudo apt-get -qy dist-upgrade
RUN sudo apt-get -qy install imagemagick
RUN sudo apt-get -qy install graphviz

COPY . ${HOME}


# RUN cd /home/sage
RUN sage  -pip install --upgrade --no-index -v .
# RUN cd /home



