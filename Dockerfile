FROM sagemath/sagemath:latest

# Ignore APT warnings about not having a TTY
ENV DEBIAN_FRONTEND noninteractive

RUN pwd

COPY . ${HOME}


# RUN cd /home/sage
RUN sage  -pip install --upgrade --no-index -v .
# RUN cd /home



