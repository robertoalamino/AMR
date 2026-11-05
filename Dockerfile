FROM debian:bookworm-slim@sha256:734a72b3d6142a71993b4bc580ca01d8ca687f77a2c422d96a961d3376b179bc

ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

# Use maintained Debian repositories.
RUN apt-get -o Acquire::Retries=5 update \
    && apt-get -o Acquire::Retries=5 install -y --no-install-recommends \
       bash \
       gcc \
       libc6-dev \
       gawk \
       coreutils \
       diffutils \
       findutils \
       grep \
       sed \
       git \
       tar \
       gzip \
       gnuplot-nox \
       fonts-liberation \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /amr

COPY . .

CMD ["bash", "scripts/check-all.sh"]