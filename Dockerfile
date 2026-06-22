FROM python:3.11-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    procps \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
    numpy==1.23.5 \
    matplotlib==3.5.1 \
    pysam==0.19.1

WORKDIR /app

COPY swave /app/swave

ENV PYTHONPATH="/app/swave"

CMD ["python3"]