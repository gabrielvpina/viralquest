FROM python:3.12-slim AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    xz-utils \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY pyproject.toml MANIFEST.in ./
COPY viralquest ./viralquest

RUN unxz viralquest/bin/*.xz && \
    chmod +x viralquest/bin/*

RUN pip install --no-cache-dir .

FROM python:3.12-slim

RUN useradd --create-home --shell /bin/bash appuser

COPY --from=builder /usr/local/lib/python3.12/site-packages /usr/local/lib/python3.12/site-packages
COPY --from=builder /usr/local/bin /usr/local/bin

RUN chown -R appuser:appuser /usr/local/lib/python3.12/site-packages
RUN chown -R appuser:appuser /usr/local/bin

USER appuser
WORKDIR /home/appuser

ENTRYPOINT ["viralquest"]
CMD ["--help"]
