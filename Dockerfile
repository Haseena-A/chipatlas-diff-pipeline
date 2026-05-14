# chipatlas-diff-pipeline – Dockerfile
# Multi-stage build: slim runtime image

# Stage 1: build wheels
FROM python:3.11-slim AS builder
WORKDIR /build
RUN apt-get update && apt-get install -y --no-install-recommends gcc g++ && rm -rf /var/lib/apt/lists/*
COPY requirements.txt .
RUN pip install --upgrade pip && pip wheel --no-cache-dir --wheel-dir /wheels -r requirements.txt

# Stage 2: runtime
FROM python:3.11-slim AS runtime
LABEL description="ChIP-Atlas differential binding/methylation analysis pipeline"
WORKDIR /app
RUN useradd --create-home --shell /bin/bash analyst
COPY --from=builder /wheels /wheels
RUN pip install --no-cache-dir --no-index --find-links /wheels /wheels/*.whl && rm -rf /wheels
COPY src/ /app/src/
COPY cli.py /app/
COPY configs/ /app/configs/
COPY examples/ /app/examples/
RUN mkdir -p /app/results /home/analyst/.chipatlas_cache && chown -R analyst:analyst /app/results /home/analyst/.chipatlas_cache
USER analyst
ENV PYTHONPATH=/app HOME=/home/analyst
ENTRYPOINT ["python", "/app/cli.py"]
CMD ["--help"]
