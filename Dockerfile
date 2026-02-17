# Use Python 3.11 slim image
FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    GUANACO_CONFIG=/app/config.json \
    GUANACO_DATA_DIR=/app/data \
    GUANACO_MAX_CELLS=20000 \
    GUANACO_BACKED_MODE=true \
    GUANACO_WORKERS=2 \
    GUANACO_THREADS=2 \
    GUANACO_TIMEOUT=300

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    && rm -rf /var/lib/apt/lists/*

# Copy the entire application
COPY . .

# Install package and runtime dependencies from pyproject.toml
RUN pip install --no-cache-dir .

# Expose the port the app runs on
EXPOSE 8080
CMD ["/bin/sh", "-c", "gunicorn --bind 0.0.0.0:8080 --workers ${GUANACO_WORKERS} --threads ${GUANACO_THREADS} --timeout ${GUANACO_TIMEOUT} guanaco.main:server"]
