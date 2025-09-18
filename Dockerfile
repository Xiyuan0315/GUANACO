# Use Python 3.11 slim image
FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    GUANACO_CONFIG=/app/config.json \
    GUANACO_DATA_DIR=/app/data

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements first for better caching
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy the entire application
COPY . .

# Install the package in development mode
RUN pip install -e .

# Create directories for config and data
RUN mkdir -p /app/data /app/config

# Create a default config file if none exists
RUN echo '{"datasets": {}}' > /app/config.json

# Expose the port the app runs on
EXPOSE 8080

# Use gunicorn for production
CMD ["gunicorn", "--bind", "0.0.0.0:8080", "--workers", "1", "--threads", "2", "--timeout", "300", "guanaco.main:server"]