# configs/docker/validation.Dockerfile
FROM python:3.10-slim

WORKDIR /app
COPY modules/validation/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY modules/validation /app
ENTRYPOINT ["python", "/app/validation.py"]