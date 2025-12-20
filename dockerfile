FROM python:3.9
WORKDIR /app
COPY requirements.txt .
# Increase timeout to 300 seconds (5 minutes) or more
RUN pip install --default-timeout=300 -r requirements.txt
# Copier tout le code source
COPY thermodynamique.py .
COPY evaporateurs.py .
COPY cristallisation.py .
COPY optimisation.py .
COPY main.py .
COPY app_flask.py .
EXPOSE 5000
CMD ["python", "app_flask.py"]