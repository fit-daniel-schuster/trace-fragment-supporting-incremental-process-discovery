FROM python:3.10

WORKDIR /usr/src/app

COPY ./requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt
RUN pip install cvxopt

# setup default for variables
ENV PYTHONPATH "${PYTHONPATH}:/usr/src/app"
ENV LOG_FILE 'BPI_Challenge_2012.xes'
ENV EVENT_LOG_DIRECTORY "/usr/data/"
ENV RESULT_DIRECTORY "/usr/results/"
ENV REL_SUPP_INIT_MODEL "0.01"
ENV NOISE_THRESHOLD_INIT_MODEL "0"

COPY ./scripts ./scripts

CMD [ "python", "./scripts/experiments_artificial.py" ]
