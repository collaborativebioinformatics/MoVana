FROM python:3

# Install required Python libraries
RUN pip install numpy pandas matplotlib

# Set the working directory
WORKDIR ~/MoVana

# Copy the WDL script into the container
COPY . ~/MoVana

# Set the entry point
ENTRYPOINT ["python3"]
