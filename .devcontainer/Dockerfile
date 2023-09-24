# Use an official Clojure image as a parent image
FROM ubuntu:22.04

# Ensure non-interactive installations (avoid prompts)
ENV DEBIAN_FRONTEND=noninteractive

# Update the repositories and install required packages
RUN apt-get update && apt-get install -y \
    # Install OpenJDK 17
    openjdk-17-jdk \
    # Install Clojure (use the official script for this)
    && apt-get install -y curl \
    && curl -O https://github.com/clojure/brew-install/releases/download/1.11.1.1413/linux-install.sh \
    && chmod +x linux-install.sh \
    && ./linux-install.sh \
    # Install libassimp-dev
    && apt-get install -y libassimp-dev rlwrap \
    # Cleanup
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && rm linux-install.sh

# Set the working directory in docker
WORKDIR /workspace

# Copy the current directory contents into the container at /workspace
COPY . /workspace

RUN mkdir -p /var/lib/plexus

# Install system dependencies (if any)

# You can add any other tools or utilities you need

# Expose the port the app runs on (this is optional and depends on your application)
EXPOSE 8080

CMD ["tail", "-f", "/dev/null"]
