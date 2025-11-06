# S1. Load the Docker image (run once only)
sudo docker load -i bootstrap_readout_v1.0.tar

# S2. Set the working directory (local path on your host machine)
parent_dir="$(dirname "$(pwd)")"
work_dir="$parent_dir"

# S3. Run the docker container
sudo docker run -it \
  -v $work_dir:/workspace \
  -w /workspace \
  bootstrap_readout:v1.0 \
  /bin/bash
