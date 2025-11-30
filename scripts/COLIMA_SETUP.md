# Using Colima Instead of Docker Desktop

Colima is a lightweight alternative to Docker Desktop that runs Docker in a VM without the heavy GUI application.

## Installation

```bash
brew install colima
```

## Starting Colima

```bash
# Start Colima (creates VM and starts Docker daemon)
colima start

# Check status
colima status

# Stop Colima
colima stop
```

## Auto-start on Login (Optional)

```bash
# Start Colima as a service (starts automatically on login)
brew services start colima

# Stop auto-start
brew services stop colima
```

## Using Docker with Colima

Once Colima is started, Docker works exactly the same:

```bash
# Verify Docker is working
docker ps

# Build containers
docker build -t myimage .

# Run containers
docker run myimage
```

## Switching Between Docker Desktop and Colima

If you have both installed, you can switch contexts:

```bash
# Use Colima
docker context use colima

# Use Docker Desktop
docker context use desktop-linux

# List all contexts
docker context ls
```

## Building Optimization Container

With Colima running, build the container:

```bash
bash scripts/optimization/build_with_docker.sh rna-map-optimization.sif
```

Then transfer to cluster:

```bash
scp rna-map-optimization.sif user@cluster:/path/to/destination/
```

## Advantages of Colima

- **Lightweight**: No GUI, minimal resource usage
- **Fast startup**: VM starts quickly
- **Command-line only**: Perfect for automation
- **Free**: No Docker Desktop subscription needed

## Troubleshooting

### Colima won't start
```bash
# Check if VM is running
colima status

# Restart Colima
colima stop
colima start
```

### Docker commands fail
```bash
# Make sure Colima is running
colima status

# Switch to Colima context
docker context use colima

# Verify Docker socket
docker ps
```

### Out of disk space
```bash
# Check Colima disk usage
colima ssh -- df -h

# Increase disk size (requires stopping Colima)
colima stop
colima start --disk 40  # Increase to 40GB
```

