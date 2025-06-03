# Bioinformatics Docker Setup: `cu_docker_v1.2`

This Docker environment provides a pre-configured Ubuntu 24.04 system with all the bioinformatics software needed for this course. It is packaged as `cu_docker_v1.2.tar`.

> ⚠️ This setup guide is tailored for macOS with Apple Silicon (ARM). Other platforms may need modifications.

---

## 📦 Requirements

- macOS (Apple Silicon preferred)
- Docker Desktop: [Download for Mac (Apple Chip)](https://www.docker.com/get-started/)
- `cu_docker_v1.2.tar` image file

---

## 🧪 Step 1: Verify the Docker Image

Verify the integrity of the image:

```bash
md5 cu_docker_v1.2.tar
# Expected: f0e18012acefece3c1643bee4d4cbb65
```

---

## 🚀 Step 2: Set Up Docker

1. **Install Docker** from the official site.
2. **Start Docker** — the whale icon should appear in the menu bar.

---

## 📥 Step 3: Load the Docker Image

```bash
docker load -i /path/to/cu_docker_v1.2.tar
```

---

## 🧱 Step 4: Create and Run Your Container

```bash
docker run -it --name mybioenv cu_docker_v1.2:latest
```

Replace `mybioenv` with any container name you prefer.

---

## 📂 Step 5: Copy Project Files into Container

From your terminal:

```bash
docker cp /path/to/local/file mybioenv:/home/projects/
```

---

## 🔧 Step 6: Use the Container

To enter the container:

```bash
docker exec -it mybioenv /bin/bash
```

From inside the container, you can:

```bash
mkdir /home/projects
cd /home/projects
# Run your bioinformatics workflows
```

---

## 🛑 Step 7: Stop the Container

```bash
docker stop mybioenv
```

---

## 🔁 Step 8: Resume the Container Later

```bash
docker start mybioenv
docker exec -it mybioenv /bin/bash
```

---

## 🧹 Step 9: Clean Up (Final Use)

```bash
docker stop mybioenv
docker rm mybioenv
docker image rm cu_docker_v1.2:latest
```

---

## 🔍 Helpful Commands

- List images: `docker images`
- List containers: `docker ps -a`
- Help: `docker --help`

---

## 📎 Credits

Created by: Zitong He  
Documented by: [Your Name]
