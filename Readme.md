# CCM Alphafold/Rosettafold Bootcamp

This is the repo contains the dockerfile to build a container that contains most of the software needed for the bootcamp. 
Data needed is not included here and will be distributed before the bootcamp.

To use this container you will need to install docker. Please see their [instructions](https://docs.docker.com/get-docker/) 
for your operating system for installing Docker. 

There are several pieces of software included in this container, theses are:

+ openstructure
+ python (and needed packages, including juypyter lab)
+ R (including IRkernel)
+ Mustang for structure alignment

All these will be accessed through jupyterlab. 

To complete the bootcamp you will also need to install pymol.


Use this repository if you have problems getting the docker container from the docker hub. For that run 

```bash
docker pull celalp/af_bootcamp:latest
```

If you encounter errors you can build the container using:

```bash
# clone the repository
git clone https://github.com/celalp/af_bootcamp.git
cd af_bootcamp
Docker build -t celalp/af_bootcamp:latest .
```

We want the container to have access to our data files so go to the directory with all the alphafold data 

```bash
cd datadir

docker run -it --rm -v "$(pwd)":/home -p 127.0.0.1:8900:8900 af_bootcamp:latest
```


## Pymol installation 

Because pymol has a graphical user interface that we will be using installing it through Docker is challenging and may need
computer specific tweaking. Therefore it is not included in the container

Linux you can use your OS's package manager such as `apt`: 
```bash
sudo apt-get install -y pymol
```

This is the recommended method, or you can download the binary and unpack it:

```bash
wget https://pymol.org/installers/PyMOL-2.3.4_121-Linux-x86_64-py37.tar.bz2
```

For Windows installatin see [here](https://pymolwiki.org/index.php/Windows_Install)

Mac Os:

The easiset way to install pymol is through [brew](https://brew.sh/), see their instructions for installing brew on your
computer

```bash
brew install pymol
```
