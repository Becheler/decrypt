# decrypt

Cryptic biodiversity or subpopulations? Check the tutorials at: https://becheler.github.io/pages/applications.html


## Installation

The following instructions allow to install decrypt on an Ubuntu environment.

### 1 Get the dependencies

We first need to install three main dependencies.
Normally, if the dependencies are not found on your system, they should be downloaded
and installed locally, but problems sometimes happen, so you may want to install
them system-wide.

The Geospatial Data Abstraction Library (GDAL) is useful to represent a
spatially explicit landscapes. Boost is an important code resource for C++ dev
that covers a wide range of problems. SQLite3 is a lightweight database software
we use to store intermediary results.

Open a terminal and type
```
sudo apt-get install libgdal-dev libboost-all-dev sqlite3
```

### 2- Get Decrypt source code

Two options are possible here, up to you:

#### You are not interested in future Decrypt features

You may just want to download the [latest release of the project](https://github.com/Becheler/decrypt/releases).

#### You want to benefit from further developments

Then you should prefer to clone the github project, so you can update the decrypt pipeline
whenever you want. To clone the project, you need git, so first be sure that git is installed in your system.

To do so, simply open a terminal, type ``git --version`` and press Enter.

If the terminal answers something like ``git version 2.17.1``, it's good: git is already installed.
If it is not the case, then check the [git website](https://git-scm.com/) for proper installation.

Once you have successfully installed git, then open a terminal, chose a suitable
folder in your file system, and type:

```
git clone https://github.com/Becheler/decrypt.git
```

### 3 - Build and install

Create a directory ```sandbox``` somewhere on your computer. We will use it
as both an install location and an application directory the time for us to package the software better.

At this point, go to the ``decrypt`` directory, build the project, run the tests and install
the project to the ```sandbox``` location

```
cd decrypt
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=path/to/the/sandbox/directory
cmake --build . --target install --config Release
```

The structure of the sandbox directory is:
```
sandbox/
    |--- decrypt/
            |--- bin/
            |--- example/
```

### 4 - Tests

In sandbox directory:
```
mkdir output

chmod u+x decrypt/bin/model_1
./decrypt/bin/model_1 --config decrypt/example/config_1.ctl --landscape decrypt/example/australia_precipitation_6032.tif

chmod u+x decrypt/bin/model_2
./decrypt/bin/model_2 --config decrypt/example/config_1.ctl --landscape decrypt/example/australia_precipitation_6032.tif

chmod +x decrypt/animate.R
./animate.R output/N.tif 100 
```
