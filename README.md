# Introduction

>Implementing a Unscented Kalman Filter to estimate the state of a moving object with noisy LIDAR and RADAR measurements. 
>Obtaining _RMSE_ values that are lower than the tolerance requirement.

__What I have: Data from 2 Noisy Sensors__
* A LIDAR sensor measures the position in Cartesian coordinate `(x,y)`
* A RADAR sensor measures the position and veloctiy in polar coordinate `(rho, phi, d_rho)`

__What my GOALs are: Estimate the Position, Velocity and Turning Rate__
* Assumption: `CTRV`(Constant Turn Rate and Velocity) dynamic model
* Sensor measurements noises value provided by the manufacturer, and tuning the Process noise standard deviation for both longitudinal and yaw acceleration
*  Get the position `(px,py)` in Cartesian coordinate, the velocity magnitude `v`, the yaw angle `yaw`(rad), and yaw rate `yaw_d` (rad/s)

# Rusult
In the demo video, LIDAR measurements are `red circles`, RADAR measurements are `blue circles` with an arrow pointing in the direction of the observed angle, and estimation markers are `green triangles`, and px, py, vx, and vy RMSE values from Unscented Kalman Filter are within __[0.11, 0.11, 0.52, 0.52]__, and even lower than the previous EKF project.

>Demo Video ☟

[![Video](http://img.youtube.com/vi/TYSbhb8Jkcs/0.jpg)](http://www.youtube.com/watch?v=TYSbhb8Jkcs "Unscented Kalman Filter")


# Download
Simulator which can be downloaded [here](https://github.com/udacity/self-driving-car-sim/releases).

This repository includes two files that can be used to set up and intall [uWebSocketIO](https://github.com/uWebSockets/uWebSockets) for either Linux or Mac systems. For windows you can use either Docker, VMware, or even [Windows 10 Bash on Ubuntu](https://www.howtogeek.com/249966/how-to-install-and-use-the-linux-bash-shell-on-windows-10/) to install uWebSocketIO. Please see [this concept in the classroom](https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/f758c44c-5e40-4e01-93b5-1a82aa4e044f/concepts/16cf4a78-4fc7-49e1-8621-3450ca938b77) for the required version and installation scripts.

Once the install for uWebSocketIO is complete, the main program can be built and ran by doing the following from the project top directory.

1. mkdir build
2. cd build
3. cmake ..
4. make
5. ./UnscentedKF

# Other Important Dependencies
* cmake >= 3.5
* All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1 (Linux, Mac), 3.81 (Windows)
* Linux: make is installed by default on most Linux distros
* Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
* Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
* Linux: gcc / g++ is installed by default on most Linux distros
* Mac: same deal as make - [install Xcode command line tools](https://developer.apple.com/xcode/features/)
* Windows: recommend using [MinGW](http://www.mingw.org/)

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF` Previous versions use i/o from text files.  The current state uses i/o
from the simulator.
