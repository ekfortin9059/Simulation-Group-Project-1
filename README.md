# Simulation
# Group Project 1
## BoxCar Ride-Sharing Simulation

**Authors:**
1. Erin Fortin 
2. Alice Locker
3. Colin McDonagh

**Date:** Tuesday, 10 March 2026

---

## Overview of Main Files
1. `run_simulation.ipynb`: Main document to run the various versions of the simulation of BoxCar's ride-sharing platform. Contains the main simulation code and performance metrics calculations for each of the model versions.  
2. `Statistical_Tests.ipynb`: Checking Distributional Assumptions of the Model by Comparing to Real Observed Data (`drivers.xlsx` and `riders.xlsx`)
3. `Simulation(Dropoff_Locations).R`, `Simulation(Initial_Locations).R`, and `Simulation(Pickup_Locations).R`: These R files contains code to generate samples from the posterior distribution for the mean vector and covariance matrix parameters corresponding to the Bayesian model that assumes the given locations are distributed from a multivariate normal distribution (where the mean vector and covariance matrix are to be estimated from the data). 
4. `uniform_testing.ipynb`: Analysis of the validity of BoxCar's assumptions about certain arrivals and locations following a Uniform distribution

## Supplementary File
1. `requirements.txt`: contains the packages used to complete this project. The user of these documents should make sure these packages are installed before running the code. 
