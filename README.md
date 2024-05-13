## Estimating Underreporting on HIV Cases in Kenya Using the Bayesian Approach
## Project Description
This project focuses on estimating underreporting rates for HIV cases in Kenya using a Bayesian approach. The methodology involves simulating data and analyzing real-world data to infer the true number of HIV cases based on reported data and certain assumptions.

## Technologies Used
R programming language
Statistical modeling and analysis techniques
Markov Chain Monte Carlo (MCMC) algorithms
Bayesian inference
Simulation Part
Simulated Data Generation
Simulated data is generated using a Poisson distribution with parameters lambda and P.

## Likelihood Calculation (L)
The likelihood function is defined to compute the likelihood of observing the simulated data given parameters P and lambda.

## Full Probability Calculation (FULLP)
A function to compute the full probability based on data, P, lambda, and additional parameters a and b.

## Metropolis-Hastings Algorithm (MHA)
Implementation of the Metropolis-Hastings algorithm to calculate acceptance ratios for MCMC.

## Markov Chain Monte Carlo (MCMC) Implementation (TRIAL)
A function implementing MCMC using the Metropolis-Hastings method to sample values of lambda and P.

## Real Data Analysis
Data Import and Cleaning
Real data related to HIV cases in Kenya is imported and cleaned for analysis.

## Likelihood Function and Optimization
The likelihood function is defined for real data, and optimization techniques are applied to estimate parameters P and lambda.

## Installation Instructions
Clone the repository to your local machine.
Ensure you have R installed along with necessary packages like readr for data import.
Run the provided R scripts in the specified order to simulate data, analyze real data, and estimate underreporting rates.
## Usage
Modify the parameters and functions as needed for different datasets or analysis scenarios.
Run the scripts in an R environment or IDE to execute the simulations and analyses.
## Results and Insights
The project provides insights into estimating underreporting rates for HIV cases based on simulated and real data analysis.
Visualizations and statistical summaries are generated to interpret the results effectively.

