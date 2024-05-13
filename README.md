# Estimating Underreporting on HIV Cases in Kenya Using the Bayesian Approach
## Introduction
This project focuses on estimating the true proportion of HIV/AIDS cases in Kenya using statistical models. The goal is to address challenges related to underreporting and improve the accuracy of population parameter estimates.

## Problem Statement
Underreporting of HIV/AIDS cases undermines the accuracy of population estimation, impacting resource allocation, targeted interventions, and progress monitoring.

## Objectives
Estimate the true proportion of HIV/AIDS cases in Kenya.
Develop and implement statistical models to improve estimation accuracy.
Determine the true population parameters using Bayesian approaches.
## Methodology
Data Source
The primary dataset is sourced from the National Syndemic Disease Control Council (NSDCC), containing HIV intervention data.

Statistical Models
Utilizing Poisson and binomial models to account for underreporting.
Bayesian approach with Markov Chain Monte Carlo (MCMC) using the Metropolis-Hastings method.
Usage
## Requirements
R programming language
Required R packages: readr

## Steps
Import the real data using the readr package.
Clean the data and divide it into yearly datasets (e.g., 2014, 2015, 2017).
Run the statistical models (e.g., L1, L2, L3) to estimate true population parameters.
Analyze results, including descriptive statistics and underreporting percentages.

## Results
Significant underreporting percentages observed across multiple years.
Consistent bias in reported HIV cases, indicating a need for adjustment in estimation methods.

## Conclusion and Recommendations
Adjustments for underreporting are crucial for accurate population estimation.
Consider integrating covariates into models and exploring alternative approaches for count data.


