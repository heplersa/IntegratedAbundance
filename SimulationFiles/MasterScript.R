######################################################################################################################
####################### Master File to Implement Simulation for Paper ################################################
######################################################################################################################

#### generate the H=100 simulated data sets
#### outputs a file called SimulatedData.Rda
source('SimulateData.R')

#### runs the MCMC algorithm for all 100 simulated data sets for the joint outcome with yearly survey data case
#### uses a 4-iteration loop with 25 cores so within each loop iteration 25 runs are done in parallel
##### outputs 4 files - JointYearly.1.Rda, JointYearly.2.Rda, JointYearly.3.Rda, JointYearly.4.Rda with each 
##### saving the mean, standard deviations, and the .025 and .975 quantiles for the quantities of interests for the 25 data sets in each
source('JointYearly.R')

#### repeat for the joint outcome sparse survey data case
source('JointSparse.R')

#### repeat for single outcome yearly survey data
source('SingleYearly.R')

#### repeat for single outcome sparse survey data
#### repeat for single outcome yearly survey data
source('SingleSparse.R')

#### Summarize the results
#### this script computes the coverage probabilities and
#### errors rates presented in the paper
#### this also recreates the figures from the paper
source('PlotResults.R')