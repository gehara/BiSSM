# BiSSM
R code to simulate a Bidimentional Stepping Stone model

The main script showing an example is BiSSM_script.R

Clone the repository to your computer and source BiSSM_script.R in your R console or open it Rstudio. 
This example will run 10 simulations for two models: IBD and Island

The example dataset that this is based on comes from 7 samples of the species **Piculus aurulentus** from the Atlantic Forest, Brazil.

IBD: In this model migrations are scaled by geographic distance and dispersal rate per generations.
Island: In this model migrations are sampled from a prior.

The script is just an example and will run only 100 simulations per model plus 10 simulations as testing data. It also runs a neural network for only 10 epochs, so results are meaningless. To have meaningfull results you should change the script and run at least 10,000 simulations and train the neural network for at least 1000 epochs.   
