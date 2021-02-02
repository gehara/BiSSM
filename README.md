# BiSSM
R code to simulate a Bidimentional Stepping Stone model and calculate the number of segregating sites, nucleotide diversity, Tajima's D and pairwise Fst between populations.

## Dependencies
 To run this code you need to have the following packages installed.
 > install.packages("devtools")
 
 > devtools::install_github("gehara/PipeMaster")" 
 
 > install.packages("reshape")
 
 > install.packages("geosphere")
 
 > install.packages("keras") 
 
 *Keras can be tricky to install. See the keras webpage (https://keras.rstudio.com/) for additional information on how to install the package.*
 

The main script showing an example is BiSSM_script.R

To try it out clone the repository to your computer and *source* BiSSM_script.R in your R console or open it in Rstudio and run.

This example will run 100 simulations for two models: IBD and Island

The example dataset that this is based on 7 samples of the species *Piculus aurulentus* from the Atlantic Forest, Brazil.

IBD: migrations (Nm) are scaled by geographic distance and dispersal rate per generations, which is sampled from a uniform prior.

Island: migrations (Nm) are sampled from a uniform prior; dispersal = Nm/distance.

The two necessary inputs are North.txt and N_pop_assignment.txt (exemples included)

It is also necessary to specify other parameters (more information inside BiSSM_script.R)

The file N_observed_3371_loci.txt is the observed number of segsites, Pi, and pairwise Fst. 

The script is just an example and will run only 100 simulations per model plus 10 simulations as testing data. It also runs a neural network for only 10 epochs, so results are meaningless. To have meaningful results you should change the script and run at least 1,000 simulations and train the neural network for at least 1000 epochs.   

Several outputs will be generated.

1) *_N: simulated training data for each model
2) test_*_N: simulated test data for each model
3) PCA.pdf: prior predictive checks, Pdfs of PCA plots and model fit.
4) *IBD.pdf: estimated vs true values for each parameter of the model.
5) *Island.pdf: estimated vs true values for each parameter of the model.
6) *_result.txt: parameter estimates and error for each model