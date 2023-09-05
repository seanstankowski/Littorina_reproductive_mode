This directory contains two jupyter notebooks and many files that correspond with simulations conducted in our study. One that was used to simulate genealogies under a four-taxon model using the program MSprime, and a second one used to analyse the distribution of topology weights in the ternary framework. Below is a slight more detailed description of the purpose of each script, but both notebooks are well annotated within.

Ternary_Simulations.v2.ipynb:
This notebook controls the coalescent simulator MS prime, allowing the user to simulate large numbers (10s to hundreds of thousands) of trees under a four-population framework. Various parameters, including the duration and timing of splits, gene flow, and effective population size can be manipulated. We used these to conduct the simulations shown in figs S9 to S19, and table S4. 

ExploringTheTernaryForrest.ipynb:
This notebook contains coded needed to produce the figures shown in fig. S10 to S19, as well at the quantitative estimates of asymmetry shown in table S4. 

.tree:
These files contain the trees produced during each simulation. The first part of the name of each file corresponds to the simulation codes in table S4. 

weights.csv: 
These files contain results of topology weighting for each simulation (topology weighting performed in Twisst). The first part of the name of each file corresponds to the simulation codes in table S4. 
