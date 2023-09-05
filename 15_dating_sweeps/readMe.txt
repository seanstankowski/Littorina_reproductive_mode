The scripts in this directory were used to estimate the ages of selective sweeps at live-bearing alleles, as reported in the main text and fig. S32.

Rough_dating_of_sweeps.ipynb:
This Jupyter notebook contains the code used to date the sweeps. 

pi_pvt:
This file contains all of the data needed for the sweep estimation and is called by Rough_dating_of_sweeps.ipynb. 
Contig_id	The ID of the assmbly contig
popgen_ID	An aribtraty ID for each region where tr ==1 
start	The start of the region
end	The end of the region
mid	The midpoint of the region
region_length	The total length of the rgion
sites	Total number of sequenced sites in the region
sites_shared	The number of shared alleles between egg-layers and live bearers that were removed from the analysis before pi_ private was calculated 
sites_corrected	The number of sites with the with the number of shared sites put back in to correct the denominator for the pi calculation
pi_Egglayer_pvt	pi egg-layer at private alleles before correcting sequence length
pi_Brooder_pvt	pi live bearer at private alleles before correcting sequence length
pi_Brooder_pvt_corrtected	pi live bearer after correction 
dxy_Egglayer_Brooder	dxy between between egg-layers and live bearers
pi_Egglayer	pi egg-layer all sites
pi_Brooder	pi-live bearer all sites
T	time since the sweep of the live_bearing alleles