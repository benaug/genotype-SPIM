# genotype-SPIM

Nimble sampler for genotype SPIM from "Spatial proximity moderates genotype uncertainty in genetic tagging studies"

https://www.pnas.org/content/117/30/17903.short

Currently, only the Poisson and negative binomial observation models are supported without any occasion-level or behavioral response effects. For the Poisson observation model, there is one version with 1 set of genotyping error rates applied to all samples (the only model for the negative binomial observation model) and a second version that allows categorical sample type covariates (e.g., low and high quality samples). See test scripts.

The negative binomial observation model will require "better data" in order to estimate the overdispersion parameter, i.e., more genotype information and/or less home range overlap. This sampler demonstrates why I've set up the custom ID update the way I have. Using the Poisson observation model, genoSPIM can be written in BUGS code without the custom ID update by specifying the *independent* sample-level likelihoods. However, once you switch to any other observation model, you cannot do this. The Metropolis-Hastings ID update will work with any observation model (assuming you switch in the correct likelihood in the custom update).

DISCLAIMER: There was been a bug before 2/17/2022 in the G.true update. It has been fixed and the sampler has been tested more extensively to make sure the correction fixed the issue.