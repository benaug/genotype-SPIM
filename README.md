# genotype-SPIM

Nimble sampler for genotype SPIM from "Spatial proximity moderates genotype uncertainty in genetic tagging studies"

https://www.pnas.org/content/117/30/17903.short

Currently, only the Poisson and Negative Binomial observation models are supported without any occasion-level or behavioral response effects. See test scripts.

The negative binomial observation model will require "better data" in order to estimate the overdispersion parameter, i.e., more genotype information and/or less home range overlap. This sampler demonstrates why I've set up the custom ID update the way I have. Using the Poisson observation model, genoSPIM can be written in BUGS code without the custom ID update by specifying the *independent* sample-level likelihoods. However, once you switch to any other observation model, you cannot do this. The Metropolis-Hastings update will work with any observation model (assuming you switch in the correct likelihood in the custom update).

