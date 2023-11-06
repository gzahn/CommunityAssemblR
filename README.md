# CommunityAssemblR

<img src="https://github.com/gzahn/CommunityAssemblR/blob/main/media/hex_transparent.png" alt="drawing" width="200"/>



### R functions for building mock species abundance tables and manipulating them

This is a modular set of tools for building and altering community matrices. The standard community matrix gives abundance counts for N taxa (columns) in N samples (rows). There are also tools for simulating transplantations from a donor community to a resident community under various theoretical situations. These functions are intended to be modular and chainable so that a unique community can be progressively built that reflects a variety of selective pressures. This community can then be used to test the detection accuracy of predictive assembly models.

This package is under active development. *caveat emptor*

Eventual goal: to simulate ecological communities shaped by various assembly pressures for testing model detection of 

	- Co-occurrence relationships
	- Priority effects
	- Competitive exclusion
 	- Antagonism
  	- Facilitation
  	- Niche pre-emption
  	- Keystone taxa
  	- Abiotic filtering
  	- Dispersal
  	- 
	- etc.

### Installation

```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("gzahn/CommunityAssemblR")
```

### Build Info

	- R version 4.3.1 (2023-06-16)
	- Platform: x86_64-pc-linux-gnu (64-bit)
	- Running under: Ubuntu 22.04.3 LTS

### Example workflow (so far):


[**Tutorial**](https://github.com/gzahn/CommunityAssemblR/blob/main/tutorial.pdf)


**Simulate microbiome transplantation where the dominant resident taxa are antagonistic to 50% of the new donor taxa**
```
# # make a resident community matrix (even)
even <- build_even_community(n.taxa = 100,n.samples = 44,n.reads = 3000, taxa.sd = 30) 

# # add some network hub taxa relationships and randomly increase some taxa abundances to simulate (known) variation
recipient <- link_taxa_abundances(even,n.taxa = 35, relationship = 'hub',link.scale = 3,n.links = 3) %>% 
  increase_abundances(prop = .2,increase.scale = 3,margin = "taxa")

# # build a donor community based on even resident community (to simulate a community for transplantation into the resident community)
donor <- build_donor_community(resident.comm = even, n.transplant.taxa = 30,overlap = .75)

# # simulate transplantation with resident community antagonism as the primary factor for success of novel taxa persistence
final <- transplant_w_antagonism(recipient, donor, antag.ubiq = .5, antag.strength = 10, antag.abundant = TRUE, transplant.only = TRUE)
```


**Simulate microbiome transplantation where all the new donor taxa encounter steep environmental filtration along a gradient**
```
recipient <- even()
donor <- build_donor_community(even,n.transplant.taxa = 30,overlap = .5)

final <- transplant_w_antagonism(recipient,donor,antag.strength = 0) # stochastic transplantation
final_w_filtration <- filter_taxa_along_gradient(final,prop = 1,transplant.only = TRUE,groups = 2,gradient.strength = 3)

# compare community structure with and without the environmental filtration
final %>% t() %>% heatmap(Colv = NA,Rowv = NA)
final_w_filtration %>% t() %>% heatmap(Colv = NA,Rowv = NA)
```

### Citing:

Geoffrey Zahn. (2023). gzahn/CommunityAssemblR: 0.1.1 (0.1.1) [Computer software]. Zenodo. https://doi.org/10.5281/ZENODO.10073694

[![DOI](https://zenodo.org/badge/714101784.svg)](https://zenodo.org/doi/10.5281/zenodo.10073694)




