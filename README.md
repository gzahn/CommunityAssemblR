# CommunityAssemblR

<img src="https://github.com/gzahn/CommunityAssemblR/blob/main/media/hex_transparent.png" alt="drawing" width="200"/>



### R functions for building mock species abundance tables and manipulating them

This is a modular set of tools for building and altering community matrices. The standard community matrix features abundance counts for N taxa (columns) in N samples (rows). There are also tools for simulating transplantations from a donor community to a resident community under various theoretical situations.

This package is under active development.

Eventual goal: to simulate various ecological assembly theories for testing model detection of 

	- Network relationships
	- Priority effects
	- Competitive exclusion
 	- Antagonism
  	- Facilitation
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

### Example workflow (sa far):

```
# # make a resident community matrix (even)
even <- build_even_community(n.taxa = 100,n.samples = 44,n.reads = 3000, taxa.sd = 30) 
# # add some network hub taxa relationships and randomly increase some taxa abundances to simulate (known) variation
recipient <- link_taxa_abundances(even,n.taxa = 35, relationship = 'hub',link.scale = 3,n.links = 3) %>% 
  increase_abundances(prop = .2,increase.scale = 3,margin = "taxa")
# # build a donor community based on even resident community (to simulate a community for transplantation into the resident community)
donor <- build_donor_community(resident.comm = even, n.transplant.taxa = 30,overlap = .75)
# # perform transplantation with resident antagonism as the primary factor for success of novel taxa transplantion
final <- transplant_w_antagonism(recipient, donor, antag.ubiq = .5, antag.strength = 10, antag.abundant = TRUE)
```
