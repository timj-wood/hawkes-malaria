# hawkes-malaria

Code for the blog post [Hawkes processes in epidemiology](https://timj-wood.github.io/posts/hawkes-processes/).

Simulates a malaria outbreak as a Hawkes process, fits the model to the case times, and estimates which cases were imported and which were local. Built on the [epihawkes](https://github.com/mrc-ide/epihawkes) R package from Unwin *et al.* (2021).

## Files

- `hawkes_malaria.R` – simulate, fit, classify, and make the figures
- `replicates.R` – repeat the fit on 20 simulated datasets
- `rep_plot.R` – plot the replicate results
- `results/replicates.csv` – output of `replicates.R`
- `figures/` – figures used in the post

## Run

```r
install.packages(c("remotes", "ggplot2"))
remotes::install_github("mrc-ide/epihawkes")

source("hawkes_malaria.R")
source("replicates.R")
source("rep_plot.R")
```

Run from the repository root.

## Reference

Unwin *et al.* (2021). Using Hawkes Processes to model imported and local malaria cases in near-elimination settings. *PLoS Comput Biol* 17(4): e1008830. https://doi.org/10.1371/journal.pcbi.1008830
