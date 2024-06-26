library(tidyverse)
library(mashr)

set.seed(1)

load("results/tox_de.mash_inputs.RData")

mdata = mash_set_data(as.matrix(betas), as.matrix(ses))

U.c = cov_canonical(mdata)

V.em = mash_estimate_corr_em(data, U.c, details = TRUE)
mdata = mash_update_data(mdata, V=V.em$V)

m.c = mash(mdata, U.c)

lfsr <- get_lfsr(m.c)
pm <- get_pm(m.c)
psd <- get_psd(m.c)
pihat <- get_estimated_pi(m.c)
save(lfsr, pm, psd, pihat, file="results/tox_de.mash_outputs.RData")