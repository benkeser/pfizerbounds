
## Motivation

Pfizer released their [protocol for the Phase III
trial](https://www.documentcloud.org/documents/7212814-C4591001-Clinical-Protocol.html)
of their preventive COVID-19 vaccine on September 18, 2020. The analysis
uses a Bayesian procedure to estimate P(VE \> 30 | data), with 4 planned
interim analyses + a final analysis. The success criteria are defined as
having an estimate of P(VE \> 30 | data) \> 0.995 at an interim analysis
or estimate of P(VE \> 30 | data) \> 0.986 at the final analysis. A
relatively disperse prior for VE is selected and the protocol states
that the type I error of the design is controlled at level 0.021.

This analysis plan is unique in several ways; both the number of interim
analyses and the Bayesian approach make the approach challenging to
compare to e.g., [Moderna’s
protocol](https://www.modernatx.com/sites/default/files/mRNA-1273-P301-Protocol.pdf),
which uses 2 interim analyses and O’Brien-Fleming-type alpha spending to
control type 1 error at level 0.025.

I wanted to answer the following question: if Pfizer’s success criteria
(detailed by their case splits at each interim analysis in Table 5 of
the protocol, page 103) were applied with a more standard frequentist
analysis, what are the implied success criteria and how does it compare
to more standard forms of sequential monitoring, like O’Brien Fleming.

## Results

The figure below shows the stopping boundaries on a Z-statistic scale at
each planned analysis. We see that overall **the Pfizer trial is largely
similar to a trial that uses the same number of interim analyses with an
exact binomial test with Pocock-style sequential monitoring**.

Notably, this style of monitoring is *more aggressive* stopping early
than O’Brien-Fleming-type monitoring, which is employed by Moderna.

![Plot of implied Pfizer bounds vs. OBF and
Pocock](README_files/figure-gfm/bounds-plot-1.png)

## Limitations

The computed O’Brien-Fleming and Pocock boundaries are directly computed
on the Z-statistic scale, but they do not account for the discreteness
of the outcome. Thus, the Z-statistic implied by an actual test with
this discrete outcome that would control type I error appropriately
under these spending functions may look a little different than what is
computed here. However, note that this is the Z-statistics used to
compute the Pfizer bounds *was* the one implied by the discrete test.

## Analysis details

For simplicity, I focus on a test based on an exact binomial test and
what that test would return given the case splits that define success
for Pfizer. Curious readers may consult [Dragalin & Fedorov
(2006)](https://www.tandfonline.com/doi/full/10.1080/10543400600721554)
for details (or to check my math\!). Here’s the function I wrote to
compute a p-value for this test given the inputs of a vaccine trial.

``` r
#' Compute p-value from exact binomial test
#' 
#' Uses the test described in section 3 of Dragalin & Fedorov (2006)
#' (doi: 10.1080/10543400600721554) which is an exact conditional
#' binomial test based on a Poisson approximation to the binomial. The
#' notation used in the function aligns with that used in that paper.
#' 
#' @param n1 number in placebo arm
#' @param n2 number in vaccine arm
#' @param null_ve null hypothesis vaccine efficacy
#' @param x1 number cases in vaccine arm
#' @param x2 number cases in placebo arm

binom_pval <- function(n1, n2, x1, x2, null_ve = 0.3){
    # randomization ratio
    h <- n1 / n2
    # null value on risk ratio scale
    R_star <- (1 - null_ve)
    # prob a case is from vax under null
    pi_star <- R_star / (h + R_star)
    # p-value from binomial CDF
    pval <- pbinom(x2, x1 + x2, pi_star)
    return(pval)
}
```

Now, I use this function to get the implied bounds on the scale of a
normally-distributed test statistic, which is the scale on which typical
stopping boundaries are computed.

``` r
# Pfizer parameters
num_analysis <- 5
num_endpt_at_analysis <- c(32, 62, 92, 120, 164)
crit_num_vax <- c(6, 15, 25, 35, 53)
num_vax <- 17600
num_plc <- 17600

# compute pval for exact binomial test and associated z-stat
pfizer_zstat_bound <- rep(NA, num_analysis)
pfizer_pval_bound <- rep(NA, num_analysis)
for(i in seq_len(num_analysis)){
    pfizer_pval_bound[i] <- binom_pval(n1 = num_vax, 
                                       n2 = num_plc, 
                                       x1 = num_endpt_at_analysis[i] - crit_num_vax[i], 
                                       x2 = crit_num_vax[i])
    pfizer_zstat_bound[i] <- -qnorm(pfizer_pval_bound[i])
}
```

Next, I use the
[`ldbounds`](https://cran.r-project.org/web/packages/ldbounds/index.html)
package to compute bounds associated with Lan-DeMets alpha spending
functions that mimic O’Brien Fleming and Pocock boundaries.

``` r
## computing z-stat boundaries based on alpha-spending
library(ldbounds)
obflem_bounds <- bounds(t = num_endpt_at_analysis / 164, iuse = 1, alpha = 0.025)
obflem_bounds_zstat <- obflem_bounds$upper.bounds
pocock_bounds <- bounds(t = num_endpt_at_analysis / 164, iuse = 2, alpha = 0.025)
pocock_bound_zstat <- pocock_bounds$upper.bounds
```

Finally, I put it all together into a plot.

``` r
## plotting the bounds
library(ggplot2)
library(wesanderson)

plot_data <- data.frame(
  `COVID-19 cases` = rep(num_endpt_at_analysis, 3),
  Bound = factor(c(rep("Pfizer", num_analysis), 
            rep("O'Brien-Fleming-type", num_analysis),
            rep("Pocock-type", num_analysis)), levels = c("O'Brien-Fleming-type", "Pocock-type", "Pfizer")),
  `Critical Z-statistic` = c(pfizer_zstat_bound, obflem_bounds_zstat, pocock_bound_zstat),
  check.names = FALSE
)

p <- ggplot(plot_data, aes(x = `COVID-19 cases`, y = `Critical Z-statistic`)) + 
        geom_point(aes(color = Bound), size = 3, alpha = 0.5) + 
        geom_line(aes(color = Bound), size = 1, alpha = 0.5) + 
        scale_color_manual(values = wes_palettes$Zissou1[c(1,3,5)]) + 
        theme_bw()
p
```
