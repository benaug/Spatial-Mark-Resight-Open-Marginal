# Spatial Mark-Resight Open Marginal

Open-population generalized spatial mark-resight models using N-prior data augmentation, marginalization over latent individual identities, 
optional integrated telemetry data, and optional integrated trap random effects.

## Overview

This repository extends the marginalized spatial mark-resight (SMR) models in [Spatial-Mark-Resight-Marginal](https://github.com/benaug/Spatial-Mark-Resight-Marginal) to open populations using the Jolly-Seber N-prior data augmentation framework in [Jolly-Seber-N-Prior-DA](https://github.com/benaug/Jolly-Seber-N-Prior-DA).
Here, I only consider the generalized SMR model ([Whittington et al. 2018](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2664.12954)) 
which was introduced for closed populations to correct for spatially non-random sampling to deploy marks. This approach has also been termed
"Spatial Capture Mark-Resight" by [Efford and Hunter (2018)](https://onlinelibrary.wiley.com/doi/full/10.1111/biom.12766),
making the capture process explicit, which is probably better terminology. More recently a second modeling approach has
been classified as generalized SMR ([Whittington et al. 2025](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecs2.70246)). 
This approach discards information in order to avoid modeling the marking process. In my opinion, this second approach does not meet the definition of
generalized SMR stated in Whittington et al. (2018) that ties it to modeling the marking process. Regardless, here, I mean generalized SMR in the 
Whittington et al. (2018) context.

Basic MR models assume that the capture probability of marked and unmarked individuals are the same. In the context of SMR, 
this assumption requires that the spatial distribution of marked and unmarked animals across the state space are the same, which 
will not be true if, for example, individuals are targeted for marking disproportionately inside or around the resighting detector
array. By modeling the marking process, we are modeling the process by which the marked and unmarked individuals come to have have different
spatial distributions, which removes the bias induced when ignoring this process. A further implication of spatially non-random sampling is
that covariates for density-habitat relationships will be estimated with bias. Moving to open populations, modeling the marking process
is required to use telemetry survival states to inform survival parameters without introducing bias. Finally, the marking process contains 
valuable information about recruitment that is not present in the sighting data. For all these reasons, I do not consider open population
SMR without the marking process data here.

Extending SMR to an open-population model also leads to a spatial integrated population model 
(IPM, [Chandler and Clark 2014](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12153)) that combines
marked detections, unmarked detections, partially marked detections (e.g., marked with no ID samples), telemetry informed 
within and between year movement, telemetry survival data, or other data streams. Chandler and Clark (2014) introduced the spatial
IPM with an example where primary sessions may contain completely marked or completely unmarked detection data--SCR and unmarked SCR primary
occasions. In the open-population models I consider here, I allow any combination of marking and sighting observation processes across 
primary occasions, and consider that marks can be carried across primary occasions. One could do marking and resighting in only a subset
of primary occasions and sighting only in the remainder, with the possibility of primary occasions with no observation process. Another option
is that marking and sighting do not have to be in the same primary session, allowing demographic change between marking and sighting, which is
not allowed in the closed population approach (Whittington et al. 2018).

There an many modeling approaches one might take for open-population SMR which will differ in MCMC efficiency and modeling flexibility. Here,
I use the model in the Jolly-Seber-N-Prior-DA repository and the marginal SMR observation model for specific reasons.
First, the Jolly-Seber model of Chandler and Clark (2014) has the desirable feature that the recruitment model is per capita as a function
of realized abundance in the previous primary occasion. The per capita formulation makes a fixed recruitment parameter through time more
ecologically realistic, and a single parameter is better identified and more efficiently estimated than using one recruitment parameter per primary
occasion. Alternatively, a trend, other continuous relationship through time, or random effects are more ecologically realistic with a per 
capita parameterization. One could opt for per capita recruitment as a function of expected abundance, but the realized parameterization
is arguably more ecologically realistic, and implies more demographic stochasticity.

While the Chandler and Clark (2014) model is ecologically realistic, it can be very slow to run due to the data augmentation design, 
particularly as the level of required population augmentation, M, increases. Further, the entry distributions are binomial approximations
of Poisson distributions that are achieved only as M goes to infinity. Therefore, I opt to use the exact Poisson specification of this model
in the Jolly-Seber-N-Prior-DA repository that also uses a more efficient MCMC algorithm to update the latent population entry and exit dates of each
individual. 

Second, I chose the marginal observation model because it mixes better than the conditional version where individual identities of latent
identity samples must be updated in the MCMC algorithm. This is especially important for open populations where we are introducing many more
latent variables that must be traversed along with the latent individual identities. One drawback of the typical marginal observation model
approach used in BUGS-based software, though, is that updates of individual inclusion indicators and activity centers are very slow because
the total expected detection rate for unidentified mark types are resummed across all individuals when updating these parameters for each
single individual. We avoid this by using the approach of [Herliansyah et al. (2024), section 4.3](https://link.springer.com/article/10.1007/s13253-023-00598-3)
inside the custom updates for individual inclusion indicators and activity centers. A second drawback is that we can only use the Poisson
distribution for the sighting model.

## Related repositories

The open-population formulation, including the population model, N-prior data augmentation, survival and recruitment processes, 
latent population-state histories, mobile activity-center models, and associated MCMC updates, is described in the Jolly-Seber repository.

- [Jolly-Seber-N-Prior-DA](https://github.com/benaug/Jolly-Seber-N-Prior-DA)

The SMR observation model used here that allows for unknown marked status observations is from the closed-population marginalized SMR 
repository. 

- [Spatial-Mark-Resight-Marginal](https://github.com/benaug/Spatial-Mark-Resight-Marginal)

For a Poisson sighting model, latent individual identity can be marginalized analytically for the unidentified observation
types. Individual identity therefore does not need to be sampled during MCMC, which improves mixing, especially
for open-population models. For sighting observation models for which individual identity cannot be marginalized analytically, e.g. 
negative binomial, conditional open-population SMR models that explicitly sample latent individual IDs are available here:

- [Spatial-Mark-Resight-Open-Conditional](https://github.com/benaug/Spatial-Mark-Resight-Open-Conditional)

An SCR implementation integrating telemetry survival and location data that allows mixing and matching 2 SCR observation models across
primary periods (another spatial IPM) is available here:

- [Spatial-IPM-Telemetry](https://github.com/benaug/Spatial-IPM-Telemetry)

## Spatial mark-resight observation model

The observation model contains separate marking and sighting processes. Individual identity is known in the marking data. During 
sighting occasions, individual identity and marked status may be fully observed, partially observed, or unknown.

### Marking process

For individual $i$, primary occasion $g$, and marking detector $j$, detection probability is modeled using a half-normal detection function,

$$
p_{i,g,j} = p_{0,g,j}\exp\left(-\frac{\|\mathbf{s}_{i,g}-\mathbf{x}_{g,j}^{M}\|^2}{2\sigma^2}\right),
$$

where $'\mathbf{s}_{i,g}'$ is the activity center of individual $i$, $'\mathbf{x}_{g,j}^{M}'$ is the location of marking detector $j$,
$p_{0,g,j}$ is baseline marking detection probability, and $\sigma$ is the spatial scale parameter. Conditional on being alive,

$$
Y_{i,g,j}^{M} \sim\mathrm{Binomial}\left(K_{g,j}^{M},p_{i,g,j}\right),
$$

where $K_{g,j}^{M}$ is marking effort. Model files are set up to share the same $\sigma$ between the marking and sighting observation processes, but this can be modified.

### Sighting process

The sighting process is generated in two stages. First, each individual generates a latent true number of sightings at each detector.
Second, those sightings are classified according to whether individual identity and marked status are observed.

For individual $i$, primary occasion $g$, and sighting detector $j$, let

$$
\lambda_{i,g,j} = z_i^{\mathrm{super}}z_{i,g}\lambda_{0,g}\exp\left(-\frac{\|\mathbf{s}_{i,g}-\mathbf{x}_{g,j}^{S}\|^2}{2\sigma^2}
\right),
$$

where $z_i^{\mathrm{super}}$ indicates membership in the superpopulation, $z_{i,g}$ indicates that the individual is alive during primary
occasion $g$, $\lambda_{0,g}$ is the baseline sighting rate, $\mathbf{s}_{i,g}$ is the activity center, $\mathbf{x}_{g,j}^{S}$ is the
detector location, and $\sigma$ is the spatial scale parameter.

The latent true number of sightings is

$$
Y_{i,g,j}^{\mathrm{true}}\sim\mathrm{Poisson}\left(K_{g,j}^{S}\lambda_{i,g,j}\right),
$$

where $K_{g,j}^{S}$ is sighting effort.

Let $m_{i,g}$ denote the mark state of individual $i$ during primary occasion $g$, with $m_{i,g}=1$ if the individual is marked and $m_{i,g}=0$ otherwise. This notation applies to the standard model versions, in which mark state is constant within a primary occasion.

In the `Interspersed` versions, mark state is instead indexed by sighting occasion,

$$
m_{i,g,k}=\begin{cases}1, & \text{if individual } i \text{ is marked at sighting occasion } k \text{ of primary occasion } g,\\
0, & \text{otherwise}.\end{cases}
$$

Thus, an individual can have $m_{i,g,k}=0$ for sighting occasions occurring before it is marked and $m_{i,g,k}=1$ for subsequent sighting
occasions within the same primary occasion, or it may lose the mark, be remarked, etc. The standard model is the special case in which

$$
m_{i,g,k}=m_{i,g}
$$

for all sighting occasions $k$ within primary occasion $g$.

For clarity, the observation-model and marginalization equations below are written using the standard, non-interspersed notation $m_{i,g}$. 
The same thinning and marginalization apply in the interspersed model separately for each sighting occasion $k$, replacing $m_{i,g}$ with
$m_{i,g,k}$ and retaining the sighting-occasion dimension in the observed counts.

For marked individuals, define

$$
\boldsymbol{\theta}^{M}=\left(\theta_{\mathrm{ID}}^{M},\theta_{\mathrm{noID}}^{M},\theta_{\mathrm{unk}}^{M}\right),
$$

where the three elements are the probabilities that a sighting is recorded as marked with ID, marked with no ID, or unknown marked status.

For unmarked individuals, define

$$
\boldsymbol{\theta}^{U}=\left(0,\theta_{\mathrm{um}}^{U},\theta_{\mathrm{unk}}^{U}\right),$$

where

$$
\theta_{\mathrm{unk}}^{U}=1-\theta_{\mathrm{um}}^{U}.
$$

The first element is zero because an unmarked individual cannot generate an individually-identifiable observation. Conditional on the latent true count,
sightings are thinned into the observed sample types using a multinomial distribution. If individual $i$ is marked,

$$
\left(Y_{i,g,j}^{\mathrm{ID}},Y_{i,g,j}^{\mathrm{noID}},Y_{i,g,j}^{\mathrm{unk},M}\right)
\mid
Y_{i,g,j}^{\mathrm{true}},m_{i,g}=1\sim\mathrm{Multinomial}\left(Y_{i,g,j}^{\mathrm{true}},
\left(\theta_{\mathrm{ID}}^{M},\theta_{\mathrm{noID}}^{M},\theta_{\mathrm{unk}}^{M}\right)\right).
$$

If individual $i$ is unmarked,

$$
\left(Y_{i,g,j}^{\mathrm{ID}},Y_{i,g,j}^{\mathrm{um}},Y_{i,g,j}^{\mathrm{unk},U}\right)
\mid Y_{i,g,j}^{\mathrm{true}},m_{i,g}=0\sim\mathrm{Multinomial}
\left(Y_{i,g,j}^{\mathrm{true}},\left(0,\theta_{\mathrm{um}}^{U},\theta_{\mathrm{unk}}^{U}\right)\right).
$$

Thus, every true sighting is assigned to exactly one observed sample type. This is the conditional observation model used in the
data-generating process and in the conditional SMR repository that updates latent individual IDs instead of marginalizing them out.

## Marginalization over sample type and individual identity

Because the latent true count is Poisson and the observed sample types are generated by multinomial thinning, the individual-level counts
in each sample type are independent Poisson random variables after marginalizing over $Y_{i,g,j}^{\mathrm{true}}$.

For a marked individual,

$$
Y_{i,g,j}^{\mathrm{ID}}\sim\mathrm{Poisson}\left(K_{g,j}^{S}\lambda_{i,g,j}\theta_{\mathrm{ID}}^{M}\right),
$$

$$
Y_{i,g,j}^{\mathrm{noID}}\sim\mathrm{Poisson}\left(K_{g,j}^{S}\lambda_{i,g,j}\theta_{\mathrm{noID}}^{M}\right),
$$

and

$$
Y_{i,g,j}^{\mathrm{unk},M}\sim\mathrm{Poisson}\left(K_{g,j}^{S}\lambda_{i,g,j}\theta_{\mathrm{unk}}^{M}\right).
$$

For an unmarked individual,

$$
Y_{i,g,j}^{\mathrm{um}}\sim\mathrm{Poisson}\left(K_{g,j}^{S}\lambda_{i,g,j}\theta_{\mathrm{um}}^{U}\right),
$$

and

$$
Y_{i,g,j}^{\mathrm{unk},U}\sim\mathrm{Poisson}\left(K_{g,j}^{S}\lambda_{i,g,j}\theta_{\mathrm{unk}}^{U}\right).
$$

Marked-with-ID sightings retain individual identity and remain in the likelihood at the individual level. For the unidentified observation
types, individual identity is not observed. Define the total underlying sighting intensity from marked individuals as

$$
\Lambda_{g,j}^{M}=\sum_{i=1}^{M}m_{i,g}\lambda_{i,g,j},
$$

and from unmarked individuals as

$$
\Lambda_{g,j}^{U}=\sum_{i=1}^{M}(1-m_{i,g})\lambda_{i,g,j}.
$$

For the interspersed models, these quantities are sighting-occasion specific,

$$
\Lambda_{g,j,k}^{M}=\sum_{i=1}^{M}m_{i,g,k}\lambda_{i,g,j},
$$

and

$$
\Lambda_{g,j,k}^{U}=\sum_{i=1}^{M}(1-m_{i,g,k})\lambda_{i,g,j}.
$$

Thus, the same individual can contribute to $\Lambda_{g,j,k}^{U}$ before it is marked and to $\Lambda_{g,j,k}^{M}$ after it is marked within
the same primary occasion. 

The observed marked-with-no-ID count is

$$
Y_{g,j}^{\mathrm{noID}}=\sum_{i=1}^{M}Y_{i,g,j}^{\mathrm{noID}}.
$$

By the superposition property of independent Poisson random variables,

$$
Y_{g,j}^{\mathrm{noID}}\sim\mathrm{Poisson}\left( K_{g,j}^S \theta_{\mathrm{noID}}^M \Lambda_{g,j}^{M} \right).
$$

Similarly, the observed unmarked count is

$$
Y_{g,j}^{\mathrm{um}}=\sum_{i=1}^{M}Y_{i,g,j}^{\mathrm{um}},
$$

with

$$
Y_{g,j}^{\mathrm{um}}\sim\mathrm{Poisson}\left(K_{g,j}^{S}\theta_{\mathrm{um}}^{U}\Lambda_{g,j}^{U}\right).
$$

Sightings with unknown marked status can arise from either marked or unmarked individuals,

$$
Y_{g,j}^{\mathrm{unk}}=\sum_{i=1}^{M}\left(Y_{i,g,j}^{\mathrm{unk},M}+Y_{i,g,j}^{\mathrm{unk},U}\right),
$$

giving

$$
Y_{g,j}^{\mathrm{unk}}\sim\mathrm{Poisson}\left(K_{g,j}^{S}\left(\theta_{\mathrm{unk}}^{M}\Lambda_{g,j}^{M}+
\theta_{\mathrm{unk}}^{U}\Lambda_{g,j}^{U}\right)\right).
$$

The marginalized likelihood therefore follows from two properties of the Poisson distribution:

1. **Poisson thinning:** multinomial classification of a Poisson count produces independent Poisson counts for the resulting sample types.
2. **Poisson superposition:** summing independent Poisson counts across latent individual identities produces another Poisson
count with rate equal to the sum of the individual rates.

No other count distribution has both of these properties.


## Mark states

This model requires that the mark status of all individuals is known in every primary session, therefore, it is most directly
applicable to studies using telemetry collars for marks because one can typically know when the collar stops working and therefore
providing individual identity information. In the standard model versions, mark status is assumed constant within a
primary occasion. The `Interspersed` versions instead retain mark status at the secondary-occasion level, allowing mark status to change
within a primary occasion as individuals are captured and marked.
The interspersed approach is also required to account for any deaths that happen during resighting. While more flexible, the 
interspersed version requires modeling the occasion dimension and therefore can be much slower, depending on the number of occasions.
One might choose to aggregate occasions to speed up run time at the cost of less accuracy of the mark states relative to sightings.

One feature that is not currently implemented here is that a dead GPS collar may only provide "marked with no ID" observation types. The proper
way to account for this is to include a second marked state that represents "dead GPS collar", and estimate the probability of "marked with
no ID" and "unknown marked status" observation types separately. There is an example of this in the marginal SMR repository for closed 
populations. Another approach that discards information is to treat these detection types as "unmarked", leaving no "marked with no ID" 
observation types.


## Integrated telemetry data

Telemetry information can contribute to both survival and the spatial observation model.

### Telemetry survival

Telemetry survival states provide direct information about the latent population state $z_{i,g}$. A telemetry observation can identify
an individual as known alive or known dead during a primary occasion, while missing or censored observations leave the corresponding 
state unknown. These observations are integrated directly with the Jolly-Seber latent survival histories. The current implementation
assumes telemetry censoring is uninformative.

### Telemetry locations

Telemetry locations contribute to the spatial likelihood and inform activity centers and $\sigma$. For fixed-activity-center models,
all telemetry locations for an individual are centered on its common activity center. For mobile-activity-center models, locations are
associated with the activity center for the corresponding primary occasion.

The activity-center and movement models themselves are described in [Jolly-Seber-N-Prior-DA](https://github.com/benaug/Jolly-Seber-N-Prior-DA).

## Integrated detector random effects

Model versions containing `TrapRE` allow additional heterogeneity among detector × sighting-session combinations.

Let $u_{g,j}$ be a multiplicative random effect on the sighting intensity at detector $j$ during primary occasion $g$, with

$$
u_{g,j}\sim\mathrm{Gamma}\left(\theta_{d,g},\theta_{d,g}\right),
$$

using the shape-rate parameterization. Therefore,

$$
E(u_{g,j})=1
$$

and

$$
\mathrm{Var}(u_{g,j})=\frac{1}{\theta_{d,g}}.
$$

Conditional on $u_{g,j}$, the same detector random effect multiplies the sighting intensities from all individuals and therefore affects
all sighting sample types at that detector. The random effects are integrated out rather than explicitly sampled during MCMC. This preserves
the marginalized individual-ID observation model while allowing extra detector-level variation in sighting rates. Let $T_{g,j}$ denote
the total number of sightings at detector $j$ and let $\mu_{g,j}$ denote the corresponding expected total count before applying the detector
random effect. Integrating over $u_{g,j}$ gives

$$
T_{g,j}\sim\mathrm{NegativeBinomial}\left(\theta_{d,g},\mu_{g,j}\right),
$$

under the mean-dispersion parameterization with

$$
E(T_{g,j})=\mu_{g,j}
$$

and

$$
\mathrm{Var}(T_{g,j})=\mu_{g,j}+\frac{\mu_{g,j}^2}{\theta_{d,g}}.
$$

The implementation evaluates the ordinary marginalized Poisson likelihood for the individual sample types and adds the likelihood correction
induced by integrating over the shared detector random effect. As $\theta_{d,g}$ increases, detector-level heterogeneity decreases and 
the model approaches the corresponding model without `TrapRE`.

## Model versions

The repository contains seven model versions. All use the same marginalized SMR observation model and open-population framework described above.

### Fixed activity centers

1. **JS-SMR-Dcov Generalized**

   Generalized open-population marginalized SMR model with fixed activity centers. Marking and sighting effort can vary among primary
   occasions, including primary occasions with marking only, sighting only, both methods, or no observation effort. Mark status is constant within each primary occasion, and sighting observations are summed over secondary occasions before fitting the model.

2. **JS-SMR-Dcov Generalized TrapRE**

   Version of model 1 with integrated detector random effects for the sighting process. A multiplicative gamma random effect is shared
   by all sighting types at a detector within a primary occasion and is integrated out analytically rather than sampled during MCMC.

3. **JS-SMR-Dcov Generalized Interspersed**

   Modification of model 1 that allows marking and sighting occasions to be interspersed within a primary occasion. Mark status can
   therefore change during the primary occasion as individuals are captured and marked. Sighting observations and mark states are retained
   at the secondary-occasion level so that each sighting is modeled using the individual's mark status at that time.

4. **JS-SMR-Dcov Generalized Interspersed TrapRE**

   Version of model 3 with integrated detector random effects for the sighting process.

### Mobile activity centers

The following models replace the fixed activity centers with the mobile activity-center model described in
[Jolly-Seber-N-Prior-DA](https://github.com/benaug/Jolly-Seber-N-Prior-DA).

5. **JS-SMR-Dcov-mobileAC Generalized**

   Mobile activity-center version of model 1. Individuals have an activity center for each primary occasion, with movement among
   primary occasions modeled using the Gaussian movement and resource-selection model described in the Jolly-Seber repository. Mark status
   is constant within each primary occasion and sighting observations are summed over secondary occasions.

6. **JS-SMR-Dcov-mobileAC Generalized TrapRE**

   Version of model 5 with integrated detector random effects for the sighting process.

7. **JS-SMR-Dcov-mobileAC Generalized Interspersed**

   Mobile activity-center version of model 3. Activity centers can move among primary occasions, while marking and sighting occasions can
   be interspersed within primary occasions and mark status can change as individuals are marked.