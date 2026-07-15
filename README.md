# Spatial-Mark-Resight-Open-Marginal
Jolly-Seber spatial mark-resight samplers marginalizing out latent individual IDs and integrated telemetry survival 

This model is an SMR version of the Jolly-Seber model here (JS-SCR-Dcov):

https://github.com/benaug/Jolly-Seber-N-Prior-DA

It also has allows integrated telemetry survival data, see here for an SCR version:

https://github.com/benaug/Spatial-IPM-Telemetry

The SMR model comes from the marginal SMR models here:

https://github.com/benaug/Spatial-Mark-Resight-Marginal

This model considers sample types: marked with ID, marked with no ID, unmarked, unknown marked status. To speed up computation, it uses
the approach of Herliansyah et al. (2024, section 4.3) in the custom N/z and activity center updates.

https://link.springer.com/article/10.1007/s13253-023-00598-3

For sighting observation models other than Poisson, we can no longer marginalize out individual ID and must sample them.
See here for conditional versions:

https://github.com/benaug/Spatial-Mark-Resight-Open-Conditional

Currently, I assume 1) the mark status does not change within a primary occasion (no interspersed marking and resighting), and 
2) the mark status of all marked individuals is known in all years. This likely limits us to marks that are telemetry (GPS) collars.
Marked individuals can carry their marks across years, but I currently assume mark loss/censoring is uninformative. 

The model can be modified for some unknown mark status scenarios. For example, a telemetry collar can die, so you wouldn't know that
individual is still in the population, but you may resight it and know it must be one of the marked individuals. So individuals with
dead collars would produce all "marked with  no ID" sample types.

Code is now set up so that each year can have marking effort, sighting effort, both, or neither. The only requirement is that
there is at least one marking and sighting session. Scenario set up in test script has been evaluated with 100 simulated data sets. 
With little marking effort mixing is worse. 25K iterations were sufficient and it took 2 hours on average to run on 2.2 GHz cpu.

6/3/26: Added a version with interspersed marking and sighting within years. It is increasingly slower than the model with no
interspersion as you add occasions within years.

6/26/26: Added versions with RSF-based activity center movement between primary sessions. These test scripts are set up with abundant
marked individuals with telemetry data, which is likely required for acceptable mixing.

7/15/26: Added fixed-activity-center model versions with a gamma random effect at the primary-session-by-detector level in 
the sighting process. This allows extra-Poisson variation among detectors within each primary session. The random effects are 
analytically integrated out rather than explicitly sampled. The implementation retains the component Poisson likelihoods and appends a likelihood-ratio correction that replaces the Poisson distribution of the detector-session total with its corresponding negative-binomial distribution.

So far, versions are available both with and without an explicit occasion dimension with fixed activity centers and without an explicit
occasion dimension with mobile activity centers. In the occasion-level version, the current
implementation assumes that the same detector effect is shared across all occasions within a primary session. An alternative
formulation could assign independent random effects at the primary session by detector by occasion level.

The marginalized model runs somewhat more slowly than the base Poisson model, but should generally be more
computationally efficient than explicitly sampling every detector-level random effect.
