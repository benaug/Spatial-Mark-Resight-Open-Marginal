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


Currently, I assume 1) the mark status does not change within a primary occasion (no interspersed marking and resighting), and 
2) the mark status of all marked individuals is known in all years. This likely limits us to marks that are telemetry (GPS) collars.
Marked individuals can carry their marks across years, but I currently assume mark loss/censoring is uninformative. 

The model can be modified for some unknown mark status scenarios. For example, a telemetry collar can die, so you wouldn't know that
individual is still in the population, but you may resight it and know it must be one of the marked individuals. So individuals with
dead collars would produce all "marked with  no ID" sample types.

Code is now set up so that each year can have marking effort, sighting effort, both, or neither. The only requirement is that
there is at least one marking and sighting session. Scenario set up in test script has been evaluated with 100 simulated data sets. 
With little marking effort mixing is worse. 25K iterations were sufficient and it took 2 hours on average to run on 2.2 GHz cpu.
