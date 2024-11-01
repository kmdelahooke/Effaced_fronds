This repository contains the code required to reproduce the analyses of:

<em>Delahooke, K.M., Liu, A.G., Stephenson, N.P., and Mitchell, E.G. 2024. Effaced Fronds: hidden patterns of early animal mortality? (in prep/unpublished PhD chapter)</em>


"surface_data_processing.R" takes a folder of .csv surface datasets and window coordinates, partitions taxa into 'well-preserved' and 'effaced' specimens and creates spatstat point patterns for future analysis.
Use instead "surface_data_read.R" when using the pre-partitioned data available on this repository.

Possible causes of effacement are tested in the following scripts:
- Post-fossilisation processes: 'erosion.R'
- Syn-fossilisation processes: 'taphonomy.R'
- Pre-burial mortality: 'living_dead.R'

The scripts:
- 'raster_processing_functions.R'
- 'spatial_analysis_functions.R'
- 'taxon_list.R'
  
are dependencies of the above scripts
