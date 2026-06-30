# GeneFlow-Mark-1
## GeneFlow: An RShiny-based metapopulation and gene flow simulation

L. W. Stearman

This repository holds the first iteration of an RShiny app to conduct metapopulation and gene flow models in a customizable watershed network. This version is planned to be superceded by a version in development, which should run considerably faster, allow simulations with far more individuals and loci, and have more options for watershed customization. 

#### Quick-start guide
The model is preconfigured to run a basic simulation with no major watershed modifications, a genome of 25 biallelic loci, and an organism which has an environmental optimum in 2nd order streams. To get started,
 - Run the app.R file
 - On the tab **Simulation Settings**, under the "Riverscape Settings: Please select an input file", use the "Browse" button to navigate to the folder /www and select one of the three watershed csv files
 - On the tab **Laboratory**, select the big *Simulate Metapopulations* button. This will start the metapopulation component of the simulation. The process typically runs in seconds, and when complete, a number of inferential graphs will appear along with a new *Simulate Gene Flow* button.
 - Select the big *Simulate Gene Flow* button. This process takes anywhere from a few minutes to 15 or more minutes, depending on your computer's speed. Once this process is complete, more inferential plots will appear.
 - Explore the graphical outputs! Under metapopulations, you'll find analyses of change through time and analyses of patterns by site. Under gene flow, you'll find analyses of pairwise Fst, a PCOA plot, and an isolation by distance plot.

## What is this thing?
GeneFlow is an RShiny app which leverages separate hybrid cohort- and individual-based metapopulation and gene flow simulations to examine demographic and evolutionary patterns among habitat patches. The simulation has an annual timestep, and it models an annual, biallelic, sexual species. The simulation is designed to work best when patches are arranged in a linear or dendritic fashion; however, it may work with reticulate or other connective patterns (currently untested). The core functions to run the simulations can be accessed and run on their own in R for more flexibility. However, the RShiny interface provides users with a simple way to control a large number of input parameters to the model. 

The GeneFlow model was originally developed to provide simulations of fragmentation and habitat instability in certain patches of river networks on headwater species. The default parameters reflect the original target species. However, depending on your processing power, you may be able to run much larger watersheds, higher numbers of individuals, more loci, or some of all of the above. 

## Model components
### Simulation Settings
This tab is where users set the parameters for the simulation. 
#### Riverscape Settings
This section allows users to make changes to parameters which affect the watershed or riverscape. It contains the following components:
##### Map input
Users can select from among three different input map files, or upload a custom map file. Map files are arranged as two-column from-to tables, where the left column indicates the donor segment, and the right column indicates the recipient segment. The row with a blank recipient segment (no downstream unit) is taken to be the downstream-most unit. Strahler Stream Orders of upstream segments are calculated using this unit. In addition to selecting a map, users can also select a maximum dispersal distance in any one year. 
##### Random Disturbance Events
This box allows users to wreak randomized havoc on habitat patches, simulating pulse disturbances. There is a slider to determine disturbance frequency, and a box where users can enter which patches they want to experience disturbances. All patches can be selected by typing "ALL" into this box. Disturbances are 100% mortality events for the patch, and occur probabilistically. If all patches are selected, only some will have any disturbances in a given iteration, unless the frequency is set to 1 (annual). 
##### Habitat Degradation Degree
Carrying capacity (*K*) is determined in part by the stream order of the segment, and in part by the species-specific settings (below). Stream size is a first-order approximation for a number of habitat variables, and many species may be specialists for certain sizes of streams or rivers. The habitat degradation box allows users to arbitrarily lower *K* for specific reaches, simulating press disturbances. Similar to pulse disturbances, there is a slider for the percent reduction in *K*, and a box to enter which patches to degrade. All patches can be selected by typing "ALL" into this box.
##### Habitat Patch Fragmentation Type
This box allows users to fragment the watershed, simulating dams or other barriers. Fragmentation types can be set to none, upstream-only (simulating dams), or total (bidirectional, simulating some culverts and other structures). Users can enter which patches have a fragment *downstream* of them in the box below.
##### Habitat patch map figure
This figure plots a rudimentary dendrogram of the patches. The algorithm attempts to find the most efficient way to display them, and they may not visually resemble an actual map of your sites at first. Stream order is shown by the size of the patch points. Fragmented, degraded, and disturbed patches are indicated by changes in line type, point color, and point shape, respectively.
##### Habitat patch connectivity figure
This figure shows a patch connectivity plot. Users can visualize which patches connect to which other patches. Notice that changing the maximum dispersal distance or applying fragmentation changes this plot. 

#### Species Settings
This section allows users to modify various parameters of the species. 
##### Population and Metapopulation Dynamics
This box contains two sliders, one which allows users to select the instantaneous rate of increase (*r*), and the other which allows users to select a proportion of individuals likely to disperse in each generation. 
##### Genome Settings
This box contains a slider for the number of loci and a slider for the power of 10 at which mutation rates should occur. Simualtions with >50 loci currently run *very* slowly.
##### Habitat Specificity
This box contains a number of user inputs to tune the value of (*K*) throughout the watershed. Users can select an optima, variance, and range of *K* values for the entire watershed. Below the input boxes, a total estiamted watershed *K* is shown with the user's selected parameter values. Note that this is best viewed as a long-term average, as the population model can and will over- and under-shoot this, mimicking actual fluctuations in population abundances over time.
##### *K* Across Stream Size Plot
This plot shows users an estimate of idealized *K* across each Strahler Stream Order.

#### Simulation Settings
This section allows users to modify other simulation parameters not covered by the sections above.
##### Select starting configuration
This toggle allows users to select whether they are running from scratch, or from a prior run. Currently, prior runs are very sensitive to input files, and users may have difficulty using prior runs. Note that currently **loading a prior run does not edit your parameters above to those used in the prior run**. Using this toggle set to *Start from Scratch* will overide any uploaded prior runs. 
##### Start from Scratch
This box covers multiple parameters needed to run a new simulation. Users can select an initial number of individuals to stock per plot, as well as plots to stock. This allows testing colonization model scenarios (i.e., headwater vs river mouth capture).
##### Start from a Prior Run
This box allows users to load the results of a prior run. The idea is that users could run a burn-in run, and then modify the watershed to test hypotheses. If prior run is loaded, the parameters of that run will be listed below this box.
##### Population Parameters
This box allows users to input the total number of generations in the simulation, and to set the intervals at which output genetic results are saved. Long runs may benefit from saving a subset of genetic runs, as these files can get very large very quickly. 

### Laboratory
This tab is where simulations are actually run, and results can be explored.
#### Simulate Metapopulations
Selecting this button will run the metapopulation portion of the simulation with the parameters you have selected on the **Simulation Settings** tab. Note that this does not automatically run the gene flow simulation. Running this portion of the simulation will generate six graphs, arranged in two columns:
##### (Left column, through-time analyses)
This column contains:
- A graph of the total number of occupied patches, through time.
- A graph of the total number of individuals in the simulation, over time.
- A graph of the number of individuals in a selected plot, through time. The plot can be changed in the box below the *Simulate Metapopulations* button on the top. Localized extinction events can be easily seen.
##### (Right column, patch analyses)
This column contains:
- A graph of the proportion of time each patch was occupied.
- A graph of the number of extinction events in each patch.
- A graph of the log10 Emigration:Immigration rates per patch, showing a proxy of whether patches are a source (high values) or a sink (low values).
#### Simulate Gene Flow
Selecting this button will run the gene flow portion of the simulation with the parameters you have selected on the **Simulation Settings** tab. This portion runs much, much slower than the metapopulation portion. It may take anywhere from a few to over 15 minutes on your machine. Once the simulation is complete, three graphs will appear:
- A graph of pairwise Fst.
- A graph of a Principle Coordinates Analysis based on the pairwise Fst values.
- A graph of isolation by distance, based on the map and the pairwise Fst values.
For all three graphs, users can select the output generation to plot. Please note that if you select a generation which is not on a save iteration, plots will fail!
#### Download Simulation
Selecting this button allows users to download the results of the simulation for further analyses, or use in the next round of simulation.

## Output Structure
Outputs are in R Data Serialization (RDS) format. They are a list object with four elements:
- riverscape: A list containing the input file and parameters for the riverscape.
- species: A list containing the parameters for the species.
- pop.hist: A list of detailed population history, containing:
  - A matrix of pre-generation abundances
  - A matrix of post-generation abundances
  - A list carrying the movement history among patches for individuals
  - A genealogy showing who reproduced with who
  - A matrix of abundance of females (sex-based simulation, effective female population size)
- GeneSim: A list of detailed genomic results, arranged as lists within each generation, containing lists of site-level genotypes in 3-D arrays (nested list structure)

## Known Issues
There are a few known issues with the simulation as it currently exists.
- Watershed with more or less than 25 patches may not plot correctly for either Fig 1 (watershed network diagram) or Fig 2 (patch connectivity diagram). 7
- Starting from a previous run occasionally blanks all parameters, for reasons currently unknown. Also, starting from a previous run sometimes results in the genomic simulation crashing. We are currently working to resolve this. The metapopulation simulation does not seem sensitive to this issue. 
- If the sampling interval ends on a number other than the last generation, the last generation will not be saved.
- Fst values are currently unrealistically high. This is likely due to small population sizes driving high rates of fixation. However, the general patterns of Fst do make biological sense, if not the magnitudes.

## Planned Improvements
GenFlow Mark 2 is currently under development. Due to the complexity of the RShiny App, the model is being built from the ground up rather than by trying piecewise substitution in existing architecture. Each planned improvement is already working in isolation or included in a different model and needs to be ported. Planned improvements include:
- Vastly improved run times, allowing for larger populations and more loci.
- A better algorithm to plot the watersheds.
- A disturbance generator which allows for 1) different degrees of disturbance, 2) tunable spatial autocorrelation of disturbance, and 3) fixed affected/unaffected patches. This should allow us to combine the currently separate press- and pulse- disturbance boxes.
- *K* calculated off of more than one habitat metric, and/or user input options.
- Improved run re-loading.
- An optional BACI (before-after control-impact) design setting, allowing users to specify the duration of the before and after periods.
- Better visualization tools, including more maps.
- A far more efficient data output format. 


