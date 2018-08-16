# Resource-Mapping-Tool
This work is covered by a GNU GPL v3.0, as provided in this repository (titled 'LICENSE').

This work was produced as a part of the Renewable Systems Engineering (RENESENG) grant, as funded by the Research Executive Agency of the Seventh Framework Programme of the European Union, as a Marie Curie Initial Training Network (FP7-PEOPLE-2013-ITN), with grant number 607415.

This project is for use with QGIS, written and tested using version 2.18.7 (I have not tested other versions, but hopefully will in the future). The work is written for python 2.7 (I have not tested it for compatibility with python 3, but again, hopefully will in the future).

This code is meant to assist in understanding the distribution of resources and producing the data sets necessary for a biomass supply chain optimisation

To use, either 1) save the code file and select the "Add script from file" tool in the Processing Toolbox, or 2) create a new blank script, copy-paste the code in, and save the script. To run the script select the "Run Algorithm" option on the script editor screen (looks like interlaced gears) and provide the required information.

This program creates a cellular grid for use in biomass supply chain optimisation. It then uses that grid to aggregate data from each raster-based data set provided by the user for use in the supply chain optimisation.

This script requires:

(1) a raster map of land use for the region of interest (2) the no data value used in the land use map (3) the list of raster values corresponding to land uses of interest (4) a raster map of constraints on the land (5) up to 5 raster maps of biomass yield for the region of interest (6) the product name for each yield map (7) a shapfile that outlines the area of interest (8) the desired number of cellular elements in the grid (9) the minimum relative size of cells to be used in the final grid (10 & 11) the x-value fraction at which the break points for the linear estimators occur


This script outputs:

(1) a matrix with cell ID numbers and the area of the cell (2) a matrix with the linear distance between all pairs of cells (3) a matrix with the cell ID and usable land area of the cell (4) a matrix containing the tortuosity between all cells (5) a matrix with the cell ID and and average bioyeld per cell (6) the shapefile of the cellular grid that was used to aggregate the data (7-9) the slope of the biomass yield approximations (10-12) the intercept of the biomass yield approximations 
