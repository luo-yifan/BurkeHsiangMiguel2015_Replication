
Readme for replicating Burke, Hsiang, and Miguel 2015

The data and code in the zipped directory BurkeHsiangMiguel2015_Replication.zip replicate all figures, tables, and estimates in Burke, Hsiang, Miguel 2015 (with the exception of Figure 1, which is constructed from figures in the literature).  

Full replication requires both R and Stata.  User-written programs for both R and Stata are required at various points, as highlighted in script headers where relevant.  Code was developed/implemented in Stata/MP 14, and R 3.2.0. 

Users of these data should cite Burke, Hsiang, Miguel 2015.  If you find a consequential error, or alternatively simply wish to congratulate us on a job well done, please email Marshall at mburke@stanford.edu.


Info on folders contained in replication data
—data/input, where input data is stored
—data/output, where generated data is stored that is then used to make figures.  
—figures/MainFigs_Input  where the panels for the figures in the main text are stored
—figures/ExtendedDataFigs_Input  where the panels for the figures in the extended data are stored
—script/  where scripts are stored

The main data file for the regressions is provided in both Stata 14 format (.dta), Stata 13 format (_Stata13.dta) as well as .csv. When importing the .csv into Stata, the option “case(preserve)” is needed to maintain consistency with the provided .do files.  

Instructions for replicating all figures/tables in main text and Extended Data
For all Stata and R code, user needs to set working directory to the directory where the files were unzipped.

1. First run all the Stata do files in the script/ folder :
	a. GenerateFigure2Data.do
	b. GenerateBoostrapData.do
	c. MakeExtendedDataFigure1.do, MakeExtendedDataFigure2.do, MakeExtendedDataFigure3.do, MakeExtendedDataTable1.do, MakeExtendedDataTable2.do

2.  Then run the R scripts generating the impact projections
	a. getTemperatureChange.R
	b. ComputeMainProjections.R
	c. ComputeDamageFunction.R

3.  Then run the desired MakeFigure or MakeExtendedDataFigure R script.


Instructions for generating specific figures in the main text
For all Stata and R code, user needs to set working directory to the directory where the files were unzipped.

To generate Figure 2, the user needs to run
1.  GenerateFigure2Data.do
2.  MakeFigure2.R

To generate Figure 3 or 4, the user needs to run:
3.  GenerateBootstrapData.do
4.  GetTemperatureChange.R  
5.  ComputeMainProjections.R
6.  And then run the desired MakeFigure3.R or MakeFigure4.R

To generate Figure 5, the user needs to run 3-5 above, and then run:
7. ComputeDamageFunction.R
8.  MakeFigure5.R


