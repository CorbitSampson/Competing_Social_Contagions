####################################################################################################################################
	This repository contains all code used for the generation of data and figures in 
		C.R. Sampson and J.G. Restrepo, Competing Social Contagions with Opinion Dependent Infectivity

	It also contains the figures as they appear in the article as .png and .pdf files and the data needed to generatate 
	figures 4, 5, 7, 10, 11, and 12.
#####################################################################################################################################
	Dependences and Acknowledgments:

		This could relies on a number of standard packages:
			random
			pandas
			numpy
			scipy
			time
			matplotlib
			IPython
		
		as well as the CompleX Group Interactions (XGI) package which is availble here: https://github.com/xgi-org/xgi

		N.W. Landry, M. Lucas, I. Iacopini, G. Petri, A. Schwarze, A. Patania, and L.Torres, 
		XGI: A Python package for higher-order interaction networks, Journal of Open Source Software 8, 5162 (2023).

		The primary code to build the figures is a .ipynb file written to be viewed using jupyter notebook.
######################################################################################################################################
	This repository includes the files simulation.py, CNO.py, and reduced_meanfield_object.py which are the primary objects 
	which implement our model. meanfield_object.py contains a non-reduced version of our mean-field model and is no longer in use. 
	But has been included for completeness. The opinion_generation.py file contains the code to generate initial opinions pulled from
	a truncated Gaussian distribution and a bimodal Dirac delta distribution, used for figures 10 and 11, respectivly. 

	The files labeled as "data_generators" create the data needed to generatate figures 4, 5, 7, and 10. All figures without 
	data generators are created directly from Figure_Builders.ipynb, which is a jupyter notebook file.

	The data generated in each of the data generator files has also been included, as several of this data generators have 
	a long run time. This allows the contents of Figure_Builders.ipynb to be run without running each of the data generator files.

	Lastly, all figures have been included either as .png or .pdf files, as used in the article.

	Any questions regarding the code, data, or figures in this repository can be directed to corbit.sampson@colorado.edu
########################################################################################################################################