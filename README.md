# Calgyp
Calgyp is a model designed to calculate the the concentration and dispersion of gypsum (CaSO4 H2O)  in a soil column. The model is given a daily rainfall series supplied by a stochastic weather generator and the flux of sulfate (SO4) and calcium (Ca).  
  
  

> gypsum accumulation in a soil profile with time. 

![asd](https://media.giphy.com/media/BGrb9vc3Eb7t9A0hug/giphy.gif)
 
## Operation

 1. Weather Generator
	 This section is responsible for creating a daily rainfall series of costume length which represents the statistical climate properties in a selected meteorological station. 
	 *evapostats\RainGen.R*: Is responsible for getting the meteorological data from the IMS using *GetImsRain()* and generate the wet and dry days probabilities using markov chain algorithm  in *CalculateProbabilities()*. The *GenerateRainSeries(numOfSeries)* function generates the costume length series.
	 
 2. Soil hydrological model
	A c++ model which simulate a soil column divided into compartments of a specified thickness. Every iteration represents a daily routine of moisture addition and evaporation according to the daily rain ane potential evaporation supplied by the weather generator. The model is located in the *Calcyp* directory
	 - CSM.cpp - The main function is  *CSM::Calculate()* responsible for calculating the daily water balance. Distribute the moisture between the compartments represent the soil profile. 
	 - Compartment.cpp - represent a compartment object. The properties of the compartment are soil parameters, moisture, and gypsum concentration. *solubility()* calculates the equilibrium and return the accumulated gypsum
3.  Model Analysis 
	Functions.R contains the function *CalcGypsum()* which wrap the *CSM::calculate()* in an R environment  
   
```
 print("run model")
Rcpp::sourceCpp('C:/Users/liorst/source/repos/evapocalc/Calcyp/CSM.cpp', verbose = TRUE, rebuild = 0);
cppModule <- new(CSMCLASS);
list = cppModule$Calculate(raindata$rain, raindata$PET, duration, Depth, thick, wieltingPoint, nArea, FieldCapacity, DustGyp, AETFactor, verbose, dustFlux / 10000 / 365, rainCa, rainSO4,withFC);
 ```
 This allows for fast calculation alongside comfortable model analysis.

# Output and results

    

    

	

    

