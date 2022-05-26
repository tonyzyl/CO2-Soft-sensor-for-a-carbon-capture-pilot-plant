# Soft-sensor-for-estimating-CO2-for-a-carbon-capture-pilot-plant
A hybirid mechanistic and data-driven (Denoise Autoencoder (DAE) - Long-short term Memory (LSTM)) model for estimating the CO2 concentration profile for the CO2 absorber of the carbon capture plant. The tested dimensionality reduction techniques include: DAE, principal component analysis (PCA) / propoer orthogonal decomposition (POD). The absorber has 6 sampling points while there is only one gas analyzer, this work aims to build a soft sensor to use the other continuous process record to infer the CO2 concentration at these points. The DAE-LSTM approach is then compared to a semisupervised estimation method. The paper related to this work was submitted to Computers in Industry and it is currently in revision.

![Process scheme](images/pilot_plant.jpg?raw=true)

Main documents:

Estimate_CO2_profile.ipynb

semi_supervised_stacked_AE/pilot_plant_ssae.ipynb

Related similar approach:

X. Zhang, et.al, A weighted auto regressive LSTM based approach for chemical processes modeling, Neurocomputing, Volume 367, 2019, Pages 64-74, ISSN 0925-2312, https://doi.org/10.1016/j.neucom.2019.08.006.

The semisupervised stacked autoencoder training strategy is proposed in:

X. Yuan, et al., A novel semi-supervised pre-training strategy for deep networks and its application for quality variable prediction in industrial processes, Chemical Engineering Science, Volume 217, 2020, 115509, ISSN 0009-2509, https://doi.org/10.1016/j.ces.2020.115509.