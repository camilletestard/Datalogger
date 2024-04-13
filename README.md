# Neural signatures of natural behavior in socializing macaques

This repository contains the code associated with our research on the neural underpinnings of natural social behavior in primates (rhesus macaques). To access to raw data, go to: https://osf.io/e2xsu/

Our work, titled "_Neural signatures of natural behavior in freely-socializing macaques_", was conducted by Testard & Tremblay et al. (2023) and can be accessed [here](https://www.biorxiv.org/content/10.1101/2023.07.05.547833v1). 

### Abstract

Our understanding of the neurobiology of primate behavior largely derives from artificial tasks in highly-controlled laboratory settings, overlooking most natural behaviors primate brains evolved to produce1. In particular, how primates navigate the multidimensional social relationships that structure daily life and shape survival and reproductive success remains largely unexplored at the single neuron level. Here, we combine ethological analysis with new wireless recording technologies to uncover neural signatures of natural behavior in unrestrained, socially interacting pairs of rhesus macaques within a larger colony. Population decoding of single neuron activity in prefrontal and temporal cortex unveiled robust encoding of 24 species-typical behaviors, which was strongly modulated by the presence and identity of surrounding monkeys. Male-female partners demonstrated near-perfect reciprocity in grooming, a key behavioral mechanism supporting friendships and alliances, and neural activity maintained a running account of these social investments. When confronted with an aggressive intruder, behavioral and neural population responses reflected empathy and were buffered by the presence of a partner. Surprisingly, neural signatures in prefrontal and temporal cortex were largely indistinguishable and irreducible to visual and motor contingencies. By employing an ethological approach to the study of primate neurobiology, we reveal a highly-distributed neurophysiological record of social dynamics, a potential computational foundation supporting communal life in primate societies, including our own.

### Repo Contents

This repository contains all the code used for the analysis and visualization of the behavioral and neural data collected in our study. Most of the code is in Matlab with some pieces in R (behavioral data visualization only). 
The code is organized into several directories which are split by figure content and topic. We also include a folder for external toolboxes used in the code and a folder for preprcossing the data. Each script takes the preprocessed and formatted data as input (starts by running the function "log_GenerateDatToRes_function.m"). Make sure to add the full github directory to your path in Matlab. 

If you have any questions about the code or data please reach out to Dr. Camille Testard at ctestard@fas.harvard.edu or camille.testard94@gmail.com
