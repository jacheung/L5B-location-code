# L5bLocationCode

## Premise: 
How are touch features represented in the brain? Specifically we seek in this project to understand how the cortex represents object location. Object location representation is a sensorimotor transformation that integrates information of WHEN an object was touched with WHERE the sensor was at touch. 

Leveraging the anatomical, behavioral, and technological tools available in the mouse whisker system we recorded single neurons in L5b of mouse S1 while tracking high speed video of whisker motion.  

Below are scripts used to characterize encoding and test touch location decoding using single neurons. For summary of all results and findings using this script see Cheung et al., 2019b. For a complete description of the behavioral paradigm during neural recordings refer to Cheung et al., 2019a. 

## Requirements: 
All code is built and tested on MATLAB 2018b. 
glmnet from Qian and Hastie 2013 - https://web.stanford.edu/~hastie/glmnet_matlab/
Dataset structures with neural recordings and stimulus variables - Request datasets by emailing jacheung6@gmail.com



## Main figure builders: 

#### hilbert_plots:
This script compiles the optimized binomial models of neural activity around touch and generates the below figures for analysis. 
![Alt text](./pictures/sampleModeledHilbert.png)