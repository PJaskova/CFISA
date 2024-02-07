# Compositional functional isotemporal substitution analysis (CFISA)

The proposed method CFISA is performed to assess how varying the weight of a specific range of density influences response variable. 

## Requirements

Calculations were performed in following software, using these libraries:
- R version: 4.2.1 (2022-06-23) -- "Funny-Looking Kid"
- Libraries:
   - robCompositions: log-ratio transformations; version 2.4.1
   - fda: statistical analyses of functions; version 6.1.4
   - car: logit transformation; version 3.1-2

## Available scripts
 
1. function.R:   
      a)  clr2density - back transformation from L^2 to B^2

      b) ZsplineBasis - basis with zero integral
   
      c)SmoothingSpline01 - smoothing spline
   
    All these functions was introduced in "Talská, R., Hron, K. & Grygar, T.M. Compositional Scalar-on-Function Regression with
    Application to Sediment Particle Size Distributions. Math Geosci 53, 1667–1695 (2021). https://doi.org/10.1007/s11004-021-09941-1"
   
      d)bpc1 - backwards pivot coordinates
   
      e) is - main function to calculate isotemporal substitution
2. CFISA.R

        
 
 
