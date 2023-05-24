# Sequentially planned probability ratio tests 

This repository contains an R module for optimal design and performance evaluation of sequentially planned probability ratio tests.

Currently, only the model of binary responses (sampling from a Bernoulli population) is covered.

The detailed description of the method and algorithms can be found in

[*Andrey Novikov. Group sequential hypothesis tests with variable group sizes: optimal design and performance evaluation, 2023. Communications in Statistics - Theory and Methods, to appear*] (https://arxiv.org/abs/2210.07203)

## Content description
* The file [SPPRT.R](SPPRT.R) contains all the functions providing the  user interface for all the tasks.

The list of functions can be seen below. 

### design_test

The function for designing an optimal sequentially planned probability ratio test

Arguments
* _l0_, _l1_ Lagrange multipliers
* _th0_, _th1_ the hypothesized parameter values
* _H_ horizon (maximum number of steps the test can use)
* _gsizes_ the vector containing eligible group sizes 
* _gamma_ weight for calculating Average Sampling Cost
* _h_ grid size 
* _step_fn_ function to be called in any stage of the test evaluation 

Returns the designed test 


### operating_characteristic 

The function calculating the operating characteristic of a test

Arguments
* _test_ test designed by the _design_test_ function
* _th_ parameter value at which the operating characteristic is calculated
* _step_fn_ function to be called in any stage of the evaluation

### sampling_report

The function calculating the test characteristics of a test
* _test_ test designed by the _design_test_ function
* _th_ parameter value at which the characteristic is calculated
* _step_fn_ function to be called in any stage of the evaluation
* _accounting_fn_ function to be defined to calculate the desired characteristic,
  by default, it is predefined as function _cost_ accounting for the Average Sampling Cost.

  Other variants: Average Number of Groups, Average Number of Observations (see  examples in [usage.R](usage.R))

 
* The file [usage.R](usage.R) is a usage example of these functions (corresponds to example 4.2 in the cited [*article*](https://arxiv.org/abs/2210.07203)).
 
