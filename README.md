# Generator-Codes
On this page, the generator codes for the paper "Benchmarking for Robust Discrete Optimization Problems" by Marc Goerigk and Mohammad Khosravi can be found, see https://arxiv.org/abs/2201.04985v1


Instruction of using generator codes:

Each folder contains three files named main.cpp, selection.cpp and sel.h. The main.cpp file introduces the input parameters. The selection.cpp file consists of all the function we used and the sel file.h, has all the classes defined. The input parameters, which should be given through the command line are different for each problem. Therefore, they are described for all combination of problems in the following sections: 

--------------------------------------------------------------------------------------------------------
1. MIN-MAX PROBLEM WITH DISCRETE UNCERTAINTY SET:

In this setting the parameters must be given the the exact following order for the selection problem:

n: number of items

p: number of items to be selected

N: number of scenarios

param: method of instance generation (param=0 for sampling and param=1 for HIRO)

R: methods of sampling or starting the HIRO (R=0 for MM_D_U) (R=1 for MM_D_1) (R=2 for MM_D_2)

timelimit: the time limit for finishing the process, in particular for the HIRO

budget: maximum permitted adjustment for each member of a scenario

random seed: this allows the user to generate different instances when all other parametrs are the same 

--------------------------------------------------------------------------------------------------------
2. MIN-MAX PROBLEM WITH BUDGETED UNCERTAINTY SET:

In this setting the parameters must be given the the exact following order for the selection problem:

n: number of items

p: number of items to be selected

gamma: number of items that can exceed their nominal value

param: method of instance generation (param=0 for sampling, param=1 for HIRO with only changing lower bound, param=2 for HIRO with only changing deviation and param=3 for HIRO with changing both lower bound and deviation)

R: methods of sampling or starting the HIRO (R=0 for MM_B_U) (R=1 for MM_B_1) (R=2 for MM_B_2)

t: the time limit for finishing the process, in particular for the HIRO

b: maximum permitted adjustment for each member of a scenario

random seed: this allows the user to generate different instances when all other parametrs are the same 

--------------------------------------------------------------------------------------------------------
3. MIN-MAX REGRET PROBLEM WITH DISCRETE UNCERTAINTY SET:

In this setting the parameters must be given the the exact following order for the selection problem:

n: number of items

p: number of items to be selected

N: number of scenarios

param: method of instance generation (param=0 for sampling and param=1 for HIRO)

R: methods of sampling or starting the HIRO (R=0 for MMR_D_U) (R=1 for MMR_D_1) (R=2 for MMR_D_2)

timelimit: the time limit for finishing the process, in particular for the HIRO

budget: maximum permitted adjustment for each member of a scenario

random seed: this allows the user to generate different instances when all other parametrs are the same

--------------------------------------------------------------------------------------------------------
4. MIN-MAX REGRET PROBLEM WITH INTERVAL UNCERTAINTY SET:

In this setting the parameters must be given the the exact following order for the selection problem:

n: number of items

p: number of items to be selected

param: method of instance generation (param=0 for sampling and param=1 for HIRO)

R: methods of sampling or starting the HIRO (R=0 for MMR_I_U) (R=1 for MMR_I_1) (R=2 for MMR_I_2)

t: the time limit for finishing the process, in particular for the HIRO

b: maximum permitted adjustment for each member of a scenario

random seed: this allows the user to generate different instances when all other parametrs are the same 

--------------------------------------------------------------------------------------------------------
5. TWO-STAGE PROBLEM WITH DISCRETE UNCERTAINTY SET:

In this setting the parameters must be given the the exact following order for the selection problem:

n: number of items

p: number of items to be selected

N: number of scenarios

param: method of instance generation (param=0 for sampling, param=1 for HIRO with only changing the first stage scenario and param=2 for HIRO with changing both the first and second stage scenarios)

R: methods of sampling or starting the HIRO (R=0 for 2ST_D_U) (R=1 for 2ST_D_1) (R=2 for 2ST_D_2)

timelimit: the time limit for finishing the process, in particular for the HIRO

budget: maximum permitted adjustment for each member of a scenario

random seed: this allows the user to generate different instances when all other parametrs are the same 

--------------------------------------------------------------------------------------------------------
6. TWO-STAGE PROBLEM WITH DISCRETE BUDGETED UNCERTAINTY SET:

In this setting the parameters must be given the the exact following order for the selection problem:

n: number of items

p: number of items to be selected

gamma: number of items that can exceed their nominal value

param: method of instance generation (param=0 for 2ST_DB_U) (param=1 for 2ST_DB_1) (param=2 for 2ST_DB_2)

random seed: this allows the user to generate different instances when all other parametrs are the same 

--------------------------------------------------------------------------------------------------------
7. TWO-STAGE PROBLEM WITH CONTINUOUS BUDGETED UNCERTAINTY SET:

In this setting the parameters must be given the the exact following order for the selection problem:

n: number of items

p: number of items to be selected

gamma: number of items that can exceed their nominal value

param: method of instance generation (param=0 for 2ST_CB_U) (param=1 for 2ST_CB_1) (param=2 for 2ST_CB_2)

random seed: this allows the user to generate different instances when all other parametrs are the same 

--------------------------------------------------------------------------------------------------------
8. RECOVERABLE PROBLEM WITH DISCRETE UNCERTAINTY SET:

In this setting the parameters must be given the the exact following order for the selection problem:

n: number of items

p: number of items to be selected

N: number of scenarios

delta: the maximum number of items wich can be changed in the second stage

param: method of instance generation (param=0 for sampling, param=1 for HIRO with only changing the first stage scenario and param=2 for HIRO with changing both the first and second stage scenarios)

R: methods of sampling or starting the HIRO (R=0 for RR_D_U) (R=1 for RR_D_1) (R=2 for RR_D_2)

timelimit: the time limit for finishing the process, in particular for the HIRO

budget: maximum permitted adjustment for each member of a scenario

random seed: this allows the user to generate different instances when all other parametrs are the same 

--------------------------------------------------------------------------------------------------------
9. RECOVERABLE PROBLEM WITH DISCRETE BUDGETED UNCERTAINTY SET:

In this setting the parameters must be given the the exact following order for the selection problem:

n: number of items

p: number of items to be selected

gamma: number of items that can exceed their nominal value

param: method of instance generation (param=0 for RR_DB_U) (param=1 for RR_DB_1) (param=2 for RR_DB_2)

random seed: this allows the user to generate different instances when all other parametrs are the same 

--------------------------------------------------------------------------------------------------------
10. RECOVERABLE PROBLEM WITH CONTINUOUS BUDGETED UNCERTAINTY SET:

In this setting the parameters must be given the the exact following order for the selection problem:

n: number of items

p: number of items to be selected

gamma: number of items that can exceed their nominal value

param: method of instance generation (param=0 for RR_CB_U) (param=1 for RR_CB_1) (param=2 for RR_CB_2)

random seed: this allows the user to generate different instances when all other parametrs are the same 

--------------------------------------------------------------------------------------------------------



This page is still under construction and will be updated regularly.
