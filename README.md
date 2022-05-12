# Generator-Codes
On this page, the generator codes for the paper "Benchmarking for Robust Discrete Optimization Problems" by Marc Goerigk and Mohammad Khosravi can be found, see https://arxiv.org/abs/2201.04985v1


Instruction of using generator codes:

Each folder contains three files named main.cpp, selection.cpp and sel.h. The main.cpp file introduces the input parameters. The selection.cpp file consists of all the function we used and the sel file.h, has all the classes defined. The input parameters, which should be given through the command line are different for each problem. Therefore, they are described for all combination of problems in the following sections: 

1. MIN-MAX PROBLEM WITH DISCRETE UNCERTAINTY SET:

In this setting the parameters must be given the the exact following order for the selection problem:

n: number of items
p: number of items to be selected
N: number of scenarios
param: method of instance generation (param=0 for sampling and param=1 for HIRO)
R: methods of sampling or starting the HIRO (R=0 for MM_D_U) (R=1 for MM_D_1) (R=2 for MM_D_2)
timelimit: the time limit for finishing the process, in particular for the HIRO
budget: maximum permitted adjustment for each member of a scenario
seed: this allows the user to generate different instances when all other parametrs are the same 


2. MIN-MAX PROBLEM WITH BUDGETED UNCERTAINTY SET:

3. MIN-MAX REGRET PROBLEM WITH DISCRETE UNCERTAINTY SET:

In this setting the parameters must be given the the exact following order for the selection problem:

n: number of items
p: number of items to be selected
N: number of scenarios
param: method of instance generation (param=0 for sampling and param=1 for HIRO)
R: methods of sampling or starting the HIRO (R=0 for MM_D_U) (R=1 for MM_D_1) (R=2 for MM_D_2)
timelimit: the time limit for finishing the process, in particular for the HIRO
budget: maximum permitted adjustment for each member of a scenario
seed: this allows the user to generate different instances when all other parametrs are the same

4. MIN-MAX REGRET PROBLEM WITH INTERVAL UNCERTAINTY SET:

5. TWO-STAGE PROBLEM WITH DISCRETE UNCERTAINTY SET:

6. TWO-STAGE PROBLEM WITH DISCRETE BUDGETED UNCERTAINTY SET:

7. TWO-STAGE PROBLEM WITH CONTINUOUS BUDGETED UNCERTAINTY SET:

8. RECOVERABLE PROBLEM WITH DISCRETE UNCERTAINTY SET:

9. RECOVERABLE PROBLEM WITH DISCRETE BUDGETED UNCERTAINTY SET:

10. RECOVERABLE PROBLEM WITH CONTINUOUS BUDGETED UNCERTAINTY SET:



This page is still under construction and will be updated regularly.
