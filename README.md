# MIXED (DFNDFL)

 ```
 usage: run_all [-h] [-v] [--constrained] [--alg {DFN_DFL,DFL}] [--max_fun MAX_FUN]
                [--outlev OUTLEV] [-M NM_MEMORY]

 Run DFO codes

 optional arguments:
   -h, --help            show this help message and exit
   -v, --version         Print version number
   --constrained         Solve constrained problems
   --alg {DFN_DFL,DFL}
                         Name of algo to be run
   --max_fun MAX_FUN     Maximum number of function avaluations
   --outlev OUTLEV       Output level
   -M NM_MEMORY, --NM_memory NM_MEMORY
                         History size for nonmonotone linesearch
 ```


##### AUTORI: T. Giovannelli<sup>1</sup>, G. Liuzzi<sup>2</sup>, S. Lucidi<sup>2</sup>, F. Rinaldi<sup>3</sup>

<sup>1</sup> Department of Industrial and Systems Engineering, Lehigh University (tog220@lehigh.edu)

<sup>2</sup> Department of Computer, Control and Management Engineering, "Sapienza" University of Rome (liuzzi@diag.uniroma1.it, lucidi@diag.uniroma1.it)

<sup>3</sup> Department of Mathematics "Tullio Levi-Civita", University of Padua (rinaldi@math.unipd.it)


Copyright 2021


## Figures in the paper

![run su problemi NON vincolati (34) con M=0](DFNDFL_DFL_NOMADnomod.png "Unconstrained (34,M=0)")
*Profiles on the 34 unconstrained problems of methods not using models.*

![run su problemi NON vincolati (34) con M=0](DFNDFL_NOMADmod_RBFOpt_MISO.png "Unconstrained (34,M=0)")
*Profiles on the 34 unconstrained problems of methods using models.*

![run su problemi vincolati (204) con M=0](DFNDFL_DFL_NOMADnomod_CON.png "Constrained (204,M=0)")
*Profiles on the 204 constrained problems of methods not using models.*

![run su problemi vincolati (204) con M=0](DFNDFL_NOMADmod_CON.png "Constrained (204,M=0)")
*Profiles on the 204 constrained problems of methods using models.*


## Instructions
1. at python prompt execute
2. python run_all.py
