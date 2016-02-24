	##################################
##########  compile the program   ############
	##################################

in the terminal:

g++ --std=c++11 -o test Source.cpp

---------------------------------------------------------------------------------------------------
	##################################
##########  Running the program   ############
	##################################

./test -inst <graph.name> -alg [BnB|Approx|LS1|LS2] -time <cuttoff_in_seconds> -seed <random_seed> 

######  example:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
./test -inst ../Data/karate.graph -alg LS2 -time 2 -seed 1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

NOTE:
  (1) Algorithm_name should be
                  1) BnB
									2) Approx
									3) LS1
									4) LS2
  (2) Before running the code, you should have created the folders as follows:
      Code
       |- Source.cpp
       |- readme.txt
      Ouput
       |- Solution
          |- BnB
          |- Approx
          |- LS1
          |- LS2
       |- Trace 
          |- BnB
          |- Approx
          |- LS1
          |- LS2

---------------------------------------------------------------------------------------------------
	##################################
  ############  result        ############
	##################################
In the folder Solution: to see the solution files
In the floder Trace: to see the trace files
