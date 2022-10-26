"""
This code is used to find stable motifs for thesis section 3.5. Adapted by Daniel Plaugher on 9/14/22
"""

import pystablemotifs as sm
import pyboolnet
import pystablemotifs.export as ex
import networkx as nx
from timeit import default_timer

# %% you can import rules as a string
Rules=''' x1* = x3 | x2
x2* = x1 & !x3
x3* = !x2 | !x1'''

PRIMES=sm.format.create_primes(Rules)

# %% Or input rules from a separate file
relative_path_to_model = "./models/SM_example.txt"
primes = sm.format.import_primes(relative_path_to_model)

# to visualize rules
print("RULES")
sm.format.pretty_print_prime_rules(primes)
print()


# to show attractors
ar = sm.AttractorRepertoire.from_primes(primes)
ar.summary()

df=ex.attractor_dataframe(ar)
df

# %%
# define the desired target
target={'x3':0}


## Brute Force
print("--- BRUTE FORCE ---")
start=default_timer()
interventions = sm.drivers.knock_to_partial_state(target,primes,max_drivers=2)
end=default_timer()
print("Time running method:",end-start)
print("Sets found:")
for x in interventions: 
    print({k:v for k,v in sorted(x.items())})
print()   
    
 
## Grasp Search
print("--- GRASP SEARCH ---")
GRASP_iterations=2000
start=default_timer()
interventions = sm.drivers.GRASP(target,ar.primes,GRASP_iterations)
end=default_timer()
print("Time running method:",end-start)
print("Sets found:")
for x in interventions: 
    print({k:v for k,v in sorted(x.items())})
print()    


## Internal History
print("--- INTERNAL HISTORY ---")
start=default_timer()
interventions = ar.reprogram_to_trap_spaces(target,
                                            target_method='history',
                                            driver_method='internal')
end=default_timer()
print("Time running method:",end-start)
print("Sets found:")
for x in interventions: print({k:v for k,v in sorted(x.items())})
print()

## Minimal History
print("--- MINIMAL HISTORY ---")
start=default_timer()
interventions = ar.reprogram_to_trap_spaces(target,
                                            target_method='history',
                                            driver_method='minimal')
end=default_timer()
print("Time running method:",end-start)
print("Sets found:")
for x in interventions: 
    print("---")
    print("One temporary intervention from each list, in order.")
    print("("+str(len(x))+" interventions in total)")
    for y in x: print(y,"\n")