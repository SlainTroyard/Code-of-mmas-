# code-of-mmas

**https://doi.org/10.1016/j.knosys.2021.108000**  

MMAS is a solver for the Maximum Set k-Covering Problem.  

## The input format

The input instance is asked to be in ascii format  

See the data folder for detailed examples.  


## Compile

MMAS is implemented in C++ and complied by g++ with '-O3' option.  

eg:   
```
g++ mmas.cpp -O3 -o mmas
```

## Usage

MMAS is a local search and break ties randomly which needs a random seed, so the command to run MMAS is:   
```
./mmas instance_name k_value random_seed_value time_limit
```

For example:  
```
./mmass scp41.txt 34 44 100
```

`scp41.txt` is the name of instance and the random seed is set to 44. The algorithm will run for 100s.  
