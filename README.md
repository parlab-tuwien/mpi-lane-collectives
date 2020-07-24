
# MPI Multi-Lane Collectives

## Citing article

**Jesper Larsson Träff** and **Sascha Hunold**. 2019. *Decomposing MPI Collectives for Exploiting Multi-lane Communication*. In CLUSTER, 2020, IEEE.

## Build Library and Test Code

The file `tuw_lanecoll_bench_all.c` contains some test code, from which we stripped the entire timing code.
(The timing code use some in-house, non-portable way of benchmarking specific functions. For the sake of portability, we removed these bits of code.)
Now, you can see how to set up the benchmarks and how to run the collectives. You would simply need to add your custom timing code.

You can build the software like this (when `gcc` and `mpicc` are in `$PATH`). 
```
cmake ./
make
```

If you need to customize you build process, see the file in `platform_files`
to modify compiler parameters. Then, you could run for example
```
cmake -DINCLUDE_PLATFORM_CONFIG_FILE=./platform_files/intel_compiler.cmake  .
make
```

## Running on local machine

### Testing various collectives

#### MPI_Bcast

```
# for the multi-lane implementation
mpirun -np 4 ./bin/lanetest_all -m bcast -t lane -c 10 

# for the hierarchical implementation
mpirun -np 4 ./bin/lanetest_all -m bcast -t hier -c 10

# for the default implementation
mpirun -np 4 ./bin/lanetest_all -m bcast -t default -c 10
```

#### MPI_Scan

```
mpirun -np 4 ./bin/lanetest_all -m scan -t lane -c 10 
mpirun -np 4 ./bin/lanetest_all -m scan -t hier -c 10 
mpirun -np 4 ./bin/lanetest_all -m scan -t default -c 10 
```

#### MPI_Allreduce

```
mpirun -np 4 ./bin/lanetest_all -m allreduce -t lane -c 10 
mpirun -np 4 ./bin/lanetest_all -m allreduce -t hier -c 10 
mpirun -np 4 ./bin/lanetest_all -m allreduce -t default -c 10 
```

#### MPI_Allgather

```
mpirun -np 4 ./bin/lanetest_all -m allgather -t lane -c 10 
mpirun -np 4 ./bin/lanetest_all -m allgather -t hier -c 10 
mpirun -np 4 ./bin/lanetest_all -m allgather -t default -c 10 
```
