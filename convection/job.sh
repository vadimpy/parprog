#!/bin/bash

#PBS -I walltime=00:10:00,nodes=7:ppn=1
#PBS -N sum
#PBS -q batch

mpirun --oversubscribe -np 7 convection

