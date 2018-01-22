#!/bin/csh

grep -A 1 "Energy initial, next-to-last, final" log.lammps | tail -n 1 | awk '{print $3}'
