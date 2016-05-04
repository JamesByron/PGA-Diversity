# Diversity in Genetic ALgorithms

## Synopsis

Genetic Algorithm for measuring diversity in a population and trying to foster it.

## Motivation

To be determined

## Running

Minimally, run ./singlenode

On the cluster, to run on multiple compute-nodes,
rocks run host compute 'cd PGA-Diversity; ./singlenode insData logfileflag `expr \`date +%s%N\` % 1000000000`'

## To Do

* reduce the volume in the log files: remove the repeated words but ensure a valid key is maintained that identifies which column represents which data
* change the seed initialization to use milliseconds or even micro or nano seconds
* figure out how to have literal back-quotes in this markdown file