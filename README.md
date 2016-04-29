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

* change the seed initialization to use milliseconds or even micro or nano seconds
* figure out how to have literal back-quotes in this markdown file