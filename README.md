# Diversity in Genetic ALgorithms

## Synopsis

Genetic Algorithm for measuring diversity in a population and trying to foster it.

## Motivation

To be determined

## Running

Minimally, run ./singlenode

### on the cluster

(Note: see the source of this README file to get the correct syntax.)

To run on all compute-nodes and include a specific log-file-tag,
rocks run host compute 'cd PGA-Diversity; ./singlenode insData logfileflagwithhostname`hostname` `expr \`date +%s%N\` % 1000000000`'
Note: the hostname may be important since otherwise the log files will all have the same name and become a mess.

To run on selected compute-nodes,
use a file that contains the compute-nodes' names, one per line.
For example, if the file chassis0.0-8 contains the lines compute-0-0 through compute-0-8,
then the following command would run singlenode on those nine nodes:
rocks run host `< chassis0.0-8` 'cd PGA-Diversity; ./singlenode insData logfileflag`hostname` `expr \`date +%s%N\` % 1000000000`'
Naturally, if your code repository and executable is in a directory named something other than PGA-Diversity,
you would substitute your path there.

## To Do

* change the seed initialization to use milliseconds or even micro or nano seconds
* figure out how to have literal back-quotes in this markdown file
* (further) reduce the volume in the log files: remove the repeated words but ensure a valid key is maintained that identifies which column represents which data; this was partially done by shortening existing words without altering the field numbering
