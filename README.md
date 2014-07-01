
# maxSatUZK

Build: 11

Authors: Alexander van der Grinten, Andreas Wotzlaw

University of Cologne

## About

maxSatUZK is a solver for MaxSAT and its weighted and partial variants.
The solver is based on binary search using our SAT solver satUZK as backend.
Constraints are encoded to CNF through sorting networks.

## How to build

Dependencies:
- satUZK (https://github.com/satuzk/satUZK)

Create a folder named `satuzk` and copy the satUZK source code into it.
Run `make` to build the solver.

## How to use

Run the solver using:
`./maxSatUZK <instance file>`

