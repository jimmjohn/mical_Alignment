#!/bin/bash

source env.sh

ulimit -s unlimited

mkdir fileOut

./madurai_rpcdata_10l $1 $2 $3

mv RPC* fileOut/
# ./a.out 1>tmp.txt 2>tmpe.txt
