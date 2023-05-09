#!/bin/bash

for i in {35..40}
do
    ./pcenter "../instances/pmed/pmed$i.txt" 1 p 1
    ./pcenter "../instances/pmed/pmed$i.txt" 1 p 2
    ./pcenter "../instances/pmed/pmed$i.txt" 1 p 3

done
