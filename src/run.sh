#!/bin/bash

# Accumulate runtime of the five runs
echo "runtimes" > res.csv
for i in {1..5}
do          # Grab the total time
    ./MD | awk '/500 timesteps / {print $4}' >> res.csv;
done