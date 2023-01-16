#!/bin/bash
for f in $(grep -li kill logs/mbl_100x100*); do
    grep Run $f | head -1 | cut --characters=1-7
    grep Run $f| tail -1 | cut --characters=1-7
    echo -e ""
done