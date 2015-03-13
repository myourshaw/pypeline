#!/bin/bash 
echo "Good night $1!"
if [[ $2 != '' ]]; then s=${2}; else s=3; fi
sleep ${s}s
echo "Good morning $1! You were asleep ${s} seconds"

