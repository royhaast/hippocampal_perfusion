#!/bin/bash

source /home/rhaast/venv/brainstat/bin/activate

SUBJECT=$1
HEMISPHERE=$2

export SUBJECT HEMISPHERE

jupyter nbconvert --ExecutePreprocessor.timeout=None --to html --output `realpath $3` --execute $4
