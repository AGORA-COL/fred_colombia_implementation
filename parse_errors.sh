#!/bin/bash

ERROR_FILE=$1
MISSING_DEPS=$(grep -oP '(?<=Try installing: ).*' $ERROR_FILE | tr -d '*' | tr '\n' ' ' | sed 's/  /, /g')

if [[ -n "$MISSING_DEPS" ]]; then
    echo "The following system dependencies may be missing based on R packages installation errors:"
    echo $MISSING_DEPS
else
    echo "No missing system dependencies were detected based on the R packages installation errors."
fi