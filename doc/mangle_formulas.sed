#!/usr/bin/env bash

gsed -e 's/^\$\$$/\\f\[/g' -e 's/^[ ]\$\$$/\\f\]/g' $1
