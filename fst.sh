#!/usr/bin/env bash

# cal
vcftools --gzvcf $v_f \
--weir-fst-pop $pop1 --weir-fst-pop $pop2 \
--out ovis 

# kk
