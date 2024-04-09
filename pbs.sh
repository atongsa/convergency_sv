#!/usr/bin/env bash
pbscan -win 1 -step 1 -div 1 -vcfp sel.vf \
    -pop1 $p_id -pop2 $n_id -pop3 $ref_id \
    -out pbs_can/$n_n

