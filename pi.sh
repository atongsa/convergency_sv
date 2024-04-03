#!/usr/bin/env bash

vcftools --gzvcf $v_f --window-pi --window-pi-step --keep $po_p --out $op_p

