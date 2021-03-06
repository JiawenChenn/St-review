#!/bin/bash
scRNA_path="./processed_data/MOB/ref/internal/ob.sc_cnt.hvg2000.1950cell.tsv"

scRNA_anno_path="./processed_data/MOB/ref/internal/ob.sc.1950cell.mta.tsv"

output_dir="MOB.Seqfish.int_hvg2000"

st_path="./processed_data/MOB/seqfish.cnt.genexrow.tsv"

mkdir -p ./scripts/destVI/result/${output_dir}

python ./scripts/destVI/cell2location.py \
--sc_cnt ${scRNA_path} --sc_labels ${scRNA_anno_path} \
--st_cnt ${st_path} \
-out ./scripts/Adroit/result/ \
-pre ${output_dir}