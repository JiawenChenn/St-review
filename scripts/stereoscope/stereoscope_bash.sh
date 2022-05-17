#!/bin/bash
scRNA_path="./processed_data/MOB/ref/internal/ob.sc_cnt.hvg2000.1950cell.tsv"

scRNA_anno_path="./processed_data/MOB/ref/internal/ob.sc.1950cell.mta.tsv"

output_dir="MOB.Seqfish.int_hvg2000"

st_path="./processed_data/MOB/seqfish.cnt.genexrow.tsv"


mkdir -p ./scripts/stereoscope/result/${output_dir}

stereoscope run --sc_cnt ${scRNA_path} --sc_labels ${scRNA_anno_path} -sce 30000  \
--st_cnt ${st_path} -stt -ste 30000 -stb 512 -scb 512 --gpu  --keep_noise   \
-o ./scripts/stereoscope/result/${output_dir}
