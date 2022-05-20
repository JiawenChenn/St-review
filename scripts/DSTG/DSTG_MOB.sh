# seqFISH+ MOB internal default gene set

module add python/3.7.9
module add r/4.0.3
module add gcc/4.9.1
module add clang/6.0

cd DSTG
# converted ./processed_data/MOB/seqfish.cnt.genexrow.tsv, ./processed_data/MOB/ref/internal/OB.sc.cnt.tsv and ./processed_data/MOB/ref/internal/OB.sc.mta.tsv to RDS
Rscript convert_data.R ./processed_data/MOB/ref/internal/scRNA_seqfish_count.RDS ./processed_data/MOB/st_mob_seqfish_count.RDS ./processed_data/MOB/ref/internal/scRNA_seqfish_label.RDS
python3 train.py 

# results will be in DSTG_Result/predict_output.csv
# cell type labels can be found in Infor_Data/ST_label/ST_label_1.csv