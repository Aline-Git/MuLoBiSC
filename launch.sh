# ----------------------------------------- #
# First run the Multi distance long-read 
# binning (MuLoBi)
# ----------------------------------------- #

export MuLoBi_dir=./tmp/test_MuLoBi # dummy path for the sake of demo
mkdir -p $MuLoBi_dir  
python3 ./MuLoBi/main.py \
    -i      ./example_test_files/test_contig_sequences.fasta \
    -ci     ./example_test_files/test_short_read_coverage_table.csv \
    -cp     ./example_test_files/test_long_read_coverage_table.csv \
    -b      ./example_test_files/test_busco_table.csv \
    -s      0.035 \
    -p      0.07,0.014,0.004,0.4 \
    -n      0.04779644,22.56103,6 \
    -o      $MuLoBi_dir \
    -m      True 

# ----------------------------------------- #
# Then run the Silhouette coefficient  (SC)
# optimization 
# ----------------------------------------- #

export input_dir=$MuLoBi_dir
export output_dir=./tmp/test_SC_binning_optimization # dummy path for the sake of demo
mkdir -p $output_dir
R --vanilla < ./SC_optimization/SC_binning_optimization.R

