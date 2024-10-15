# RZ725 Pools
for i in {1..5}; do
    awk -F'\t' '{count[$13]++} END {for (word in count) print word, count[word]}' /scratch/ei56/pa3687/ild-pathwayAnalysis-bal-scRNA/data/01-raw/01-demultiplexed/RZ725_Pool${i}/combined_results_w_combined_assignments.tsv > /scratch/ei56/pa3687/ild-pathwayAnalysis-bal-scRNA/data/01-raw/01-demultiplexed/RZ725_Pool${i}/sample_counts.txt
done

# RZ734 Pools
for i in {1..6}; do
    awk -F'\t' '{count[$13]++} END {for (word in count) print word, count[word]}' /scratch/ei56/pa3687/ild-pathwayAnalysis-bal-scRNA/data/01-raw/01-demultiplexed/RZ734_Pool_${i}/combined_results_w_combined_assignments.tsv > /scratch/ei56/pa3687/ild-pathwayAnalysis-bal-scRNA/data/01-raw/01-demultiplexed/RZ734_Pool_${i}/sample_counts.txt
done
