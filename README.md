# Beta-diversity-of-mixtures-over-time
This repo computes a few different beta diversity metrics (eg. Jensen-Shannon divergence) for the following kind of 
experimental setup: two environments are mixed in different proportions, and the resulting community composition is
measured over time. The goal is to understand how the beta diversity of the mixture changes over time. The beta 
diversity comparisons are between each of the two environments and the mixture at each time point.

## Installation
First, clone the repo:
```bash
git clone https://github.com/KoslickiLab/Beta-diversity-of-mixtures-over-time.git
cd Beta-diversity-of-mixtures-over-time
```
Then, install the required packages:
```bash
pip install -r requirements.txt
# or (prefferable) using conda:
conda create -n beta_diversity python=3.8
conda activate beta_diversity
conda install --file requirements.txt
```

## Usage
The main script is `samples_versus_mixture.py` which requires the following arguments:
- `--metadata`: path to the metadata file
- `--count_table`: path to the OTU/ASV count table
- `--output_prefix`: path to the output directory and the prefix desired for output files

There is an optional flag of `--plot` which will generate plots of the beta diversity metrics over time, as well as 
a clustering heatmap of the beta diversity values decorated with metadata information.

An example usage is:
```bash
python .\samples_versus_mixture.py -c example_data\count_data_real.csv -m example_data\metadata_real.csv -o example_data\real -p
```
This will produce the plots, raw beta diversity values of the environments vs. the mixtures, as well as their 
averages (averaged over replicates) in the `example_data` directory. All files will have the prefix `real_`.
## Metadata format
The metadata file should be a CSV file with the following **required** columns:
- `ID`: unique identifier for each sample (matching those in the `count_table`)
- `Environment`: the environment of each sample. **Must** use the names `env1`, `env2`, and `mix`
- `Time`: the time point of each sample (numerical values; eg. 1, 2, 5, 10, etc.)
- `Rep`: the replicate number of each sample (numerical values; eg. 1, 2, 3, etc.)

Please see the [metadata.csv](example_data/metadata.csv) file for an example.

Additional columns only affect the formatting of the optional clustering heatmap plot.

## Count table format
This should be a CSV file with the first column being the unique identifiers for each OTU/ASV, and the remaining columns
representing OTUs/ASVs in each sample. The first row should be the OTU/ASV IDs (arbitrary format), and the remaining 
entries are numeric count data. 

Please see the [count_table.csv](example_data/count_table.csv) file for an example.



