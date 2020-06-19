# AttenGen


## Steps to Recreate the Paper
### Install and Setup
Install [Blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) and [Anaconda or Miniconda](https://docs.anaconda.com/anaconda/install/). Run the commands in `setup.sh` to create your conda environment and install all dependencies.

Unzip data1.zip, data2.zip, data3.zip, and feature_vectors.zip. Merge the contents of data1, data2, and data3 into a new folder called data (data was split into three zip files to circumvent github's max file size).

### Filtering proteins of >30% similarity
Proteins with greater than 30% similarity are often homologous proteins. To prevent any unidentified virulence factors or protective antigens from making it into our dataset of negative examples, we need to filter out anything that could potentially be a homologous protein. To do this, run the file `blast.py` to run a blast search and then run `filter_blast_matches.py` to filter out highly similar proteins from the uniprot database. This will generate our set of non-virulence factors and our set of non-protective antigens needed to train the ML model.

### Training the model
The ML model trains and makes predictions using a 747-dimensional feature vector of chemical descriptors for the protein encoded by each gene. Chemical features are calculated by the propy3 library and include amino acid composition descriptors, Normalized Moreau-Broto autocorrelation descriptors, Geary autocorrelation descriptors, Composition, Transition, Distribution descriptors (CTD), and quasi-sequence order descriptors (QSO). 

To get our data in a format usable by the ML model, run the file `calculate_features.py`, which will take as input the filtered protein files and the victors / protegen data, and will return pickled feature vectors.

Next, run `xgboost_pipeline.py` to split our data into a holdout set and a training/testing set.

Now the model is ready to be trained. To do this, run `xgboost.py`