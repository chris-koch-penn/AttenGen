# AttenGen

AttenGen uses machine learning and genetic algorithms to generate candidate vaccines. It produces live attenuated vaccine candidates by reducing the number of virulence factors and strengthening or maintaining the number of protective antigens while making as few edits to the genome as possible. An improvement in the quantitative fitness of vaccine candidates can be seen with only a few mutations, and large improvements can be seen after tens of mutations. The vaccines can then be synthesized using techniques from [synthetic genomics](https://en.wikipedia.org/wiki/Synthetic_genomics) and recombinant DNA methods and then experimentally tested to see if they display attenuation and an immunogenic response.

"Modulating the abundance of viral gene expression to achieve suitable immunogenicity while limiting virus replication, dissemination, and injury is an essential element of an optimally attenuated virus." - Parks et al., DOI: 10.1128/JVI.75.2.910-920.2001

## Steps to Recreate the Paper

If you want to recreate the paper, follow all of the steps below. If you just want to generate vaccines for any virus or bacteria of interest, follow the [Install and Setup](#install) instructions then skip to the [Generating Vaccine Candidates](#vaxcandidates) section.

### <a name="install"></a>Install and Setup
Install [Blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) and [Anaconda or Miniconda](https://docs.anaconda.com/anaconda/install/). Run the commands in `setup.sh` to create your conda environment and install all dependencies.

Unzip data1.zip, data2.zip, data3.zip, data4.zip, and feature_vectors.zip. Merge the contents of data1, data2, data3, and data4 into a new folder called data (data was split into multiple zip files to circumvent github's max file size).

### Filtering proteins of >30% similarity
Proteins with greater than 30% similarity are often homologous proteins. To prevent any unidentified virulence factors or protective antigens from making it into our dataset of negative examples, we need to filter out anything that could potentially be a homologous protein. To do this, run the command `python blast.py` to run a blast search and then run `python filter_blast_matches.py` to filter out highly similar proteins from the uniprot database. This will generate our set of non-virulence factors and our set of non-protective antigens needed to train the ML model.

### Training the model
The ML model trains and makes predictions using a 747-dimensional feature vector of chemical descriptors for the protein encoded by each gene. Chemical features are calculated by the propy3 library and include amino acid composition descriptors, Normalized Moreau-Broto autocorrelation descriptors, Geary autocorrelation descriptors, Composition, Transition, Distribution descriptors (CTD), and quasi-sequence order descriptors (QSO). 

To get our data in a format usable by the ML model, run the command `python calculate_features.py`, which will take as input the filtered protein files and the victors / protegen data, and will return pickled feature vectors.

Next, run `python xgboost_pipeline.py` to split our data into a holdout set and a training/testing set.

Now the model is ready to be trained. To do this, run `python xgboost_train.py`. By default, the protegen and victors models should be trained on a dedicated GPU. On an RTX 2080, training takes roughly 4 hours. On a CPU, this would probably take over a day or possibly much longer. To train on a CPU, run the script as `python xgboost_train.py cpu`.

### <a name="vaxcandidates"></a>Generating Vaccine Candidates

Run `python genetic_algorithm.py` to start generating Covid-19 vaccine candidates. A file graphing the fitness vs generation will be saved in this directiory and the most fit samples will be saved as a pickled object. To generate vaccines with your own data for a different virus or bacteria, copy a FASTA file containing all DNA coding-sequences of the pathogen into the data folder. The following function can be called from your own script using custom parameters if you want to experiment with different population and generation sizes. Note that the number of generations is the maximum number of mutations any given sample will have. Limiting the number of generations can be used to maintain genetic similarity to the original pathogen.

```python
from genetic_algorithm import run_GA

victors_scores = "./saved_models/victors_xgboost_scores.joblib"
protegen_scores = "./saved_models/protegen_xgboost_scores.joblib"
victors_model_path = "./saved_models/victors_xgboost_model.joblib"
protegen_model_path = "./saved_models/protegen_xgboost_model.joblib"
genome_path = "PATH TO CODING SEQUENCES FOR YOUR PATHOGEN"
run_GA(victors_scores, protegen_scores, victors_model_path, 
protegen_model_path, genome_path, num_generations, pop_size)
```
