# AttenGen


## Steps to Recreate the Paper
### Install and Setup
Install [Blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) and [Anaconda or Miniconda](https://docs.anaconda.com/anaconda/install/). Run the commands in `setup.sh` to create your conda environment and install all dependencies.

Unzip data1.zip, data2.zip, data3.zip, data4.zip, and feature_vectors.zip. Merge the contents of data1, data2, and data3 into a new folder called data (data was split into multiple zip files to circumvent github's max file size).

### Filtering proteins of >30% similarity
Proteins with greater than 30% similarity are often homologous proteins. To prevent any unidentified virulence factors or protective antigens from making it into our dataset of negative examples, we need to filter out anything that could potentially be a homologous protein. To do this, run the file `python blast.py` to run a blast search and then run `python filter_blast_matches.py` to filter out highly similar proteins from the uniprot database. This will generate our set of non-virulence factors and our set of non-protective antigens needed to train the ML model.

### Training the model
The ML model trains and makes predictions using a 747-dimensional feature vector of chemical descriptors for the protein encoded by each gene. Chemical features are calculated by the propy3 library and include amino acid composition descriptors, Normalized Moreau-Broto autocorrelation descriptors, Geary autocorrelation descriptors, Composition, Transition, Distribution descriptors (CTD), and quasi-sequence order descriptors (QSO). 

To get our data in a format usable by the ML model, run the file `python calculate_features.py`, which will take as input the filtered protein files and the victors / protegen data, and will return pickled feature vectors.

Next, run `python xgboost_pipeline.py` to split our data into a holdout set and a training/testing set.

Now the model is ready to be trained. To do this, run `python xgboost_train.py`. By default, the protegen and victors models should be trained on a dedicated GPU. On an RTX 2080, training takes roughly 4 hours. On a CPU, this would probably take over a day or possibly much longer. To train on a cpu, run the script as `python xgboost_train.py cpu`.

### Evaluate the model

Run `python evaluate_model.py` to print the parameters of the best trained model and the best AUC achieved using 5-fold cross-validation. This script will also test the model on the heldout dataset and produce 4 files containing metrics and confusion matrices. Note that although the Protegen WF1 score was 0.61 for viruses, the precision for the virus class was high at 0.93, meaning most of the error came from saying things were not protective viral antigens when they were rather than mislabeling things as protective viral antigens when they were not. This is the type of error we would rather have. This means the classifier is likely still usable for generating vaccine candidates because if it says something is a protective viral antigen, it likely is (93% precision) so the fitness function of the genetic algorithm is more likely to underestimate the fitness of a sample rather than overestimate it. In short, samples deemed fit by the genetic algorithm are likely even stronger vaccine candidate than the algorithm says they are.

### Generating Vaccine Candidates

Run `python genetic_algorithm.py` to start generating Covid-19 vaccine candidates. To generate vaccines with your own data for your own virus or bacteria, copy a FASTA file containing all DNA coding-sequences into the data folder. The following function can be run from your own script using your own parameters if you want to try different population and generation sizes.


``` 
from genetic_algorithm import run_GA

victors_scores = "./saved_models/victors_xgboost_scores.joblib"
protegen_scores = "./saved_models/protegen_xgboost_scores.joblib"
victors_model_path = "./saved_models/victors_xgboost_model.joblib"
protegen_model_path = "./saved_models/protegen_xgboost_model.joblib"
genome_path = "PATH TO CODING SEQUENCES"
run_GA(victors_scores, protegen_scores, victors_model_path, protegen_model_path, genome_path, num_generations, pop_size)
```