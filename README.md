# AttenGen
Unzip data1.zip, data2.zip, data3.zip, and feature_vectors.zip. Merge the contents of data1, data2, and data3 into a new folder called data (data was split into three zip files to circumvent github's max file size). Download the uniprot database from [ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz] by opening the link or using cURL/wget, and copy and paste the fasta file into the data folder.

# Steps to Recreate the Paper
Install Blast+ and Anaconda or Miniconda. Run the commands in `setup.sh` to create your conda environment and install all dependencies.