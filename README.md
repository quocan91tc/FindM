# FindM

Python script helps to correct the starting M position of CDS. Tested on Apple Macbook Pro M1 Pro 10-core CPU with Emihu1 dataset from JGI which combines ~228k genes, it takes about more-or-less 2 hours to finish the job.


## Table of content
1. [Overview](#motivation)
2. [Requirements](#req)
3. [Instruction](#instruction)
4. [Licensing, Authors, and Acknowledgements](#licensing)


### Overview<a name="motivation"></a>
This simple Python script here whill help you, the biologician to correct the M position of CDS in the protein sequence from the JGI Database. I developped this tool for my project of studying the targetting plastidic proteins of Haptophyta. Since many genes in the database were annotated errously, the prediction of Signal Peptide is disappointing. Thus, I have to clean up the data for better prediction performance. I've tested this script on my Apple Macbook Pro M1 Pro 10-core CPU with Emihu1 dataset from JGI which is combined ~228k genes and it takes about more-or-less 2 hours to finish the task.

### Requirements <a name="req"></a>
- Python3.9.20
- [Anaconda](https://anaconda.org/anaconda/conda) or [PyEnv](https://github.com/pyenv/pyenv)
- [Pandas](https://github.com/pandas-dev/pandas)
- [Biopython](https://github.com/biopython/biopython/tree/master)


### Instructions <a name="instruction"></a>:

1. Run the following commands to set up the requirements:
	- To setup the environment:
 		- To setup, it's required either `Anaconda` (MacOS/Linux/Window compatible) or `Pyenv` (Linux/MacOs compatible).
     		- Check if you have already had `Anaconda` installed:
   			`conda --version`.
      		- It should be `conda xx.xx.x`, it shows the version of conda installed.
        	- If it's not, then please install `Anaconda` on your computer by clicking the link in the [Requirements](#req) and follow the instructions. Then after having `Anaconda` installed, you should be able to check the version.
      		- Same as `Anaconda`, you can check wether `Pyenv` is installed or not by using:
        	```
			pyenv --version
         	```
			- You can choose which one suit to you.
      		- Then run these commands to create a virtual environment with `Python3.9.20`. The advantage of using the virtual environment is that you can isolate whatever libraries and Python version installed in this environment to avoid the conflicts with your actual system environment. It's highly recommended to create virtual environment for each project:
        		- Anaconda:
	              	```
			        conda create -n findM python=3.9.20 -c conda-forge
			        conda activate findM
	               	```
           		- Pyenv:
			```
			        pyenv install 3.9.20
			        pyenv versions
			        pyenv virtualenv 3.9.20 myenv
			        pyenv activate myenv
   			```
     		- Then once you finish the task and want to quit, use: `conda deactivate` or `pyenv deactivate`
    - To install the other requirements
      	``` 	   
        pip install -r requirements.txt
       ```

2. To run perfectly the script, you must have the species_all_genes.gff and species_masked_scaffolds.fasta files. They can be downloaded from JGI [here](https://genome.jgi.doe.gov/portal/pages/dynamicOrganismDownload.jsf?organism=haptophyta).

3. The structure of the command:
```
    python findM.py  -s 'species name' -g 'path/to/file.gff' -f 'path/to/file.fasta' -t 'fasta output type'
```
- The option `-t` that let you choose the format of fasta file that you want to get: DNA sequences `'dna'` or Acide Amine `'protein'` sequences or Both `'both'`. Kindly replace those word right after the -t option. 

### Licensing, Authors, Acknowledgements<a name="licensing"></a>
This tool was created by Tran Quoc An Truong ([LinkedIn](https://www.linkedin.com/in/tran-quoc-an-truong/)).
