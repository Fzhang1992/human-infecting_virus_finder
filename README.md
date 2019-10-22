# human-infecting virus finder
Version: 1.2

Authors: Zheng Zhang, Zena Cai, Zhiying Tan, Taijiao Jiang, Guihua Zhang, Yousong Peng

Maintainer: Zheng Zhang

Base on python 3

Dependencies
-----------
Python packages "pandas","numpy","sklearn","imblearn" are needed to be installed before Installation of train_mod.py and predResult.py
human-infecting virus finder is tested to work under Python 3.6+, the dependency requirements release:

* scipy( >= 1.2)
* numpy( >= 0.16)
* scikit-learn( >= 0.21)
* pandas( >= 0.24)
* imbalanced-learn( >= 0.4.3)
* joblib( >= 0.13)


Installation
-----------
To install "pandas", "numpy" and "sklearn", open terminal and input,

	conda install pandas
	or
	pip install pandas

	conda install numpy
	or
	pip install numpy

	conda install scikit-learn
	or
	pip install -U scikit-learn

To install "imblearn", open terminal and input,

	conda install -c glemaitre imbalanced-learn
	or
	pip install -U imbalanced-learn


Usage
-----------
### decompression the virus_genome.zip

	unzip virus_genome.zip

### example command

	python train_mod.py --infect human_infecting_virus --other Other_viruses --file mod_data --kmer 4


Arguments
-----------
--infect  The file of the human-infecting virus sequences.

--other The file of the other viruses sequences.

--file data The folder of saving generate viral contig k-mer.

--kmer  The length of the k-mer.


Description
-----------
Before running the train_mod.py, please download and extract the virus_genome.zip to the current directory at first. 

The viral genomes containing human-infecting virus and other viruses were split into non-overlapping contigs of 500, 1,000, 3,000, 5,000 and 10,000 nucleotides long.



Usage
-----------
### example command

	python predResult.py --query test --output predict_result --file mod_data --kmer 4 --kjob 1 --bjob 1


Arguments
-----------
--query The file of query sequences.

--output  The name of the predict result file.

--file data The folder of using viral contig k-mer.

--kmer  The length of the k-mer.

--kjob  The number of parallel jobs to run KNeighborsClassifier.

--bjob  The number of parallel jobs to run BalancedBaggingClassifier.


Description
-----------
The function identifies human-infecting virus contig in the input file using the model trained based on the viral genomes.

The rows correspond to sequences, and the columns are from the left to the right, sequence name (name), using model (model), prediction result (label).
When prediction result is "1", mean the model predict this query sequence is human-infecting viral sequence; when prediction result is "0", mean the model predict this query sequence is other viral sequence.

For a query sequence of length L: if L<1kb, the model trained by 500bp sequences is used to predict; if 1kb<=L<3kb, the model trained by 1000bp sequences is used to predict; if 3kb<=L<5kb, the model trained by 3000bp sequences is used to predict; if 5kb<=L<10kb, the model trained by 5000bp sequences is used by predict; if 10kb<=L<15kb, the model trained by 10000bp sequences is used to predict; if L>=15kb, the model trained by the viral genome is used to predict.

The different between kjob and bjob: bjob will faster the kjob, but it will take more memory when using same processors.


Time and memory consumption of the tool
-----------
The time and memory consumption of the software versus the number of sequences inputted in three commonly used operating systems (OS) was listed in Table 1-3 below when testing the software using two threads. More threads can be used to accelerate the prediction with the cost of more memory consumption.

Table 1. The time consumption and memory peak of the software versus the numbers of contigs of various lengths inputted in mac OS (Mojave). The configuration details of the computer: 1.6 GHz Intel Core i5, 8 GB memory.

<p align="center">
  <img src="macOS_Time and memory consumption.png"/>
</p>


Table 2. The time consumption and memory peak of the software versus the numbers of contigs of various lengths inputted in Windows OS (Windows7). The configuration details of the computer: Intel® Core™ is-3740 CPU @ 3.20GHz, 6 GB memory.

<p align="center">
  <img src="Window_Time and memory consumption.png"/>
</p>


Table 3. The time consumption and memory peak of the software versus the numbers of contigs of various lengths inputted in Linux OS (ubuntu 16.04 LTS). The configuration details of the workstation: Intel® Xeon(R) CPU E5-2630 v4 @ 2.20GHz × 40, 125 GB memory.

<p align="center">
  <img src="Linux_Time and memory consumption.png"/>
</p>




