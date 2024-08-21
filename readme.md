<div id="top"></div>


<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->
<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href=https://git.bam.de/tholstei/pepgm/>
    <img src="images/peptonizer.jpg" alt="Logo"  height="300">
  </a>

<h3 align="center">The Peptonizer 2000</h3>

  <p align="center">
    Integrating PepGM and Unipept for probability-based taxonomic inference of metaproteomic samples
    <br />
  </p>
</div>


<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
      </ul>
    </li>
    <li><a href="#input">Input</a></li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
        <li><a href="#preparation">Preparation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

Introducing the Peptonizer2000 - a tool that combines the capabilities of Unipept and PepGM to analyze
metaproteomic mass spectrometry-based samples. Originally designed for taxonomic inference of viral
mass spectrometry-based samples, we've extended PepGM's functionality to analyze metaproteomic samples by
retrieving taxonomic information from the Unipept database.

PepGM is a probabilistic graphical model developed by the eScience group at BAM (Federal Institute for Materials
Research and Testing) that uses belief propagation to infer the taxonomic origin of peptides and taxa in viral samples.
You can learn more about PepGM on our eScience group at BAM (Federal Institute for Materials Research and Testing).
Please refer to our [GitHub](https://github.com/BAMeScience/PepGM) page.

Unipept, on the other hand, is a web-based metaproteomics analysis tool that provides taxonomic information for
identified peptides. To make it work seamlessly with PepGM, we've extended Unipept with new functionalities that
restrict the taxa queried and provide all potential taxonomic origins of the peptides queried. Check out more
information about Unipept [here](https://unipept.ugent.be/).

With the Peptonizer2000, you can look forward to a comprehensive and streamlined workflow that simplifies
the process of identifying peptides and their taxonomic origins in metaproteomic samples.

The Peptonizer2000 workflow is comprised of the following steps:

1. Query all identified peptides, provided by the user in a .tsv file, in the Unipept API,
   and restrict the taxonomic range queried based on any prior knowledge of the sample.
2. Assemble the peptide-taxon associations provided by Unipept into a bipartite graph,
   where peptides and taxa are represented by different nodes, and an edge is drawn between a peptide and a taxon
   if the peptide is part of the taxon's proteome.
3. Transform the bipartite graph into a factor graph using convolution trees and conditional probability table
   factors (CPD).
4. Run the belief propagation algorithm multiple times with different sets of CPD parameters until convergence,
   to obtain posterior probabilities of candidate taxa.
5. Use an empirically deduced metric to determine the ideal graph parameter set.
6. Output the top scoring taxa as a results barchart. The results are also available as comma-separated files
   for further downstream analysis or visualizations.


<div align="center">
    <img src="images/workflow.png" alt="workflow scheme" width="500">
</div>

<br>



<p align="right">(<a href="#top">back to top</a>)</p>

<!-- INPUT -->

## Input

* A .tsv file of your peptides output from any protoemic peptide search method. The first column should be the peptide, the second column it's score attributed by the search engine. An example is provided in test files. <br>
* A config file with your parameters for the peptonizer2000. A more detailed description of the configuration file can be found below. Additionally, an exemplary config file is provided in this repository.

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- GETTING STARTED -->
## Getting Started

### Prerequisites

Make sure you have git installed and clone the repo:
   ```sh
   git clone https://github.com/BAMeScience/Peptonizer2000.git
   ```
The Peptonizer relies on a snakemake workflow developed with snakemake 5.10.0. <br>
Installing snakemake requires mamba.

To install mamba:
  ```sh
conda install -n <your_env> -c conda-forge mamba
  ```

Alternatively, if you do not have conda installed, you can download mamba directly together with miniforge(intructions from the [mamba installation guide](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)):
```sh
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```

To install snakemake:
```sh
conda activate <your_env>
mamba install -c conda-forge -c bioconda -n <your_snakemake_env> snakemake
```
In accordance with the Snakemake recommendations, we suggest to save your sample data 
in `resources` folder. All outputs will be saved in `results`.

Additional dependencies necessary are Java and GCC.

The Peptonizer2000 is tested for Linux OS. <br>

All necessary binaries are autmatically installed using conda.


### Configuration file

The Peptonizer2000 relies on a configuration file in `yaml` format to set up the workflow.
An example configuration file is provided in `config/config.yaml`. <br>
Do not change the config file location.

<details> 
   <details > <summary> Peptonizer parameter </summary>
   <ul>
      <li> DataDir:  Relative path to raw spectra </li>
      <li> ResultsDir: Relative path to results </li>
      <li> ResourcesDir: Relative path to resources </li> 
      <li> ExperimentName: Name of subfolder in results </li>
      <li>TaxaInPlot: # of inferred taxa that appear in the barplot that is created of the results csv</li>
      <li>Alpha: Grid search increments for alpha </li>
      <li>Beta: Grid search increments for beta </li>
      <li>prior: grid search increments for prior </li>
   </ul>
   </details>

   <details > <summary> Sample specific parameter </summary>
   <ul>
      <li> PeptidesAndScores: path to you .tsv file of input peptides</li>
      <li> SampleName: wildcard for spectra file and folder name </li>
   </ul>
   </details>

   </details>

   <details > <summary> UniPept parameter </summary>
   <ul>
       <li>TaxaNumber: # of taxa </li>
       <li>targetTaxa: Comma separated list of taxa compromised in the UniPept query. If querying all of Unipept, use '1'</li>
   </ul> 
   </details>
</details>

### Output files

All Peptonizer2000 output files are saved into the results folder and include the following: <br>

Main results: <br>

- Peptonizer_Results.csv: Table with values ID, score, type (contains all taxids under 'ID' and all probabilities under '
  score' tosterior probabilities of n (default: 15) highest scoring taxa <br>
  <br>

Additional (intermediate): <br>
- Intermediate results folder sorted by their prior value for all possible grid search parameter combinations
- TaxaWeights.csv: csv file of all taxids that had at least one protein map to them and their weight 
- PepGM_graph.graphml: graphml file of the graphical model (without convolution tree factors). Useful to visualize the graph structure and peptide-taxon connections <br>
- paramcheck.png: barplot of the metric used to determine the graphical model parameters for n (default: 15) best performing parameter combinations <br>
- additional .csv files resulting from the clustering of taxa by peptidome
- log files for bug fixing

<p align="right">(<a href="#top">back to top</a>)</p>


## Testing the Peptonizer
<!-- Testing -->

To test the Peptonizer2000 and see if it is set up correctly on your machine, we provide a test file under resources/test_files. This should be dowloaded automatically if you follow the installation instructions above. The test file is a .tsv resulting from the sample S03 of the [CAMPI study](https://www.nature.com/articles/s41467-021-27542-8) searched against a sample specific database using X!Tandem and MS2Rescore. The original file are available through [PRIDE under PXD023217](https://www.ebi.ac.uk/pride/archive/projects/PXD023217/). 

To execute a test run of the Peptonizer2000 using the provided files: 
 
 1. Follow the installation instructions above
 2. In you terminal, go to the folder resources/test_files
 3. execute the following code to move config file to the right directory
 ```sh
 cp ./config.yaml ../../config/
 ```
 4. You need to make some alterations to the provided example config file.
    - input the path to the S03 .tsv file . It should be something like 'path_to_workflow_directory/resources/SampleData/S03_test.tsv'


You should now me all set up to run the Peptonizer2000 on the test files. In your terminal, run
```sh
snakemake --use-conda --cores <n>
````
<n> is the number of cores available on your machine to run this workflow. Make sure your mamba environment, to which you downloaded snakemake, is active.





<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#top">back to top</a>)</p>


<!-- CONTACT -->
## Contact

Tanja Holstein - [@HolsteinTanja](https://twitter.com/HolsteinTanja) - tanja.holstein@ugent.be <br>
Pieter Verschaffelt - pieter.verschaffelt@ugent.be

<p align="right">(<a href="#top">back to top</a>)</p>


<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/BAMeScience/repo_name.svg?style=for-the-badge
[contributors-url]: https://github.com/BAMeScience/repo_name/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/BAMeScience/repo_name.svg?style=for-the-badge
[forks-url]: https://github.com/BAMeScience/repo_name/network/members
[stars-shield]: https://img.shields.io/github/stars/BAMeScience/repo_name.svg?style=for-the-badge
[stars-url]: https://github.com/BAMeScience/repo_name/stargazers
[issues-shield]: https://img.shields.io/github/issues/BAMeScience/repo_name.svg?style=for-the-badge
[issues-url]: https://github.com/BAMeScience/repo_name/issues
[license-shield]: https://img.shields.io/github/license/BAMeScience/repo_name.svg?style=for-the-badge
[license-url]: https://github.com/BAMeScience/repo_name/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/linkedin_username
[product-screenshot]: images/screenshot.png
[Next.js]: https://img.shields.io/badge/next.js-000000?style=for-the-badge&logo=nextdotjs&logoColor=white
[Next-url]: https://nextjs.org/
[React.js]: https://img.shields.io/badge/React-20232A?style=for-the-badge&logo=react&logoColor=61DAFB
[React-url]: https://reactjs.org/
[Vue.js]: https://img.shields.io/badge/Vue.js-35495E?style=for-the-badge&logo=vuedotjs&logoColor=4FC08D
[Vue-url]: https://vuejs.org/
[Angular.io]: https://img.shields.io/badge/Angular-DD0031?style=for-the-badge&logo=angular&logoColor=white
[Angular-url]: https://angular.io/
[Svelte.dev]: https://img.shields.io/badge/Svelte-4A4A55?style=for-the-badge&logo=svelte&logoColor=FF3E00
[Svelte-url]: https://svelte.dev/
[Laravel.com]: https://img.shields.io/badge/Laravel-FF2D20?style=for-the-badge&logo=laravel&logoColor=white
[Laravel-url]: https://laravel.com
[Bootstrap.com]: https://img.shields.io/badge/Bootstrap-563D7C?style=for-the-badge&logo=bootstrap&logoColor=white
[Bootstrap-url]: https://getbootstrap.com
[JQuery.com]: https://img.shields.io/badge/jQuery-0769AD?style=for-the-badge&logo=jquery&logoColor=white
[JQuery-url]: https://jquery.com 