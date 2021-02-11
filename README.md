# FIGARO

*The beauty of science is to make things simple*.  It in this spirit that we developed and present to you FIGARO: An efficient and objective tool for optimizing microbiome rRNA gene trimming parameters.  FIGARO will quickly analyze error rates in a directory of FASTQ files to determine optimal trimming parameters for high-resolution targeted microbiome sequencing pipelines, such as those utilizing [DADA2](https://github.com/benjjneb/dada2 "Github") and [Deblur](https://github.com/biocore/deblur "Github").  The mission of this application is identical to the [ZymoBIOMICS](https://www.zymoresearch.com/pages/zymobiomics-portfolio) mission: increasing the reproducibility and standardization of microbiome analysis.

#### Publication
Please see [FIGARO: An efficient and objective tool for optimizing microbiome rRNA gene trimming parameters](https://www.biorxiv.org/content/10.1101/610394v1 "Preprint version")

## Quick Start Guide

#### Docker
```
git clone https://github.com/Zymo-Research/figaro.git
cd figaro
docker build -t figaro .
docker container run --rm -e AMPLICONLENGTH=[amplicon length] -e FORWARDPRIMERLENGTH=[forward primer length] \
    -e REVERSEPRIMERLENGTH=[reverse primer length] -v /path/to/fastqs:/data/input \
    -v /path/to/output:/data/output figaro
```

#### Command line
```
git clone https://github.com/Zymo-Research/figaro.git
cd figaro
# depending on your system configuration, the following commands
# are either python3/pip3 or python/pip
python3 setup.py bdist_wheel
pip3 install --force-reinstall dist/*.whl
figaro -i /path/to/fastq/directory -o /path/to/output/files -a [amplicon length] \
    -f [forward primer length] -r [reverse primer length]
```

## Getting Started
As with all bioinformatics applications, the key to getting started with FIGARO is to have the right data set up in the right way.  Before getting started you will need:
- A directory of paired-end FASTQ files
    - Both paired-ends should be in the same folder
    - Reads should be from the same sequencing run using the same library preparation method
        - Reads from different sequencing runs or different library prep methods should have their trimming parameters generated separately, as their optimal parameters may be different
    - The Illumina naming standard is currently used to identify forward- and reverse-direction reads. Other naming standards will be supported very soon and can be added easily by a user with a little python knowledge. If you need quick support for your FASTQ naming scheme, please send an email to **mweinstein @t zymoresearch .com**.
- A bit of knowledge about your sequences
    - Expected amplicon size
        - **Amplicon size should only be the length of the targets sequence and should not include primers**
        - For amplicons with some expected biological variation in length, use the longest expected size
    - Forward and reverse primer lengths
    - Desired overlap length
        - Often between 10 and 20 for most applications.
        - DADA2 often recommends 20
        - Remember: longer overlaps means more 3' end must be kept at a cost of decreased overall sequence quality


### Prerequisites

The only prerequisites to running this application is a Python3 interpreter if running outside of Docker on the command line or Docker if using the containerized version.  This application was designed Docker-first and the containerized version is the official recommendation for running this application.


### Installation

This application is designed to be lightweight and simple to use.  The intended use is via its Docker, but can also be run using your Python interpreter (this was written for Python 3.6.7).  For either case, you will need to  clone the FIGARO repository.  This can either be done by downloading the packaged code from Github or running the following commands in the directory where you wish to have FIGARO live:

#### Dockerized version

Clone the FIGARO repository

```
git clone https://github.com/Zymo-Research/figaro.git
```

Descend into the directory containing FIGARO
```
cd figaro
```

Build the container

```
docker build -t figaro .
```

#### Command line version

Clone the FIGARO repository

```
git clone https://github.com/Zymo-Research/figaro.git
```

Descend into the directory containing FIGARO
```
cd figaro
```

Install Python packages
```
pip3 install -r requirements.txt
```

## Running FIGARO

#### Dockerized version

FIGARO is designed to take in data through mounted directories and parameters through environment variables.  As such, the command should appear as follows to run FIGARO on a 450 base amplicon generated using 20 base primers in either direction:

```
docker container run --rm -e AMPLICONLENGTH=450 -e FORWARDPRIMERLENGTH=20 -e REVERSEPRIMERLENGTH=20 -v /path/to/fastqs:/data/input -v /path/to/output:/data/output figaro
```

The user can set several parameters using environment variables passed into the container at runtime. The environment variables that can be passed are as follows:

| Variable        | Type           | Default  | Description |
| --------------- |:--------------:|:--------:|-------------|
AMPLICONLENGTH | integer | **REQUIRED** | The length of the amplified sequence target **not including primers**. User is required to set this.
FORWARDPRIMERLENGTH | integer | **REQUIRED** | The length of the forward primer. User is required to set this.
REVERSEPRIMERLENGTH | integer | **REQUIRED** | The length of the reverse primer. User is required to set this.
OUTPUTFILENAME | string | trimParameters.json | The desired name of the JSON list of trim parameters and their scores
INPUTDIRECTORY | string | /data/input | Directory **inside** the container with FASTQ files. You generally shouldn't have to change this.
OUTPUTDIRECTORY | string | /data/output | Directory **inside** the container for writing output files. You generally shouldn't have to change this.
MINIMUMOVERLAP | integer | 20 | How much you want your paired end sequences to overlap in the middle for merging
SUBSAMPLE | integer | *See description* | What fraction of reads to analyze (1/x) from the FASTQ files. Default value will call a function that sets this based upon the size of the fastq files for a sliding scale.
PERCENTILE | integer | 83 | The percentile to target for read filtering.  The default value of 83 will remove reads that are about 1 standard deviation worse than the average read for that direction in that position. You can generally expect a few percentage points below your percentile value of reads to pass the filtering.
FILENAMINGSTANDARD | string | illumina | Naming convention for files. Currently supporting Illumina and Zymo Services (zymo). Others can be added as requested.

#### Command line version

FIGARO will analyze an entire directory of FASTQ files on the system. Run the following command from FIGARO's directory:

```
python3 figaro.py -i /path/to/fastq/directory -o /path/to/output/files -a 450 -f 20 -r 20
```

The user can set several parameters using environment variables passed into the container at runtime. The environment variables that can be passed are as follows:

| Flag            | Short | Type           | Default  | Description |
|:---------------:|:-----:|:--------------:|:--------:|-------------|
--ampliconLength | -a | integer | **REQUIRED** | The length of the amplified sequence target **not including primers**. User is required to set this.
--forwardPrimerLength | -f | integer | **REQUIRED** | The length of the forward primer. User is required to set this.
--reversePrimerLength | -r | integer | **REQUIRED** | The length of the reverse primer. User is required to set this.
--outputFileName | -n | string | trimParameters.json | The desired name of the JSON list of trim parameters and their scores
--inputDirectory | -i | string | *current working directory* | Directory with FASTQ files.
--outputDirectory | -o | string | *current working directory* | Directory writing output files.
--minimumOverlap | -m | integer | 20 | How much you want your paired end sequences to overlap in the middle for merging
--subsample | -s | integer | *See description* | What fraction of reads to analyze (1/x) from the FASTQ files. Default value will call a function that sets this based upon the size of the fastq files for a sliding scale.
--percentile | -p | integer | 83 | The percentile to target for read filtering.  The default value of 83 will remove reads that are about 1 standard deviation worse than the average read for that direction in that position. You can generally expect a few percentage points below your percentile value of reads to pass the filtering.
--fileNamingStandard | -F | string | illumina | Naming convention for files. Currently supporting Illumina and Zymo Services (zymo). Others can be added as requested.

#### As Python package

FIGARO will analyze an entire directory of FASTQ files on the system. After installation, FIGARO analysis can be run from Python as follows:

```

from figaro import figaro
resultTable, forwardCurve, reverseCurve = figaro.runAnalysis(sequenceFolder, ampliconLength, forwardPrimerLength, reversePrimerLength, minimumOverlap, fileNamingStandard, trimParameterDownsample, trimParameterPercentile)
```

|Parameter        | Type           | Default  | Description |
|:---------------:|:--------------:|:--------:|-------------|
sequenceFolder| string | **REQUIRED** | The folder containing the sequences to analyze. User is required to set this.
minimumCombinedReadLength| integer | **REQUIRED** | The length of the amplified sequence target plus overlap **not including primers**. User is required to set this.
forwardPrimerLength | integer | **REQUIRED** | The length of the forward primer. User is required to set this.
reversePrimerLength | integer | **REQUIRED** | The length of the reverse primer. User is required to set this.
minimumOverlap | integer | 20 | The minimum length of overlap desired for read merging
fileNamingStandard | string | illumina | Naming convention for files. Currently supporting Illumina and Zymo Services (zymo). Others can be added as requested.--outputFileName | -n | string | trimParameters.json | The desired name of the JSON list of trim parameters and their scores
subsample | integer | *See description* | What fraction of reads to analyze (1/x) from the FASTQ files. Default value will call a function that sets this based upon the size of the fastq files for a sliding scale.
percentile | integer | 83 | The percentile to target for read filtering.  The default value of 83 will remove reads that are about 1 standard deviation worse than the average read for that direction in that position. You can generally expect a few percentage points below your percentile value of reads to pass the filtering.

Output from this will be three values in this order: a list of results, ranked by score, an exponential curve object describing the error model for the forward reads, and the same kind of object describing the error model for the reverse reads.




## FIGARO output
Trimming parameter candidates will have the following information:
- Forward trim position
- Reverse trim position
- Forward expected error value
- Reverse expected error value
- Percent read retention
    - Percentage of reads expected to pass filtering with these parameters
- Parameter score
    - Determined by the percentage of reads expected to pass filtering with a penalty for allowed errors that grows exponentially. See our publication for more details on this formula.
    
After a successful run of FIGARO, you can expect the following outputs:
- A list of trimming parameter candidates
    - See above for data associated with trimming parameter sets
    - This list will be output to both the console (STDOUT) and a file.
        - The console will have one set per line, sorted by score, with the highest-scoring set first
        - The file will have a list of parameters in JSON format (spaced for easy reading by humans), sorted by score, with the highest-scoring set first.
    - Plots for expected error values over the course of the amplicon for both forward and reverse directions
        - These will have the exponential regression model included as well as its r-squared
            - This model appears to generally run a very high r-squared value. If the model has an r-squared below 0.95, something is not right.
        - Figures will be in PNG format. I will be open to supporting additional formats in the future if there is adequate demand.
        - These can provide a useful method to monitor sequence qualities over time.

## Contributing

We welcome and encourage contributions to this project from the microbiomics community and will happily accept and acknowledge input (and possibly provide some free kits as a thank you).  We aim to provide a positive and inclusive environment for contributors that is free of any harassment or excessively harsh criticism. Our Golden Rule: *Treat others as you would like to be treated*.

## Versioning

We use a modification of [Semantic Versioning](https://semvar.org) to identify our releases.

Release identifiers will be *major.minor.patch*

Major release: Newly required parameter or other change that is not entirely backwards compatible
Minor release: New optional parameter
Patch release: No changes to parameters

## Authors

- **Michael M. Weinstein** - *Project Lead, Programming and Design* - [michael-weinstein](https://github.com/michael-weinstein)
- **Aishani Prem** - *Testing, Design* - [AishaniPrem](https://github.com/AishaniPrem)
- **Mingda Jin** - *Testing, Code Review* - [jinmingda](https://github.com/jinmingda)
- **Shuiquan Tang** - *Design* - [shuiquantang](https://github.com/shuiquantang)
- **Jeffrey Bhasin** - *Design, Code Review* - [jeffbhasin](https://github.com/jeffbhasin)

See also the list of [contributors](https://github.com/Zymo-Research/figaro/contributors) who participated in this project.

## License

This project is licensed under the GNU GPLv3 License - see the [LICENSE](LICENSE) file for details.
This license restricts the usage of this application for non-open sourced systems. Please contact the authors for questions related to relicensing of this software in non-open sourced systems.

## Acknowledgments

We would like to thank the following, without whom this would not have happened:
* The Python Foundation
* The staff at Zymo Research
* The microbiomics community for making us aware of this need
* Our customers

---------------------------------------------------------------------------------------------------------------------

#### If you like this software, please let us know at info@zymoresearch.com.
#### Please support our continued development of free and open-source microbiomics applications by checking out the latest microbiomics offerings from [ZymoBIOMICS](https://www.zymoresearch.com/pages/zymobiomics-portfolio)
