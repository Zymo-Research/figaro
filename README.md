# FIGARO

The beauty of science is to make things simple.  It is with that motto in mind that we developed and present to you FIGARO: An efficient and objective tool for optimizing microbiome rRNA gene trimming parameters.  FIGARO will quickly analyze error rates in a directory of FASTQ files to determine optimal trimming parameters for high-resolution targeted microbiome sequencing pipelines, such as those utilizing [DADA2](https://github.com/benjjneb/dada2 "Github") and [Deblur](https://github.com/biocore/deblur "Github").  The mission of this application is identical to the ZymoBIOMICS mission: increasing reproducibility and standardization of microbiome analysis.

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

####Dockerized version

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

####Command line version

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

End with an example of getting some data out of the system or using it for a little demo

## Running FIGARO

#### Dockerized version

FIGARO is designed to take in data through mounted directories and parameters through environment variables.  As such, the command should appear as follows:

```
docker container run --rm -e AMPLICONLENGTH=450 -e FORWARDPRIMERLENGTH=20 -e REVERSEPRIMERLENGTH=20 -v /path/to/fastqs:/data/input -v /path/to/output:/data/output figaro
```

The user can set several parameters using environment variables passed into the container at runtime. The environment variables that can be passed are as follows:

| Tables        | Are           | Cool  |
| ------------- |:-------------:| -----:|
| col 3 is      | right-aligned | $1600 |
| col 2 is      | centered      |   $12 |
| zebra stripes | are neat      |    $1 |


Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc

