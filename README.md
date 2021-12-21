<p align="center">
<img src="/Help/OpenIPDM.png" height="110">

<p align="center">
Open-Source Toolbox for Infrastructures Probabilistic Deterioration Modelling 
</p>

OpenIPDM is a MATLAB open-source platform that stands for infrastructures probablstic deterioration model. This software is developed to perform analyses on a network-scale visual inspection data, while accounting for the uncertainty associated with each inspector.
The main application window in OpenIPDM enables assessing the structural deterioration behaviour and the effect of interventions at different levels starting from the structural element level up to the network level.
OpenIPDM also include several toolboxes that facilitate performing verification and validation analyses on visual inspection data, in addition to learning model parameters.
Furthermore, OpenIPDM has the capacity to handle missing data such as, missing interventions or missing structural attributes.

<p align="center">
<img src="/Help/OpenIPDMMain.png" height="500">
    
For tutorials, see: [YouTube channel](https://www.youtube.com/).

## How to cite

*OpenIPDM, an Open-Source Software for Infrastructures Probablstic Deterioration Modelling*

### Prerequisites

- Matlab (version 2020a or higher) installed on Mac OSX or Windows

- The Matlab Statistics and Machine Learning Toolbox is required.

### Installation

1. Download and extract the ZIP file or clone the git repository in your working directory. 
2. The working directory should include the following folders:
    - Scripts
    - Tools
    - Parameters
    - Network Data
    - Figures
    - ExtractedData
    - Help
3. Double-click OpenIPDM.mlapp file to start MATLAB App Designer, and from the top ribbon in App Designer, click *Run*


### Getting started
After starting OpenIPDM, the main user interface will open along with a message box to load the database. Note that the message box will not show up, if pre-processed data already exist in the folder Network Data. 
If you do not see anything except Matlab errors verify your Matlab version, and your Matlab path.

## Input

`OpenIPDM` takes as an input two types of file formats

1. '.csv': this file format is generally considered for the raw database.
2. '.mat': for files containing model paramters and/or pre-processed database.

## Output

`OpenIPDM` generally provides the following outputs:

1. Deterioration state estimates. 
2. Service-life of an intervention.
3. Effect of interventions.

Further details about the outputs can be found in the OpenIPDM documentation manual.


## Remarks


## Built With

* [Matlab](https://www.mathworks.com/products/matlab.html) - Main Development
* [matlab2tikz](https://github.com/matlab2tikz/matlab2tikz) - Figures for LaTeX

## Contributing



## Versioning

To be done...

## Authors

* **Zachary Hamida** - *Initial code and development* - [webpage](http://profs.polymtl.ca/jagoulet/Site/Goulet_web_page_ZHAMIDA.html)
* **Blanche Laurent** - *Development* - [webpage](http://profs.polymtl.ca/jagoulet/Site/Goulet_web_page_BLAURENT.html)

Note: OpenIPDM has been originally developed based on the inspection and interventions database of the Transportation Ministry of Quebec (MTQ).

## License

This project is licensed under the MIT license - see the [LICENSE.txt](LICENSE.txt) file for details

## Acknowledgments

Some parts of the project have greatly benefited from existing work:

* Kevin Murphy (Kalman filter, [Kalman filter toolbox for Matlab](https://www.cs.ubc.ca/~murphyk/Software/Kalman/kalman.html#other))
* Dan Simon (Constrained Kalman filtering, [Constrained Kalman filter](https://academic.csuohio.edu/simond/ConstrKF/))

