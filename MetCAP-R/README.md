# MetCAP-R

<!-- Put your badges here: -->

<!-- If your software has been certified by MetaboHUB, put the following badge here : -->
[![MetaboHUB Logo](logos/metabohub_logo-20x20.png)![MetaboHUB title](https://img.shields.io/badge/MetaboHub-Software-0066cc?style=flat-square)](https://www.metabohub.fr)

## Metadata

- authors: <jean-marie.savignac@inrae.fr>, <emilien.jamin@inrae.fr>
- creation date: `2025-03-19`
- main usage: Automatize the processing of metabolomic data obtained by XCMS (W4M)

## Description

<!-- NOTE: this section is required -->

The goal of this project is to automatize the different process applied to metabolomic data obtained by the use of XCMS and CAMERA packages from W4M.
A first script is used to curate the datas (variableMetadata), then a second allows us to match the datas with our database, and the last automatize the process of propagation of annotations to unidentified datas.

## Features

<!-- NOTE: this section is required -->

- :sparkles: Curation: Suppression of unusable data for annotation (solvent frontline, NaCL) and CAMERA PcGroup correction
- :rocket: Database creation: Create the internal DB file at R Data format
- :wrench: Shiny_matching: R Shiny application to match the data with the internal DB
- :bar_chart: Automatic_annotation: Detect the specific pattern of parent ion for annotation ([M+H]+ for positive and [M-H]- for negative ion mode) and propagate the annotation following different rules depending of the retention time.


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 
See deployment for notes on how to deploy the project on a live system.

### Prerequisites

<!-- NOTE: this section is required; list what things you need to install the software and how to install them -->

Before you begin, ensure you have met the following requirements:

- **Operating System**: Windows 10+, macOS 10.14+, or Linux
- **Programming Language**: [R] version 4.4.1 or higher
- **Package Manager**: Package versions :

  - readr version R 4.4.2
  - readxl version R 4.4.3
  - haven version R 4.4.3
  - jsonlite version R 4.4.3
  - data.table version R 4.4.3
  - openxlsx version R 4.4.3
  - dplyr version R 4.4.2
  - tidyselect version R 4.4.2
  - tidyverse version R 4.4.2
  - stringr version R 4.4.2
  - tidyr version R 4.4.2
  - shiny version R 4.4.3

- **Other dependencies**: RStudio 2025.05.1+513 "Mariposa Orchid" Release (ab7c1bc795c7dcff8f26215b832a3649a19fc16c, 2025-06-01) for windows
Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) RStudio/2025.05.1+513 Chrome/132.0.6834.210 Electron/34.5.1 Safari/537.36, Quarto 1.6.42

### Installing

<!-- NOTE: this section is required -->

Install R and Rstudio in your computer and open the scripts.

Run scripts by selecting the whole file (CTRL+A) and click run or CTRL+Entry

### Running the tests

<!-- NOTE: this section is required -->

With Rstudio open, open the "Database Creation.R " File and run the script, choose the "Database_Test.xlsx" (or internal database using the file as template) while asked and you will obtain the two databases at Rdata format "DBQTOFNEG_XXXXXX.Rdata" and "DBQTOFPOS_XXXXXX.Rdata" where XXXXXX indicate the date of creation.

Then you can open the "shiny_matching_V6.R" adn run the app on RStudio, a window will open and you can match your file obtained by XCMS with your database and obtain a .txt file similar to the file "Annotated_File_Curation_test_POS.txt"

The next steps depends on the ionisation mode used for analysis, select "POS" or NEG" version of the scripts.

For the next step the "Annotated_File_Curation_test_POS.txt" is used as example.
Open the "Curation_CsvFiles_POS.R " and run the script selecting the annotated file to curate it. You will then obtaine a .xlsx file similar to "Curation_Test_POS_Results.xlsx".

Open the "Automatic_annotation POS.R" and run it, select the "Curation_Test_POS_Results.xlsx" and you will obtain the "data_annoted.xlsx" file.

## Authors

<!-- NOTE: list all authors, maintainers, relevant contributors, ... with their MAIN ROLES and MAIN AFFILIATIONS. -->

- **Jean-Marie Savignac** - *Initial work* - MTH, INRAE, MetaToul
- **Jean-François Martin** - *Developer* - MTH, INRAE, MetaToul
- **Emilien Jamin** - *project management, support* - MTH, INRAE, MetaToul
- **Ludovic Cottret** - *project maintainer, support* - MTH, INRAE, MetaToul
- **Théo Perion** - *support, data provider* - MTH, INRAE, MetaToul
- **Léa Phegnon** - *support, data provider* - MTH, INRAE, MetaToul

<!-- optional: add all contributors in a specific markdown file -->

## License

<!-- NOTE: this section is required; we recommend Cecill License -->

**MetCAP-R** is distributed under the open license CeCILL-2.1 (compatible GNU-GPL).
Please refer to [LICENSE.md](LICENSE.md) file for further details.

## Support &amp; External resources

<!-- NOTE: this section is facultative; customize / remove useless / add relevant / ... items in the list below -->

- :page_facing_up: **CheatSheet** - [MTH CheatSheet](https://todo-link-to-this-resource)
- :scroll: **Wiki**
  - [Wiki for Users](https://todo-link-to-this-resource)
  - [Wiki for Developers](https://todo-link-to-this-resource)
- :question: **Frequently Asked Questions** - [FAQ](https://todo-link-to-this-resource)
- :book: **Documentation** - [Full documentation](https://docs.projectname.com)
- :bug: **Bug Reports** - [GitHub Issues](https://github.com/username/project-name/issues)
- :speech_balloon: **Discussions** - [GitHub Discussions](https://github.com/username/project-name/discussions)
- :email: **Email** - <support@projectname.com>

## Acknowledgments

<!-- NOTE: this section is required -->

- Thanks to [Savignac Jean-Marie](https://forge.inrae.fr/jean-marie.savignac), [Jamin Emilien](https://forge.inrae.fr/emilien.jamin), [Cottret Ludovic](https://forge.inrae.fr/ludovic.cottret) & [Perion Théo](https://forge.inrae.fr/theo.perion) for their valuable contributions
- Built with [R](https://www.r-project.org/) and [RStudio](https://posit.co/download/rstudio-desktop/)
- Special thanks to the open-source community


## Additional Resources

<!-- NOTE: this section is optional -->

- [Project Website](https://projectname.com)
- [API Documentation](https://api.projectname.com/docs)

---
