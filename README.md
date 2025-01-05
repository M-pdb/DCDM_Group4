## DCDM Group 4 - IMPC Data Management Project

# Project Overview
This project is part of the Data Cleaning and Data Management Coursework
The goal is to clean, manage and visualise phenotypic data from the International Mouse Phenotyping Consortium (IMPC)

The project involves:
**Egressing data to be cleaned and collated** contains phenotypic data from knockout mice experiments
**Designing a MySQL database** for storing phenotypic data
**Developing an RShiny Dashboard**  for interactive visualisations of genotype-phenotype association

# Folder Structure
**egressed_data/** – Contains raw `.csv` data retrieved from the TRE
**cleaned_data/** – Processed and cleaned datasets ready for database import
**scripts/** – R scripts for data cleaning, collation, and analysis
**sql/** – SQL scripts for creating and querying the MySQL database 
**rshiny/** – RShiny files to create the interactive dashboard
**database_backup/** – SQL dump files to restore the database

# Clone the Repository
git clone https://github.com/M-pdb/DCDM_Group4.git

# Set up Database
in the sql/ directory, the following code setups the script in MySQL: mysql -u [username] -p < database_setup.sql

# RShiny Dashboard
dashboard includes three visualisations to visualise phenotypic data by selecting specific knockout mice or phenotypes
**Phenotype Scores for Selected Knockout Mice** Displays statistical significance of phenotypes associated with a specific gene knockout
**Phenotype Comparison Across Knockouts** Visualises phenotype scores across all knockout mice for a selected phenotype
**Gene Clusters** Shows clusters of genes with similar phenotype scores to identify related patterns

# Team Members
- Mohammod Maahi Hamza
- Naol Duguma
- Gyumin Oh
- Moritz Ziewer
- Warnakula Perera
