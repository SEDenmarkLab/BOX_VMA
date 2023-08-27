
Run in mollienv

1. Data in duplicate scraped manually from the original raw data into 'data/unprocessed_data_initial_screening.csv'. The original raw spreadsheets from Merck can be found in 'data/raw_data_spreadsheets/'. Note that the numbering in the spreadsheets changed a bit to make the paper more organized. Refer to the numbering in paper figures. The structure is in the upper left corner of each raw data file spreadsheet. Each spreadsheet provides 48 data points. There are 12 total spreadsheets for 576 individual data.

2. Run python homogenize_normalize_data.py

This will create both "processed" and "pretty unprocessed" data sheets. The processed sheets have the correction where the sign of the 4,4'-(S,S) is inverted. Both output data sheets are the average of duplicate, or in some cases triplicate, runs (234 data points). The catalyst_label_key.csv file is used here to correct and assign 4-digit barcodes

3. Run python initial_data_analysis.py > initial-data-analysis.log

This does the Anderson-Darling, F test, etc.

   
All in silico structure have the 4S configuration. No corrections are made in the processed data for the actual
configuration of the catalysts.

The configurations of the catalysts in the experimental data are also tracked in the catalyst_label_key.
You can view the actual catalyst structure configurations in the raw data files and the chemdraw. Note that several of
the structure in Merck's raw spreadsheets are drawn incorrectly, as they were later corrected by XRD. For example 
185_2_1_10 was found to be 4,5-syn and of the 4S configuration by XRD. Everything is corrected in the
chemdraw and the catalyst label key.
