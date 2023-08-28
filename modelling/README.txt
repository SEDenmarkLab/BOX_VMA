
This is all done in mollienv environment (see insilico libraries for dependencies)

The data in preprocess_x is the same in the descriptors folder. Will need unzipped to run. Y data is from ee_data folder. Aggregated features after substrate-dependent RFE are included in aggregated_features. The full library aggregated features are also in here.


1. in preprocess_y, run make_LOCO_splits.py

this converts corrected data to ddg (kcal/mol), transforms the Y data with Yeo-Johnson, makes the LOCO splits and saves files, which is important for parallelized workflows later.

2. run RFE_LOCO_ridge.py

this does the substrate-agnostic RFE

3. run model_RFE_splits_subs_dependent.py

this does substrate-dependent RFE

4. run aggregate_features.py

this does the voting scheme and aggregates the features. Note that the whole insilico library features are too large to add to git repo, so I have only included the reduced feature space.

5. run test_aggregated_features.py

this tests our modelling architectures on the aggregated feature spaces. It also saves pickled models for each one. Helper functions to do random and one-hot controls are also in this file.

6. run prediction.py

this runs our predictions on the nn's for out-of-sample model test