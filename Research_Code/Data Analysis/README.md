# R_ADNI_Analysis
Code used to analyze the [ADNI](https://adni.loni.usc.edu/) data set.


## The three steps of the application.

- **Step 1. Data preparation** 

Removing patients from the analysis, converting to the correct format, creating long data frames. This process is described in the file `1_prepare_data.R`.

- **Step 2. Data conversion to mFData**

As we use the forked [MFPCA](https://github.com/d-gomon/MFPCA) package, we first need to convert the data to the mFData format (see [funData](https://github.com/ClaraHapp/funData)). This process is described in `2_data_to_mFData.R`.

- **Step 3. Application of methods to the resulting data set**

The six described methods are applied to the data in `3_apply_methods_mFData.Rmd`.

## Two additional steps

- **Step 0. Function utility**

Functions used in step 1 to prepare the data. See `0_function_utility.R`.

- **Step 4. Data Summary**

Some summaries used to describe the data set. See `4_Data_Summaries.R`

## Figures

The four resulting figures from Step 3 are attached in the subdirectory `Figures`.

