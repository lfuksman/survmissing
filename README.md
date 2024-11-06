# survmissing

The codes in the repository are algorithms for data-driven nonparametric estimators of marginal and conditional density and hazard rate for censored data with different missing mechanisms. In survival analysis lifetimes of interest are often not observed directly due to right censoring. For instance, in a cancer study where the lifetime of interest is from the time of surgery to cancer relapse, that time is right censored if the patient dies from unrelated conditions before the relapse happens. Moreover, 
because this type of event data comes from real-world data collection, missing data is a common complication. Estimation of density and hazard rate, taking the incompleteness of the data into account, is of interest in the survival analysis field as it can allow experts to have more accurate insights into the underlying risk.
Using Rubin’s terminology, there are three types of missing: missing completely at random (MCAR), missing at random (MAR) and missing not at random (MNAR). The various R codes in the repository are supplements to the following papers: 

1.	Efromovich, S.* and Fuksman, L.* Study of imputation procedures for nonparametric density estimation based on missing censored lifetimes. Computational Statistics & Data Analysis, 2024, Volume 198. https://doi.org/10.1016/j.csda.2024.107994.
    - paper_1_code: univariate density estimation under MCAR missing

2.	Efromovich, S.* and Fuksman, L.* Nonparametric estimation for missing right censored lifetimes. Submitted, in review.
    - paper_2_code: marginal and conditional density and hazard rate estimation under MAR missing	 

3.	Efromovich, S.* and Fuksman, L.* Missing in Survival Analysis. In: Ansari, J., et al. Combining, Modelling and Analyzing Imprecision, Randomness and Dependence. SMPS 2024. Advances in Intelligent Systems and Computing, 2024, vol 1458. Springer, Cham. https://doi.org/10.1007/978-3-031-65993-5_16.
	- paper_3_code: bivariate density estimation under MCAR missing 

