# Roiico
Roiico (which stands for Identifying <ins>**Regions of interest**</ins> in imaging data based on <ins>**interval-censored outcomes**</ins>) is a package that performs semiparametric estimation of the group penalized regression for images with irregular boundaries and interval-censored outcomes proposed in Lee et al. (2024+). It requires the installation of the [BPST package](https://first-data-lab.github.io/blogs/docs/BPST.html#few-notes-to-consider) in `R`.

# How to import the Functions #
> install.packages("devtools")<br />
> library(devtools) <br /> 
> library(BPST) <br />
> source_url("https://github.com/lcyjames/Roiico/blob/main/Roiico.R?raw=TRUE")

# Usage #
The package contains 2 functions:
|Functions  | Description|
|------------- | -------------|
RoiicoSIM  | Generate a data set according to the simulation study in Lee et al. (2024+)
RoiicoEST  | Perform the semiparametric estimation methods of Lee et al. (2024+)

|Rds file  | Description|
|------------- | -------------|
VT62.Rds | It records the triangulation grids with 62 triangles stored in an object `VT`, where `VT$V` pertains to the vertices, `VT$Tr` pertains to which vertices are used to form the triangles, and `VT$Z_grid` contains all the Cartesian coordinates of the dots sampled within the boundaries.
VT118.Rds | It shares a similar structure with VT62.Rds, but comprises 118 triangles instead of 62.

<ins>**RoiicoSIM**</ins>

```
RoiicoSIM(seed = NA, n, deg = 2, a, beta, gamma, rho, VT, pattern)
```
This function generates a data set according to the model of the simulation study in Lee et al. (2024+) that takes the following arguments:
>- `n` is the sample size
>- `deg` is the degree of the Bernstein polynomial basis functions over triangulation, set to be 2 by default 
>- `a` is the proportion of triangles associated with the survival outcomes, set to be 0.01, 0.05, or 0.1 in the study
>- `beta` is the coefficient of the non-image covariates `Z`
>- `gamma` is the non-zero coefficient of the image covariates of length `(deg+2)choose(2)` which is common to all selected triangles
>- `rho` is the transformation parameter, set to be 0, 0.5, or 1 in the study
>- `pattern` is a vector containing the (unique) indexes of the triangles to be selected, with length less than or equal to `nrow(VT$Tr)`. A total number of floor(`nrow(VT$Tr)` times `a`) triangles will be selected beginning with the first element of `pattern`.

Example:
```
#This is the setting with sample size 500, number of triangles 62, and transformation parameter = 0
VT   <-readRDS("VT62.rds")
Data <-RoiicoSIM(seed = 1234, n = 500, a = 0.05, beta = c(0.5,-0.5), gamma = c(0.1,0.2,0.3,0.5,0.6,0.4), rho = 0, VT = VT, pattern = c(50,26,25,9,8,7))

head(Data$surv.data)
#  id       Li        Ri DL DI Z1         Z2
#1  1 1.700018 2.6307966  0  1  0 -0.1202458
#2  2 0.000000 0.9046819  1  0  1  0.1724521
#3  3 3.000000       Inf  0  0  1 -0.1113737
#4  4 0.000000 0.2796044  1  0  0 -0.2530162
#5  5 3.000000       Inf  0  0  0 -1.9084928
#6  6 3.000000       Inf  0  0  1  0.5512539

head(Data$selected.true)
#[1] 25 26 50
```

The data structure is as follows:
>- `Data$Y.data` is a matrix with the rows referring to subjects, and the columns referring to the pixel values of the sampled dots within the boundaries
>- `Data$surv.data` is a data frame which contains `Li`: left-endpoint of the interval, `Ri`: right-endpoint of the interval, `DL`: left-censoring indicator, `DI`: interval-censoring indicator; `Z1` and `Z2` are the non-image covariates
>- `Data$selected.true` indexes which triangles are associated with the survival outcome


<ins>**RoiicoEST**</ins>
```
RoiicoEST(Yi, Zi, Li, Ri, DL, DI, rho, VT, deg1 = 2, deg2 = 3, J = 7, tolerance = 10^{-4}, lambda.grid = 10^seq(6,-6,-0.1), TRACE = FALSE)
```
This function performs the semiparametric estimation methods of Lee et al (2024+). The details of the arguments are as follows:
>- `Yi` is a matrix shown above, referring to the structure in `Data$Y.data`, with the number of rows = `n` and the number of columns = total number of dots sampled within the boundaries
>- `Zi` is a matrix of non-image covariates with the number of rows = `n` and the number of columns = number of non-image covariates
>- `Li` is a vector of left-endpoints of the intervals with length `n`, which takes `0` for left-censored subjects; 
>- `Ri` is a vector of right-endpoints of the intervals with length `n`, which takes `Inf` for right-censored subjects
>- `DL` is the left-censoring indicator with length `n`
>- `DI` is the interval-censoring indicator with length `n`
>- `rho` is the transformation parameter, typically set to be 0, 0.5, or 1
>- `deg1` is the degree of the Bernstein polynomial basis functions over triangulation, set to be 2 by default
>- `deg2` is the degree of the I-spline basis functions for the approximation of the cumulative hazard function, set to be 3 by default
>- `J` is the dimension of parameter `omega`; the order of the I-splines equal to `(J-deg2+1)`; `J` must be greater than or equal to `deg2`
>- `tolerance` is the stopping criterion for the EM algorithm, set to 10^{-4} by default
>- `lambda.grid` is a vector of candidate tunning parameter values for penalty term `lambda`
>- `TRACE` is an option for tracking the solution path of the parameter estimation, set to FALSE by default

Example:
```
VT   <-readRDS("VT62.rds")
Data <-RoiicoSIM(seed = 1234, n = 500, a = 0.05, VT = VT, beta = c(0.5,-0.5), gamma = c(0.1,0.2,0.3,0.5,0.6,0.4), rho = 0, pattern = c(50,26,25,9,8,7))
Result<-RoiicoEST(Yi = Data$Y.data, Zi = cbind(Data$surv.data$Z1, Data$surv.data$Z2), Li = Data$surv.data$Li, Ri = Data$surv.data$Ri, 
                  DL = Data$surv.data$DL, DI = Data$surv.data$DI, rho = 0, VT = VT, tolerance = 10^{-3}, lambda.grid = 10^seq(3, 0.5, length.out=100), TRACE = TRUE)
Result
#$selected
#[1] 25 26 50

```
# Contact #
Lee Chun Yin, James <<james-chun-yin.lee@polyu.edu.hk>>

# Reference #
Lee, C. Y., Shi, H., Ma D., Beg M. F., and Cao, J. (2024+). Identification of regions of interest in neuroimaging data with irregular boundaries based on semiparametric transformation models and interval-censored outcomes.
