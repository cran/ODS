ODS
============
***
Outcome-dependent sampling (ODS) schemes are cost-effective ways to enhance study efficiency. In ODS designs, one observes the exposure/covariates with a probability that depends on the outcome variable. Popular ODS designs include case-control for binary outcome, case-cohort for time-to-event outcome, and continuous outcome ODS design (Zhou et al. 2002). Because ODS data has biased sampling nature, standard statistical analysis such as linear regression will lead to biases estimates of the population parameters. This package implements four statistical methods related to ODS designs: (1) An empirical likelihood method analyzing the primary continuous outcome with respect to exposure variables in continuous ODS design (Zhou et al., 2002). (2) A partial linear model analyzing the primary outcome in continuous ODS design (Zhou, Qin and Longnecker, 2011). (3) Analyze a secondary outcome in continuous ODS design (Pan et al. 2018). (4) An estimated likelihood method analyzing a secondary outcome in case-cohort data (Pan et al. 2017).

The references are the following: 

Zhou H, Weaver M, Qin J, Longnecker M, Wang M. (2002). A semiparametric empirical likelihood method for data from an outcome‐dependent sampling scheme with a continuous outcome. *Biometrics*, 58(2):413-421.

Zhou H, Qin G, Longnecker M. (2011). A partial linear model in the outcome‐dependent sampling setting to evaluate the effect of prenatal PCB exposure on cognitive function in children. *Biometrics*, 67(3):876-885.

Pan Y, Cai J, Kim S, Zhou H. (2017). Regression analysis for secondary response variable in a case‐cohort study. *Biometrics*. 

Pan Y, Cai J, Longnecker M, Zhou H. (2018). Secondary outcome analysis for data from an outcome‐dependent sampling design. *Statistics in medicine*, 37(15):2321-2337.

### Linear model in ODS ###

We assume that in the population, the primary outcome variable $Y$ follows the linear model:
$$
Y = \beta_{0} + \beta_{1}X + \epsilon
$$
where $X$ are the covariates, and $\epsilon\sim N(0, \sigma^2)$. In continuous ODS design, a simple random sample is taken from the full cohort, then two supplemental samples are taken from tails of the $Y$ distribution, i.e. $(-\infty, \mu_{Y} - a*\sigma_{Y})$ and $(\mu_{Y} + a*\sigma_{Y}, +\infty)$. As ODS data is not a simple random sample of the overall population, naive regression analysis will yield to invalid estimates of the population parameters. Zhou et al. (2002) develops a semiparametric empirical likelihood estimator (MSELE) for conducting inference on the parameters in the linear model. 

Function **odsmle** provides the parameter estimates, and function **se.spmle** calculates the standard error for MSELE estimator.   

### Partial linear model in ODS ###
We assume that in the population, the primary outcome variable $Y$ follows the partial linear model:
$$
E(Y|X,Z)=g(X)+Z^{T}\gamma
$$
where $X$ is the expensive exposure, $Z$ are other covariates. $g(\cdot)$ is an unknown smooth function. Zhou, Qin and Longnecker (2011) considers a penalized spline method to estimate the nonparamatric function $g(\cdot)$ and other regression coefficients $\gamma$ under the ODS sampling scheme. 

Function **Estimate_PLMODS** computes the parameter estimates and standard error in the partial linear model. Function **gcv_ODS** calculates the generalized cross-validation (GCV) for selecting the smoothing parameter. The details can be seen in Zhou, Qin and Longnecker (2011). 

### Secondary analysis in ODS design ###
We assume that in the population, the primary outcome $Y_1$ and the secondary outcome $Y_2$ satisfy the following conditional mean model:
$$
E(Y_1|X,Z)=\beta_0+\beta_1X+\beta_2Z
$$
$$
E(Y_2|X,Z)=\gamma_0+\gamma_1X+\gamma_2Z
$$
Pan et al. (2018) proposed an augmented inverse probability weighted estimating equation to analyze the secondary outcome (parameters: $\gamma_0, \gamma_1, \gamma_2$) for data obtained from the continuous ODS design. Function **secondary_ODS** computes the parameter estimates and standard error for $(\beta, \gamma)$.  

### Secondary analysis of case-cohort data ###
When the primary outcome is survival time, case-cohort design is commonly used to enhance study efficiency. We assume that the primary outcome (survival time) follows the Cox model:
$$
\lambda(t|X,Y_2,Z)=\lambda_0(t)\exp(\gamma_1X+\gamma_2Y_2+\gamma_3Z)
$$
$Y_2$ is a secondary outcome that satisfy the following linear model:
$$
Y_2 = \beta_{0} + \beta_{1}X + \beta_2Z + \epsilon
$$
where $\epsilon\sim N(0, \sigma^2)$. Pan et al. (2017) proposed a nonparametric estimated likelihood approach for analyzing the secondary outcome $Y_2$ when the data is obtained from a case-cohort study. Function **secondary_casecohort** computes the parameter estimates and standard error for $(\beta, \gamma)$. 

### Package installation ###

~~~
install.packages("devtools")
devtools::install_github("Yinghao-Pan/ODS")
~~~
