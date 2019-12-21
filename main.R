rm(list=ls())

#sd=bootstrapping
setwd("C:/Users/Yan Liu/Documents/UMD research/GHGEs research/2 Integrated model code/code for model estimation/library_dcm_may31_2015-1-26/library_dcm")

# read the script file that calibrates the DC4
source("utils/probitUtils.R")
source("utils/create_X.R")
source("utils/models.R")
source("utils/dc_utils.R")
source("models/logit.R")
source("models/dc4.R")

# read your data. if your data is separated by spaces just use read.table

D = read.table("car_ownership_dc4_input_data.csv", header = TRUE, sep = ",")
#D_Primary = D[108:1032,]

# common is a list whose elements are vector of the same length than the
# number of alternatives. The elements of each vector corresponds to
# variables that must have the same coefficients for each alternative
#common = list(c("zero","logsum1","logsum2","logsum3"))
common = list()

# specific is a list that contains as many vectors as there are alternatives
# each vector contains variables that have their own coefficient for the
# corresponding alternative
#specific = list(c(), c("const", "INCOME_LOW", "INCOME_MID", "INCOME_HIGH", "DRVRCNT",   
#                       "HHR_SEX", "HTRESDN_INCOME_LOW", "HTRESDN_INCOME_MID", "HTRESDN_INCOME_HIGH"), c("const", "INCOME_LOW", "INCOME_MID", "INCOME_HIGH",
#                     "DRVRCNT", "HHR_SEX", "HTRESDN_INCOME_LOW", "HTRESDN_INCOME_MID", "HTRESDN_INCOME_HIGH"), c("const", "INCOME_LOW", "INCOME_MID", 
#                      "INCOME_HIGH", "DRVRCNT", "HHR_SEX", "HTRESDN_INCOME_LOW", "HTRESDN_INCOME_MID", "HTRESDN_INCOME_HIGH"))

specific = list(c(), c("const", "DRVRCNT", "HHR_SEX"), c("const", "DRVRCNT", "HHR_SEX"), c("const", "DRVRCNT", "HHR_SEX"))

# reg is the vector that contains the variables of your regression
#reg = c("const", "HHFAMINC", "HHR_SEX", "HTRESDN_1000", "MEAN_COST")
reg = c("const", "HHFAMINC", "MEAN_COST")

# YPrString and YRegString contain the name of the Discrete and the Continuous variable
YPrString = "HHVEHCNT"
YRegString = "MILE_10K"

#build model frame
specPmv = list(
  probit = list(common = common, specific = specific),
  reg = reg,
  YPr = YPrString,
  YReg = YRegString,
  method = "pmvnorm",
  delta = 1e-5,
  reltol = 1e-5,
  SD = "bootstrap",    # none or hessian
  #nboot = 100,
  nboot = 10,
  verbose = TRUE
  #,
  #start = "C:/Users/Huiyun Liu/Desktop/code try/start_dc4_2009_noPTe-8.txt"
  #start = NULL
  #output = "~/Dropbox/dc/results_dc4/PmvPlots_dc4_2009_noPTe-5.pdf",
  #plotrange = 0.8,
  # number of points we use for the plots
  # needs to be relatively big
  #numpoints = 200
)

#estimate start value by logit model directly
#spec=c(list("Y" = specPmv$YPr), specPmv$probit)
#modelFns=logit
#D=D[1:1289,]
#start= model(modelFns, spec, D)

modelFns = dc4

Sys.time()
result_dc4_2009_nb100 = model(modelFns, specPmv, D) 
Sys.time()

