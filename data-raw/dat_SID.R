## code to prepare `dat_SID` dataset goes here
# simulated ordinal longitudinal data
# helper function
calculate_multinomial_probabilities <- function(lnORs,p0,common=TRUE) {
  # number of ordinal categories
  k <- length(p0)
  if(common) lnORs = rep(lnORs[1],k-1)
  # cumulative logit of p0
  logit0 <- -qlogis(cumsum(p0))[-k]
  # cumulative logit of p1 is shifted by supplied log ORs
  logit1 <- logit0 + lnORs
  # get cumulative probabilities for p1
  cump1 <- plogis(-logit1)
  # get multinomial probabilities for all categories
  p1 <- c(cump1,1)-c(0,cump1)
  return(p1)
}
set.seed(123)
# sample size
N <- 100
# base category probabilities when all covariates are zero
baseprobs <- c(0.02,0.06,0.11,0.12,0.47,0.21,0.01)
# number of visits, unit is week
visits <- paste0("week",c(0,1,2,4,8,12,26)) # baseline is visit = 0
# simulate data from a proportional odds model
# because lower GBS-DS scores are better
# and we model the logit of P(Y >= j)
# so negative coefficients means beneficial effect
# (OR < 1 of a worse GBS-DS score)
# positive time trend towards higher score (self-healing)
time_effects <- -0.05*c(1,2,4,8,12,26)
# covariate effects estimated from SID trial
covs_effects <- c("age"=0.005,"pre_diarrhea"= -0.23)
# time by treatment interaction
time_trt_effects <- -0.05*c(1,2,4,8,12,26)
# correlation structure, the simstudy package generates correlated values from
# the logistic distribution using a standard normal copula-like approach
# with supplied correlation matrix
rho <- 0.6
n_visits <- length(visits)
exponent <- abs(matrix(1:n_visits - 1, nrow = n_visits, ncol = n_visits, byrow = TRUE) -
                  (1:n_visits - 1))
corMatrix <- rho^exponent

# calculate marginal probabilities for trt/control at each visit, without covariate effects
probs_control_allvisits <- t(sapply(time_effects, function(x)
{
  p0 <- calculate_multinomial_probabilities(x,baseprobs)
  return(p0)
}))
probs_control_allvisits <- rbind(baseprobs,probs_control_allvisits)

probs_trt_allvisits <- t(sapply(time_effects+time_trt_effects, function(x)
{
  p1 <- calculate_multinomial_probabilities(x,baseprobs)
  return(p1)
}))
probs_trt_allvisits <- rbind(baseprobs,probs_trt_allvisits)

# define covariate distribution
# 2 prognostic covariates
# age, mean 60 sd 10, from SID-GBS trial
def <- simstudy::defData(varname = "age",dist = "normal",formula = 60,variance = 10^2)
# 40% patients have preceding diarrhea, from SID-GBS trial
def <- simstudy::defData(def,"pre_diarrhea",dist = "binomial",formula = 0.4,variance = 1)
# 1:1 randomization
def <- simstudy::defData(def,varname = "trt",dist = "trtAssign",formula = "1;1")
# relationship specification, for continuous variables, mean-centered
def <- simstudy::defData(def,varname = "z",
                         formula = paste0(
                           ifelse(covs_effects['age']>=0,"+",""),covs_effects['age'],'*age',
                           ifelse(covs_effects['pre_diarrhea']>=0,"+",""),covs_effects['pre_diarrhea'],'*pre_diarrhea'),
                         dist="nonrandom")
data_cov <- simstudy::genData(N,def)
data_cov$age <- round(data_cov$age,0)
data_trt_cov <- data_cov[data_cov$trt==1,]
data_trt <- simstudy::genOrdCat(data_trt_cov,adjVar = "z",
                                baseprobs = probs_trt_allvisits,
                                corMatrix = corMatrix,
                                asFactor = FALSE
)
data_control_cov <- data_cov[data_cov$trt==0,]
data_control <- simstudy::genOrdCat(data_control_cov,adjVar = "z",
                                    baseprobs = probs_control_allvisits,
                                    corMatrix = corMatrix,
                                    asFactor = FALSE)

# change the visit names and order by id
data_wide <- rbind(data_trt,data_control) |>
  dplyr::rename_with(~visits, starts_with("grp")) |>
  dplyr::arrange(id)
# combine and convert to long format
data_long <- reshape2::melt(data_wide,measure.vars = visits,
                            variable.name = "visit",value.name = "GBS_DS") |>
  dplyr::arrange(id)
data_long <- data_long |>
  dplyr::mutate(time = as.numeric(stringr::str_extract(visit,"\\d+$")))
# change variable names to match the SID trial data set
colnames(data_wide)[which(colnames(data_wide)=="id")] <- "patient_ID"
colnames(data_wide)[which(colnames(data_wide)=="trt")] <- "treat_ITT"
colnames(data_wide)[6:12] <- paste0("GBS_DS_week",c(0,1,2,4,8,12,26))
colnames(data_long)[which(colnames(data_long)=="id")] <- "patient_ID"
colnames(data_long)[which(colnames(data_long)=="time")] <- "week"
colnames(data_long)[which(colnames(data_long)=="trt")] <- "treat_ITT"
colnames(data_long)[which(colnames(data_long)=="y")] <- "GBS_DS"

# change GBS-DS range from 1-7 to 0-6 scale
data_wide[,
         paste0("GBS_DS_week",
                c(0,1,2,4,8,12,26))] <- data_wide[,
                                                 paste0("GBS_DS_week",
                                                        c(0,1,2,4,8,12,26))] - 1
data_long$GBS_DS <- data_long$GBS_DS - 1
data_long$GBS_DS_baseline <- rep(data_wide$GBS_DS_week0,
                                as.vector(table(data_long$patient_ID)))
# delete baseline as an outcome
data_long <- data_long[data_long$week>0,]

# delete z
data_wide <- data_wide |>
  dplyr::select(-z)
data_long <- data_long |>
  dplyr::select(-z)

dat_SID <- list("wide_format"=data_wide,
                "long_format"=data_long)



usethis::use_data(dat_SID, overwrite = TRUE)
