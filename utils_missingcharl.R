# estimating the parameters for the Charlson missingness model (mm):

# inputs: vector V of n_visits, vector Yobs of 1 for 'code', 0 for no code (0 or NA
# in original data)

# theta = (p,phi) 

estimate_mm <- function(Y_obs,V){
  # initialize
  theta_old <- init_mm(Y_obs,V,phi=FALSE)
  # optimization parameters
  eps <- 1e-4
  Kmax <- 100
  step_scale <- 1
  # loop
  for(kk in 1:Kmax){
    theta_new <- newton_mm(theta_old,Y_obs,V,step_scale)
    if(mean((theta_old - theta_new)^2) < eps){
      break
    }
    theta_old <- theta_new
  }
  # return
  return(theta_new)
}

init_mm <- function(Y_obs,V,phi=FALSE){
  p_init <- mean(Y_obs[V > 200])
  if(phi){
    phi_init <- mean(na.omit(sapply(1:200,function(k){1-(1 - (mean(Y_obs[V==k])/p_init))^(1/k)})))
  }
  else{
    phi_init <- 0.02
  }
  return(c(p_init,phi_init))
}

newton_mm <- function(theta,Y_obs,V,step_scale=1){
  p <- theta[1]
  phi <- theta[2]
  # sub components of gradient and Hessian
  c1 <- 1 - (1-phi)^V
  c2 <- V*((1-phi)^(V-1))
  c3 <- V*(V-1)*((1-phi)^(V-2))
  # components of gradient
  g1 <- sum((1/p)*(Y_obs==1)) + sum(((-c1)/(1-p*c1))[(Y_obs==0)])
  g2 <- sum((c2/c1)[((Y_obs==1) & (V > 0))]) + sum(((-p*c2)/(1-p*c1))[((Y_obs==0) & (V > 0))])
  grad <- c(g1,g2)
  # components of Hessian
  H11 <- sum(((-1)/(p^2))*(Y_obs==1)) + sum(((-(c1)^2)/((1-p*c1)^2))*(Y_obs==0))
  H12 <- sum((c2/((1-(p^2)*c1)^2))*(Y_obs==0))
  H22 <- sum(((-c3*c1 - c2^2)/(c1^2))[((Y_obs==1) & (V > 0))]) + sum((((-c3*((1/p) - c1)) - c2^2)/(((1/p) - c1)^2))[((Y_obs==0) & (V > 0))])
  hess <- matrix(c(H11,H12,H12,H22),2,2)
  # newton step
  theta_new <- theta - step_scale*solve(hess,grad)
  return(theta_new)
}

# impute using theta estimated with the missingness model above
impute_mm <- function(Y_obs,V,theta){
  # copy
  Y_new <- Y_obs
  # calculate expectations
  Q <- (1-theta[2])^V
  Y_imp <- (theta[1]*Q) / (1-theta[1]*(1-Q))
  # fill imputed values
  Y_new[Y_new==0] <- Y_imp[Y_new==0]
  # return imputed indicator
  return(Y_new)
}

# define dictionaries to process the new charlson indicators
# 15 diseases (no 'cancer' or 'met_car'), plus 'gerd' = 16 new indicators

# order from HOSEA data dictionary Sample page
charl_names <- c('chf','ctd','dem','diab_c',
                 'gerd','hiv','mld','msld',
                 'para','rd','cd','copd',
                 'diab_nc','mi','pud','pvd')

# translated from Jenny's SQL code 
charl_checks <- list(
  'chf'=list(
    'icd9'=function(s){
      as.integer(pmax(substr(s,1,6)=='398.91',
                      substr(s,1,6)=='402.01',
                      substr(s,1,6)=='402.11',
                      substr(s,1,6)=='402.91',
                      substr(s,1,6)=='404.01',
                      substr(s,1,6)=='404.03',
                      substr(s,1,6)=='404.11',
                      substr(s,1,6)=='404.13',
                      substr(s,1,6)=='404.91',
                      substr(s,1,6)=='404.93',
                      substr(s,1,5)=='425.4',
                      substr(s,1,5)=='425.5',
                      substr(s,1,5)=='425.7',
                      substr(s,1,5)=='425.8',
                      substr(s,1,5)=='425.9',
                      substr(s,1,4)=='428.'))
    },
    'icd10'=function(s){
      as.integer(pmax(substr(s,1,5)=='I09.9',
                      substr(s,1,5)=='I11.0',
                      substr(s,1,5)=='I13.0',
                      substr(s,1,5)=='I13.2',
                      substr(s,1,5)=='I25.5',
                      substr(s,1,5)=='I42.0',
                      substr(s,1,5)=='P29.0',
                      substr(s,1,3)=='I43',
                      substr(s,1,3)=='I50',
                      substr(s,1,5)=='I42.5',
                      substr(s,1,5)=='I42.6',
                      substr(s,1,5)=='I42.7',
                      substr(s,1,5)=='I42.8',
                      substr(s,1,5)=='I42.9'))
    }),
  'ctd'=list(
    'icd9'=function(s){
      as.integer(pmax(substr(s,1,5)=='446.5',
                      substr(s,1,5)=='710.0',
                      substr(s,1,5)=='710.1',
                      substr(s,1,5)=='710.2',
                      substr(s,1,5)=='710.3',
                      substr(s,1,5)=='710.4',
                      substr(s,1,5)=='714.0',
                      substr(s,1,5)=='714.1',
                      substr(s,1,5)=='714.2',
                      substr(s,1,5)=='714.8',
                      substr(s,1,3)=='725'))
    },
    'icd10'=function(s){
      as.integer(pmax(substr(s,1,3)=='M05',
                      substr(s,1,3)=='M06',
                      substr(s,1,5)=='M31.5',
                      substr(s,1,5)=='M35.1',
                      substr(s,1,5)=='M35.3',
                      substr(s,1,5)=='M36.0',
                      substr(s,1,3)=='M32',
                      substr(s,1,3)=='M33',
                      substr(s,1,3)=='M34'))
    }),
  'dem'=list(
    'icd9'=function(s){
      as.integer(pmax(substr(s,1,4)=='290.',
                      substr(s,1,5)=='294.1',
                      substr(s,1,5)=='331.2'))
    },
    'icd10'=function(s){
      as.integer(pmax(substr(s,1,3)=='F00',
                      substr(s,1,3)=='F01',
                      substr(s,1,3)=='F02',
                      substr(s,1,3)=='F03',
                      substr(s,1,5)=='F05.1',
                      substr(s,1,3)=='G30',
                      substr(s,1,5)=='G31.1'))
    }),
  'diab_c'=list(
    'icd9'=function(s){
      as.integer(pmax(substr(s,1,5)=='250.4',
                      substr(s,1,5)=='250.5',
                      substr(s,1,5)=='250.6',
                      substr(s,1,5)=='250.7'))
    },
    'icd10'=function(s){
      as.integer(substr(s,1,5) %in% paste0('E1',apply(expand.grid(c(0:4),c(2:5,7)),1,paste0,collapse='.')))
    }),
  'gerd'=list(
    'icd9'=function(s){
      as.integer(pmax(substr(s,1,6)=='530.11',
                      substr(s,1,6)=='530.81'))
    },
    'icd10'=function(s){
      as.integer(substr(s,1,3)=='K21')
    }),
  'hiv'=list(
    'icd9'=function(s){
      as.integer(pmax(substr(s,1,4)=='042.',
                      substr(s,1,4)=='043.',
                      substr(s,1,4)=='044.'))
    },
    'icd10'=function(s){
      as.integer(pmax(substr(s,1,3)=='B20',
                      substr(s,1,3)=='B21',
                      substr(s,1,3)=='B22',
                      substr(s,1,3)=='B24'))
    }),
  'mld'=list(
    'icd9'=function(s){
      as.integer(pmax(substr(s,1,6)=='070.22',
                      substr(s,1,6)=='070.23',
                      substr(s,1,6)=='070.32',
                      substr(s,1,6)=='070.33',
                      substr(s,1,6)=='070.44',
                      substr(s,1,6)=='070.54',
                      substr(s,1,5)=='070.6',
                      substr(s,1,5)=='070.9',
                      substr(s,1,4)=='570.',
                      substr(s,1,4)=='571.',
                      substr(s,1,5)=='573.3',
                      substr(s,1,5)=='573.4',
                      substr(s,1,5)=='573.8',
                      substr(s,1,5)=='573.9',
                      substr(s,1,5)=='V427.'))
    },
    'icd10'=function(s){
      as.integer(pmax(substr(s,1,3)=='B18',
                      substr(s,1,5)=='K70.0',
                      substr(s,1,5)=='K70.1',
                      substr(s,1,5)=='K70.2',
                      substr(s,1,5)=='K70.3',
                      substr(s,1,5)=='K70.9',
                      substr(s,1,5)=='K71.7',
                      substr(s,1,5)=='K76.0',
                      substr(s,1,5)=='K76.8',
                      substr(s,1,5)=='K76.9',
                      substr(s,1,5)=='Z94.4',
                      substr(s,1,5)=='K71.3',
                      substr(s,1,5)=='K71.4',
                      substr(s,1,5)=='K71.5',
                      substr(s,1,3)=='K73',
                      substr(s,1,3)=='K74',
                      substr(s,1,5)=='K76.2',
                      substr(s,1,5)=='K76.3',
                      substr(s,1,5)=='K76.4'))
    }),
  'msld'=list(
    'icd9'=function(s){
      as.integer(pmax(substr(s,1,5)=='456.0',
                      substr(s,1,5)=='456.1',
                      substr(s,1,5)=='456.2',
                      substr(s,1,5)=='572.2',
                      substr(s,1,5)=='572.3',
                      substr(s,1,5)=='572.4',
                      substr(s,1,5)=='572.8'))
    },
    'icd10'=function(s){
      as.integer(pmax(substr(s,1,5)=='I85.0',
                      substr(s,1,5)=='I85.9',
                      substr(s,1,5)=='I86.4',
                      substr(s,1,5)=='I98.2',
                      substr(s,1,5)=='K70.4',
                      substr(s,1,5)=='K71.1',
                      substr(s,1,5)=='K72.1',
                      substr(s,1,5)=='K72.9',
                      substr(s,1,5)=='K76.5',
                      substr(s,1,5)=='K76.6',
                      substr(s,1,5)=='K76.7'))
    }),
  'para'=list(
    'icd9'=function(s){
      as.integer(pmax(substr(s,1,5)=='334.1',
                      substr(s,1,4)=='342.',
                      substr(s,1,4)=='343.',
                      substr(s,1,5)=='344.0',
                      substr(s,1,5)=='344.1',
                      substr(s,1,5)=='344.2',
                      substr(s,1,5)=='344.3',
                      substr(s,1,5)=='344.4',
                      substr(s,1,5)=='344.5',
                      substr(s,1,5)=='344.6',
                      substr(s,1,5)=='344.9'))
    },
    'icd10'=function(s){
      as.integer(pmax(substr(s,1,5)=='G04.1',
                      substr(s,1,5)=='G11.4',
                      substr(s,1,5)=='G80.1',
                      substr(s,1,5)=='G80.2',
                      substr(s,1,3)=='G81',
                      substr(s,1,3)=='G82',
                      substr(s,1,5)=='G83.0',
                      substr(s,1,5)=='G83.1',
                      substr(s,1,5)=='G83.2',
                      substr(s,1,5)=='G83.3',
                      substr(s,1,5)=='G83.4',
                      substr(s,1,5)=='G83.9'))
    }),
  'rd'=list(
    'icd9'=function(s){
      as.integer(pmax(substr(s,1,6)=='403.01',
                      substr(s,1,6)=='403.11',
                      substr(s,1,6)=='403.91',
                      substr(s,1,6)=='404.02',
                      substr(s,1,6)=='404.03',
                      substr(s,1,6)=='404.12',
                      substr(s,1,6)=='404.13',
                      substr(s,1,6)=='404.92',
                      substr(s,1,6)=='404.93',
                      substr(s,1,4)=='582.',
                      substr(s,1,5)=='583.0',
                      substr(s,1,5)=='583.1',
                      substr(s,1,5)=='583.2',
                      substr(s,1,5)=='583.4',
                      substr(s,1,5)=='583.6',
                      substr(s,1,5)=='583.7',
                      substr(s,1,4)=='585.',
                      substr(s,1,4)=='586.',
                      substr(s,1,5)=='588.0',
                      substr(s,1,5)=='V420.',
                      substr(s,1,5)=='V451.',
                      substr(s,1,4)=='V56.'))
    },
    'icd10'=function(s){
      as.integer(pmax(substr(s,1,5) %in% paste0('N0',apply(expand.grid(c(3,5),c(2:7)),1,paste0,collapse='.')),
                      substr(s,1,5)=='I12.0',
                      substr(s,1,5)=='I13.0',
                      substr(s,1,3)=='N18',
                      substr(s,1,3)=='N19',
                      substr(s,1,5)=='N25.0',
                      substr(s,1,5)=='Z49.0',
                      substr(s,1,5)=='Z49.1',
                      substr(s,1,5)=='Z49.2',
                      substr(s,1,5)=='Z94.0',
                      substr(s,1,5)=='Z99.2'))
    }),
  'cd'=list(
    'icd9'=function(s){
      as.integer(pmax(substr(s,1,6)=='362.34',
                      substr(s,1,4)=='430.',
                      substr(s,1,4)=='431.',
                      substr(s,1,4)=='432.',
                      substr(s,1,4)=='433.',
                      substr(s,1,4)=='434.',
                      substr(s,1,4)=='435.',
                      substr(s,1,4)=='436.',
                      substr(s,1,4)=='437.',
                      substr(s,1,4)=='438.'))
    },
    'icd10'=function(s){
      as.integer(pmax(substr(s,1,3) %in% paste0('I6',c(0:9)),
                      substr(s,1,3)=='G45',
                      substr(s,1,3)=='G46',
                      substr(s,1,5)=='H34.0'))
    }),
  'copd'=list(
    'icd9'=function(s){
      as.integer(pmax(substr(s,1,5)=='416.8',
                      substr(s,1,5)=='416.9',
                      substr(s,1,4)=='490.',
                      substr(s,1,4)=='491.',
                      substr(s,1,4)=='492.',
                      substr(s,1,4)=='493.',
                      substr(s,1,4)=='494.',
                      substr(s,1,4)=='495.',
                      substr(s,1,4)=='496.',
                      substr(s,1,4)=='500.',
                      substr(s,1,4)=='501.',
                      substr(s,1,4)=='502.',
                      substr(s,1,4)=='503.',
                      substr(s,1,4)=='504.',
                      substr(s,1,4)=='505.',
                      substr(s,1,5)=='506.4',
                      substr(s,1,5)=='508.1',
                      substr(s,1,5)=='508.8'))
    },
    'icd10'=function(s){
      as.integer(pmax(substr(s,1,3) %in% paste0('J',apply(expand.grid(c(4,6),c(0:7)),1,paste0,collapse='')),
                      substr(s,1,5)=='I27.8',
                      substr(s,1,5)=='I27.9',
                      substr(s,1,5)=='J70.1',
                      substr(s,1,5)=='J70.3'))
    }),
  'diab_nc'=list(
    'icd9'=function(s){
      as.integer(pmax(substr(s,1,5)=='250.0',
                      substr(s,1,5)=='250.1',
                      substr(s,1,5)=='250.2',
                      substr(s,1,5)=='250.3',
                      substr(s,1,5)=='250.8',
                      substr(s,1,5)=='250.9'))
    },
    'icd10'=function(s){
      as.integer(substr(s,1,5) %in% paste0('E1',apply(expand.grid(c(0:4),c(0,1,6,8,9)),1,paste0,collapse='.')))
    }),
  'mi'=list(
    'icd9'=function(s){
      as.integer(pmax(substr(s,1,4)=='410.',
                      substr(s,1,4)=='412.'))
    },
    'icd10'=function(s){
      as.integer(pmax(substr(s,1,3)=='I21',
                      substr(s,1,3)=='I22',
                      substr(s,1,5)=='I25.2'))
    }),
  'pud'=list(
    'icd9'=function(s){
      as.integer(pmax(substr(s,1,4)=='531.',
                      substr(s,1,4)=='532.',
                      substr(s,1,4)=='533.',
                      substr(s,1,4)=='534.'))
    },
    'icd10'=function(s){
      as.integer(pmax(substr(s,1,3)=='K25',
                      substr(s,1,3)=='K26',
                      substr(s,1,3)=='K27',
                      substr(s,1,3)=='K28'))
    }),
  'pvd'=list(
    'icd9'=function(s){
      as.integer(pmax(substr(s,1,5)=='093.0',
                      substr(s,1,5)=='437.3',
                      substr(s,1,4)=='440.',
                      substr(s,1,4)=='441.',
                      substr(s,1,5)=='443.1.',
                      substr(s,1,5)=='443.2',
                      substr(s,1,5)=='443.8',
                      substr(s,1,5)=='443.9',
                      substr(s,1,5)=='447.1',
                      substr(s,1,5)=='557.1',
                      substr(s,1,5)=='557.9',
                      substr(s,1,5)=='V434.'))
    },
    'icd10'=function(s){
      as.integer(pmax(substr(s,1,5)=='I73.1',
                      substr(s,1,5)=='I73.8',
                      substr(s,1,5)=='I73.9',
                      substr(s,1,5)=='I77.1',
                      substr(s,1,5)=='I79.0',
                      substr(s,1,5)=='I79.2',
                      substr(s,1,5)=='K55.1',
                      substr(s,1,5)=='K55.8',
                      substr(s,1,5)=='K55.9',
                      substr(s,1,5)=='Z95.8',
                      substr(s,1,5)=='Z95.9',
                      substr(s,1,3)=='I70',
                      substr(s,1,3)=='I71'))
    }))


