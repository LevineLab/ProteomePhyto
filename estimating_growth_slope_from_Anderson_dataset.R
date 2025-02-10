# Code by Arianna I Krinos

pacman::p_load(ggplot2,data.table,bbmle,patchwork,tidyr,dplyr)

## Input file
"
sample_data.csv - CSV file with three columns; one column is the growth curve number
e.g., growth curve #1 (column name: `Curve`), 
the next is the temperature (column name: `Temp`), and the third is the growth
rate (column name: `GR`). 
Additional columns can be provided that provide more data detail -
e.g., strain name or experiment number.
"

## Functions ##
nbcurve <- function(x,opt,w,a,b){
  res<-a*exp(b*x)*(1-((x-opt)/(w/2))^2)
  res
}

LL1 <- function (y, x, a, b, w, o){
    N = nbcurve(x=x,a=a,b=b,w=w,opt=o)
    N[N<=0]=0.1
    N=N # eliminate missing data from loglikelihood
    y=y
    return(-sum(dnorm(y,N,log = TRUE),na.rm=TRUE)) # the negative log likelihoods: the order of N and y don't matter)
}

revised_df = data.frame()
derived_traits_all = read.csv("sample_data.csv")
if ("Curve" %in% colnames(derived_traits_all)) {
  group_numbers = unique(derived_traits_all$Curve)
} else {
  print("Column `Curve` is missing.")
  quit(status=1)
}

for (group_no in group_numbers) {
  derived_traits = derived_traits_all %>% dplyr::filter(Curve==group_no)
    
  ## SET UP GUESS VALUES FOR MODEL PARAMETERS ##
  guess_opt = as.numeric(unique((derived_traits %>% 
                            dplyr::filter(GR==max(GR,na.rm=T)))$Temp))
  a = 0.1 # scale param 1
  b = 0.01 # scale param 2
  o_guess = guess_opt # optimum temperature
  w_guess = 15
    
  ## ESTIMATE GROWTH CURVE PARAMETERS USING LOG LIKELIHOOD ## 
  m1 = bbmle::mle2(minuslogl = LL1, start = list(a = a, b = b,
                                            o = o_guess,w = w_guess),
              data = list(y=as.numeric(derived_traits$GR),
                          x=as.numeric(derived_traits$Temp)),
              control=list(maxit=1000000,abstol=1e-20))
    
  a=summary(m1)@coef["a","Estimate"]
  b=summary(m1)@coef["b","Estimate"]
  w=summary(m1)@coef["w","Estimate"]
  opt=summary(m1)@coef["o","Estimate"]
    
  derived_traits = derived_traits %>%
        dplyr::mutate(Estimates = nbcurve(Temp,opt,w,a,b))
  opt_val = optimize(nbcurve,interval=c(-20,40),maximum=TRUE,
                     a=a,
                     b=b,
                     w=w,
                     opt=opt)$maximum

  wid_eq_1 = w/2+opt
  wid_eq_2 = -w/2+opt
  tol_gr=0.01 ## When GR<this, width is not valid.
  test_temp=wid_eq_1
  while ((nbcurve(test_temp,opt,w,a,b) < tol_gr)) {
    test_temp=test_temp+sign(test_temp)*-0.001
  }
  wid_eq_1=test_temp
  test_temp=wid_eq_2
  while ((nbcurve(test_temp,opt,w,a,b) < tol_gr)) {
    test_temp=test_temp+sign(test_temp-wid_eq_1)*-0.001
  }
  wid_eq_2=test_temp

  
  tolerance_val=nbcurve(opt_val,opt,w,a,b)*0.2
  temps_test=seq(from = wid_eq_2+sign(wid_eq_2)*wid_eq_2*0.01, 
                 to = wid_eq_1+sign(wid_eq_2)*wid_eq_2*0.01, by = 0.001)
  window_1_tot = -100
  window_2_tot = -100
  for (temp in temps_test) {
    if ((window_1_tot==-100)&(window_2_tot==-100)&(nbcurve(temp,opt,w,a,b) >= tolerance_val)) {
      window_1_tot = temp
    }else if ((window_2_tot==-100)&(window_1_tot!=-100)&(nbcurve(temp,opt,w,a,b) <= tolerance_val)){
      window_2_tot=temp
    }
  }
  ## if we did not successfully find the window
  if ((window_1_tot==-100)|(window_2_tot==-100)){
    next
  }
  ## if the difference between the two temperatures is unrealistically large
  if (abs(wid_eq_2-wid_eq_1) > 50) {
    next
  }
    
  target = nbcurve(opt_val,opt,w,a,b)*0.2
  print(derived_traits)
  revised_df = rbind(revised_df,
                     data.frame(Revised_opt=opt_val,
                                low_temp_0=window_1_tot,high_temp_0=window_2_tot,
                                slope_20_perc=target,
                                max_GR=nbcurve(opt_val,opt,w,a,b),opt=opt,
                                rmse=sqrt(sum((derived_traits$GR-derived_traits$Estimates)^2)/nrow(derived_traits)),
                                t_opt=opt_val,a=a,b=b,w=w,opt=opt,
                                Group=group_no))
}
revised_df$left_slope=(revised_df$max_GR-revised_df$slope_20_perc)/(revised_df$t_opt-revised_df$low_temp_0)
revised_df$right_slope=(revised_df$max_GR-revised_df$slope_20_perc)/(revised_df$high_temp_0-revised_df$t_opt)

write.csv(revised_df %>% dplyr::select(Group,a,b,w,Revised_opt,
                                       left_slope,right_slope,rmse) %>%
          dplyr::rename(T_opt=Revised_opt),
          "slopes_and_parameters.csv")
