# Code by Arianna I Krinos

pacman::p_load(ggplot2,data.table,bbmle,patchwork,tidyr,dplyr)

## Data files
"
1. `lat_phenotypes` - Su says that the previous estimates for the left slopes of these curves were too low compared to Anderson et al. (2021). Need to reestimate these - the independent variable is solely `E_a`
2. `heatmap_simulations` - narrower thermal window, but here we want to derive left slopes according to `E_a`, `ctag_max` and `DIN` using some simulation output
3. Subsample `heatmap_simulations` to see if taking just a subset of points - 5, 10, 15, 20 temperatures - changes significantly the error in determining `T_opt`
"

lat_phenotypes = read.csv("data/lat_phenotypes.csv")
heatmap_simulations = read.csv("data/heatmap_simulations.csv")

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

revised_df=data.frame()
derived_traits_all = heatmap_simulations
unique_vals=derived_traits_all %>% dplyr::distinct(DIN,light,Ea,ctag_max)

for (curr in c(1:nrow(unique_vals))) {
  derived_traits = derived_traits_all %>% dplyr::right_join(unique_vals[curr,])
  strain_curr=paste(as.character(unique_vals[curr,]),collapse="_")
  group_curr="MODEL"
  x = sort(runif(100,-5,40))
    
  ## SET UP GUESS VALUES FOR MODEL PARAMETERS ##
  guess_opt = as.numeric(unique((derived_traits %>% 
                            dplyr::filter(mu==max(mu,na.rm=T)))$temp))
  a = 0.1 # scale param 1
  b = 0.01 # scale param 2
  o_guess = guess_opt # optimum temperature
  w_guess=15
    
  ## ESTIMATE GROWTH CURVE PARAMETERS USING LOG LIKELIHOOD ## 
  m1 = bbmle::mle2(minuslogl = LL1, start = list(a = a, b = b,
                                            o = o_guess,w = w_guess),
              data = list(y=as.numeric(derived_traits$mu),
                          x=as.numeric(derived_traits$temp)),
              control=list(maxit=1000000,abstol=1e-20))
    
  a=summary(m1)@coef["a","Estimate"]
  b=summary(m1)@coef["b","Estimate"]
  w=summary(m1)@coef["w","Estimate"]
  opt=summary(m1)@coef["o","Estimate"]
    
  derived_traits = derived_traits %>%
        dplyr::mutate(Estimates = nbcurve(temp,opt,w,a,b))
  opt_val = optimize(nbcurve,interval=c(-20,40),maximum=TRUE,
                     a=a,
                     b=b,
                     w=w,
                     opt=opt)$maximum

  wid_eq_1 = w/2+opt
  wid_eq_2 = -w/2+opt
  tol_gr=0.01 ## we don't want to say the width is valid when GR<this
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
    
  revised_df = rbind(revised_df,
                     data.frame(DIN_light_Ea_ctagmax=strain_curr,Revised_opt=opt_val,
                                low_temp_0=window_1_tot,high_temp_0=window_2_tot,
                                slope_20_perc=target,
                                max_GR=nbcurve(opt_val,opt,w,a,b),opt=opt,
                                rmse=sqrt(sum((derived_traits$mu-derived_traits$Estimates)^2)/nrow(derived_traits)),
                                t_opt=opt_val,a=a,b=b,w=w,opt=opt))
}
revised_df$left_slope=(revised_df$max_GR-revised_df$slope_20_perc)/(revised_df$t_opt-revised_df$low_temp_0)
revised_df$right_slope=(revised_df$max_GR-revised_df$slope_20_perc)/(revised_df$high_temp_0-revised_df$t_opt)

write.csv(revised_df %>% dplyr::select(DIN_light_Ea_ctagmax,Revised_opt,
                                       a,b,w,opt,left_slope,right_slope),
          "lat_phenotypes_parameters.csv")

## instead, sample a subset of points (n=points_to_sample)
for (points_to_sample in c(10,15,20,25)) {
    for (curr in c(1:nrow(unique_vals))) {
      derived_traits = derived_traits_all %>% 
            dplyr::right_join(unique_vals[curr,]) %>%
            slice_sample(n=points_to_sample)
      strain_curr=paste(as.character(unique_vals[curr,]),collapse="_")
      group_curr="MODEL"
      x = sort(runif(100,-5,40))

      ## SET UP GUESS VALUES FOR MODEL PARAMETERS ##
      guess_opt = as.numeric(unique((derived_traits %>% 
                                dplyr::filter(mu==max(mu,na.rm=T)))$temp))
      a = 0.1 # scale param 1
      b = 0.01 # scale param 2
      o_guess = guess_opt # optimum temperature
      w_guess=15

      ## ESTIMATE GROWTH CURVE PARAMETERS USING LOG LIKELIHOOD ## 
      m1 = bbmle::mle2(minuslogl = LL1, start = list(a = a, b = b,
                                                o = o_guess,w = w_guess),
                  data = list(y=as.numeric(derived_traits$mu),
                              x=as.numeric(derived_traits$temp)),
                  control=list(maxit=1000000,abstol=1e-20))

      a=summary(m1)@coef["a","Estimate"]
      b=summary(m1)@coef["b","Estimate"]
      w=summary(m1)@coef["w","Estimate"]
      opt=summary(m1)@coef["o","Estimate"]
        
        
      derived_traits = derived_traits %>%
        dplyr::mutate(Estimates = nbcurve(temp,opt,w,a,b))
      opt_val = optimize(nbcurve,interval=c(-20,40),maximum=TRUE,
                         a=a,
                         b=b,
                         w=w,
                         opt=opt)$maximum
      wid_eq_1 = w/2+opt
      wid_eq_2 = -w/2+opt
      tol_gr=0.01 ## we don't want to say the width is valid when GR<this
      test_temp=wid_eq_1
      if (!is.na(nbcurve(test_temp,opt,w,a,b))) {
          iter_limit=100
          iters = 0
          while ((nbcurve(test_temp,opt,w,a,b) < tol_gr) & (iters<iter_limit)) {
            test_temp=test_temp+sign(test_temp)*-0.001
            if (is.na(nbcurve(test_temp,opt,w,a,b))) {
                break
            }
            iters = iters + 1
          }
          wid_eq_1=test_temp
          test_temp=wid_eq_2
          iters = 0 
          while ((nbcurve(test_temp,opt,w,a,b) < tol_gr) & (iters<iter_limit)) {
            test_temp=test_temp+sign(test_temp-wid_eq_1)*-0.001
            if (is.na(nbcurve(test_temp,opt,w,a,b))) {
                break
            }
            iters = iters + 1
          }
          wid_eq_2=test_temp
      }
        
      ## if these ended up flipped in order, switch them.
      if (wid_eq_1 < wid_eq_2) {
          wid_eq_2 = wid_eq_1
          wid_eq_1 = test_temp
      }

      tolerance_val=nbcurve(opt_val,opt,w,a,b)*0.8
      temps_test=seq(from = wid_eq_2+sign(wid_eq_2)*wid_eq_2*0.5, 
                     to = wid_eq_1+sign(wid_eq_2)*wid_eq_2*0.5, by = 0.01)
      window_1 = -100
      window_2 = -100
      for (temp in temps_test) {
        if ((window_1<0)&(window_2<0)&(nbcurve(temp,opt,w,a,b) >= tolerance_val)) {
          window_1 = temp
        }else if ((window_2<0)&(window_1>0)&(nbcurve(temp,opt,w,a,b) <= tolerance_val)){
          window_2=temp
        }
      }
        
      ## if we did not successfully find the window
      if ((window_1==-100)|(window_2==-100)){
        next
      }
      ## if the difference between the two temperatures is unrealistically large
      if (abs(wid_eq_2-wid_eq_1) > 50) {
        next
      }

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

      x1=seq(min(temps_test),opt_val, by=0.001)
      x2=seq(opt_val,max(temps_test), by=0.001)
      lowerbound_0.2 <- x1[which(abs(nbcurve(x1,opt,w,a,b)-target)==
                                 min(abs(nbcurve(x1,opt,w,a,b)-target)))]
      upperbound_0.2 <- x2[which(abs(nbcurve(x2,opt,w,a,b)-target)==
                                 min(abs(nbcurve(x2,opt,w,a,b)-target)))]  

      revised_df = rbind(revised_df,
                         data.frame(DIN_light_Ea_ctagmax=strain_curr,Revised_opt=opt_val,
                                    low_temp_0=window_1_tot,high_temp_0=window_2_tot,
                                    slope_20_perc=target,
                                    rmse=sqrt(sum((derived_traits$mu-derived_traits$Estimates)^2)/nrow(derived_traits)),
                                    max_GR=nbcurve(opt_val,opt,w,a,b),opt=opt,
                                    t_opt=opt_val,a=a,b=b,w=w,num_sampled=points_to_sample,
                                    temp_points_included = paste(derived_traits$temp,collapse="_"),
                                    mu_points_included = paste(derived_traits$mu,collapse="_")))
    }
}
revised_df$left_slope=(revised_df$max_GR-revised_df$slope_20_perc)/(revised_df$t_opt-revised_df$low_temp_0)
revised_df$right_slope=(revised_df$max_GR-revised_df$slope_20_perc)/(revised_df$high_temp_0-revised_df$t_opt)

write.csv(revised_df %>% dplyr::filter((left_slope < 100) & (rmse<1)) %>%
          dplyr::select(DIN_light_Ea_ctagmax,Revised_opt,
                        a,b,w,opt,left_slope,right_slope,rmse),
          "heatmap_simulations_parameters.csv")
