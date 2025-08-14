compute_extinction_prob <- function(waifw, S, beta0, gamma, N, 
                                    tol = 1e-8, max_iter = 1000) {
  S = S/N
  NGM <- beta0 / gamma * waifw %*% diag(S)
  q <- rep(0.5, nrow(NGM))  # start at full extinction
  for (i in 1:max_iter) {
    q_new <- exp(NGM %*% (q - 1))
    if (max(abs(q_new - q)) < tol) break
    q <- q_new
  }
  return(data.frame(extinction_prob = q, outbreak_prob = 1 - q))
}

IC_manual2 = c(1-0.2, 0, 0.2)
names(IC_manual2) = c("S", "I", "R")
IC2 = setup_IC(start_pop = paras["N"], age_classes, c("S", "I", "R"), fert = fert, 
              mort = mort, IC_type = "manual", IC_manual = IC_manual2, dt)
IC2[which(IC2 == max(IC2))] = IC2[which(IC2 == max(IC2))] - 100
IC2[length(IC2)] = 101

eprob = vector("list", length(waifw))
for(i in 1:length(waifw)){
  eprob[[i]] = compute_extinction_prob(
    waifw[[i]], 
    IC2[seq(1, length(IC), 3)], 
    paras["beta0"], paras["gamma"], paras["N"], 
  ) %>% 
    mutate(age_id = 1:dplyr::n(), 
           waifw_id = i) %>% 
    left_join(data.frame(age = age_classes, 
                         bin_width = diff(c(0, age_classes)),
                         age_id = 1:length(age_classes)))
}
