function beta = calcula_beta(beta_hub, beta_tip, R_hub, R_tip,r)


k = (beta_hub - beta_tip)/ (length(r)-1);

beta = (beta_hub: -k : beta_tip )

end