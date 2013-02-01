;computes the modified blackbody for the given parameters
;Current as of 6/18

function mod_bb,temp,nu,beta

  h = 6.626d-34
  k = 1.38065d-23
  c = 2.99792d8

  denom = (exp(h * nu / (k * temp)) - 1)*c^2
  result = 2 * h * nu^(3+beta) / denom
  return,result

end
