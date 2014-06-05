function blast_errors,n_elems,band

  case band of
     250:begin
        error = 9.5d-2
     end
     350:begin
        error = 8.7d-2
     end
     500:begin
        error = 9.2d-2
     end
  endcase

  dist = randomn(seed,n_elems)
  dist = dist*error
  dist = 1+dist
  return,dist

end
