function spire_errors,n_elems,band

  case band of
     250:begin
        error = 5.8d
     end
     350:begin
        error = 6.3d
     end
     500:begin
        error = 6.8d
     end
  endcase

  dist = randomn(seed,n_elems)
  dist = dist*error
  return,dist

end
