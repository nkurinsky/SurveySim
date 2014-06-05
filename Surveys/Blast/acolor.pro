;Function to ensure all colors are computed correctly
;Current as of 6/22

function acolor,f1,f2,b1,b2

  c = alog10(f1/f2)/alog10(b1/b2)
  return,c

end
