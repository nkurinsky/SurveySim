function gstring,fvar

tvar=fix(fvar)
if(fvar gt 30000) then tvar=fvar
if(fvar le 9) then svar=string(tvar,format='(i1)')
if((fvar gt 9) and (fvar le 99)) then svar=string(tvar,format='(i2)')
if((fvar gt 99) and (fvar le 999)) then svar=string(tvar,format='(i3)')
if((fvar gt 999) and (fvar le 9999)) then svar=string(tvar,format='(i4)')
if((fvar gt 9999) and (fvar le 99999)) then svar=string(tvar,format='(i5)')
if((fvar gt 99999) and (fvar le 999999)) then svar=string(tvar,format='(i6)')
if((fvar gt 999999) and (fvar le 9999999)) then svar=string(tvar,format='(i7)')
if((fvar gt 9999999) and (fvar le 99999999)) then svar=string(tvar,format='(i8)')

return,svar
end
