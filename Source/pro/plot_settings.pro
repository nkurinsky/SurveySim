pro plot_settings,plot_type=plot_type
  cleanplot,/silent
  
  if not keyword_set(plot_type) then plot_type='reset'
  
  if strlowcase(plot_type) eq 'ps' then begin
     !p.thick=5
     !x.thick=5
     !y.thick=5
     !p.charthick=5
     !p.charsize=1.5
     cleanplot,/showonly
  endif else if strlowcase(plot_type) eq 'x' then begin
     !p.background=255
     !p.color=0
     cleanplot,/showonly
  endif else begin
     print,"Returning to default plot settings"
  endelse
end
