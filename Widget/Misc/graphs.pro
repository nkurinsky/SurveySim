pro graphs

gmain = widget_base(title='Simulation Output',/column)
r1 = widget_base(gmain,/row)
r2 = widget_base(gmain,/row)
r3 = widget_base(gmain,/row)

xdim = 500
ydim = 400

lumfunct = widget_draw(r1,xsize=xdim,ysize=ydim)
redshift = widget_draw(r1,xsize=xdim,ysize=ydim)
models = widget_draw(r1,xsize=xdim,ysize=ydim)
dcount1 = widget_draw(r2,xsize=xdim,ysize=ydim)
dcount2 = widget_draw(r2,xsize=xdim,ysize=ydim)
dcount3 = widget_draw(r2,xsize=xdim,ysize=ydim)

rerun = widget_button(r3,uvalue='do',value='Redo Simulation')
close = widget_button(r3,uvalue='close',value='Close')
quit = widget_button(r3,uvalue='quit',value='Quit')

widget_control,gmain,/realize
xmanager,'graphs',gmain,/no_block

set_plot,'x'

dists = mrdfits('output.fits',3,head)

gpts = where(dists.f3 gt 0)
f1 = dists[gpts].f1
f2 = dists[gpts].f2
f3 = dists[gpts].f3
z = dists[gpts].z
m = dists[gpts].m

print,n_elements(gpts)

widget_control,redshift,get_value=index
wset,index
h = histogram(z,nbins=50,locations=xh)
plot,xh,h,psym=10,xrange=[0,10],xstyle=1,xtitle='z',ytitle='dN/dz'

widget_control,models,get_value=index
wset,index

h = histogram(m,nbins=50,locations=xh)
plot,xh,h,psym=10,xstyle=1,xtitle='m',ytitle='dN/dm'

widget_control,dcount1,get_value=index
wset,index

h = histogram(f1,nbins=50,locations=xh)
pts = where(h gt 0)
plot,xh[pts],h[pts],/xlog,/ylog,xrange=[1e-2,1],yrange=[1e-1,1e4],ystyle=1,xstyle=1,xtitle='Flux (Log)',ytitle='dN/dS (Log)'

widget_control,dcount2,get_value=index
wset,index

h = histogram(f2,nbins=50,locations=xh)
pts = where(h gt 0)
plot,xh[pts],h[pts],/xlog,/ylog,xrange=[1e-2,1],yrange=[1e-1,1e4],ystyle=1,xstyle=1,xtitle='Flux (Log)',ytitle='dN/dS (Log)'

widget_control,dcount3,get_value=index
wset,index

h = histogram(f3,nbins=50,locations=xh)
pts = where(h gt 0)
plot,xh[pts],h[pts],/xlog,/ylog,xrange=[1e-2,1],yrange=[1e-1,1e4],ystyle=1,xstyle=1,xtitle='Flux (Log)',ytitle='dN/dS (Log)'

end

pro graphs_event,ev

widget_control,ev.id,get_uvalue=uvalue

case uvalue of
   'do': print,'Restarting Simulation'
   'close': widget_control,ev.top,/destroy
   'quit': begin
      widget_control,ev.top,/destroy
      print,'quitting'
   end
endcase

end
