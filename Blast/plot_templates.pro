;Makes model templates, modified blackbody B(v)v^beta
;creates a fits file as output
;Current as of 6/22

pro plot_templates

  !p.thick=5
  !x.thick=5
  !y.thick=5
  !y.margin = [5,2]
  !x.margin = [7,3]
  !p.charthick=4
  !p.charsize=1.5
  !p.multi=[0,3,4]

  set_plot,'ps'
  device,filename='plots/template.eps',xsize=16,ysize=12,/inches,/times,/color,/encapsulated
  wave = (findgen(991)+10)*1.d4
  true_wave = wave*1d-10
  freq = 2.9979d8/true_wave
  plot_wave = wave*1d-4

  beta_high = 2.d0
  beta_low = 1.d0
  beta_step = 0.1d0
  temp_high = 100.d0
  temp_low = 10.d0
  temp_step = 10.d0

  t_it = 0
  betas = textoidl('\beta:')
  bcolors = 0

  for i=beta_low,beta_high+beta_step,beta_step do begin
     betas = [betas,strmid(strcompress(string(i),/remove_all),0,3)]
     bcolors = [bcolors,230.0*i-210.d0]
  endfor
  
  loadct,0,/silent
  for temp=temp_low,temp_high,temp_step do begin
     print,temp
     lowest = mod_bb(temp,2.9979e12,beta_low)
     highest = max(mod_bb(temp,freq,beta_high))
     multiplot
     switch t_it of
        0:
        3: 
        6: begin
           plot,/nodata,xstyle=1,ystyle=1,[1000,10],[1d-6,1d14],/xlog,/ylog,ytitle=textoidl('B_{\nu}(T[K])\nu^{\beta}')
           break
        end
        1:
        2: 
        4: 
        5: begin
           plot,/nodata,xstyle=1,ystyle=1,[1000,10],[1d-6,1d14],/xlog,/ylog
           break
        end
        7:
        8: begin
            plot,/nodata,xstyle=1,ystyle=1,[1000,10],[1d-6,1d14],/xlog,/ylog,xtitle=textoidl('\lambda [\mum]')
           break
        end
        else:  plot,/nodata,xstyle=1,ystyle=1,[1000,10],[1d-6,1d14],/xlog,/ylog,xtitle=textoidl('\lambda [\mum]'),ytitle=textoidl('B_{\nu}(T[K])\nu^{\beta}')      
     endswitch

     loadct,39,/silent
     xyouts,10^2.4,1d11,'T = '+strmid(strcompress(string(temp+0.0005),/remove_all),0,4)+'K'
     for i=beta_low,beta_high+beta_step,beta_step do begin
        res = mod_bb(temp,freq,i)
        oplot,plot_wave,res,color=230.0*i-210.d0
     endfor
     t_it++
  endfor
  
  multiplot
  plot,/nodata,xstyle=4,ystyle=4,[0,1],[0,1]
  legend,betas,textcolors=bcolors,charsize=1.5,spacing=0,/horizontal,pos=[0.1,0.5]

  device,/close

end
