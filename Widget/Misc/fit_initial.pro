PRO fit_initial_event,ev
  COMMON alldata,data,t1a,t1b,t2,t3,file,file2,oname

  widget_control,ev.id,get_uvalue=uvalue

  CASE uvalue OF
     'go'  : begin
        print,'GO button'
        widget_control,t1a,get_value=value
        print,value
        widget_control,t1b,get_value=value
        print,value
        widget_control,t2,get_value=value
        print,value
        widget_control,t3,get_value=value
        print,value
        widget_control,file,get_value=value
        print,value
        widget_control,file2,get_value=value
        print,value
        wigdet_control,oname,get_value=value
        print,value
     end
     'draw': print,'draw event',ev.x,ev.y,ev.press,ev.release
     'quit': begin
        widget_control,ev.top,/destroy
     end
     'tab':
     'tab1':
     'tab2':
     'info': fit_info
     ELSE: print,'unidentified'
  ENDCASE
END

PRO fit_initial
  COMMON alldata,data,t1a,t1b,t2,t3,file,file2,oname

  main = widget_base(title='Initial Fitting Parameters',/row)

  p_main = widget_base(main,/column)
  g_main = widget_base(main,/column)

  r1 = widget_base(p_main,/row)
  r3 = widget_base(p_main,/column,/frame)
  r2 = widget_base(p_main,/row)

  btn = widget_button(r2,uvalue='go',value='Go',xsize=50,ysize=25)
  btn2 = widget_button(r2,uvalue='quit',value='Quit',xsize=50,ysize=25)
  btn3 = widget_button(r2,uvalue='info',value='Info',xsize=50,ysize=25)
  file = fsc_fileselect(r3,LabelName='Observation Save File: ',/mustexist,filter='*.save',filename='observation.save')
  file2 = fsc_fileselect(r3,LabelName='Model Save file: ',/mustexist,filter='*.save',filename='model.save')
  oname = fsc_fileselect(r3,LabelName='Output FITS file: ',filename='output.fits')

  zdat = {min:0.0,max:1.0,mean:0.5,sigma:0.1}
  pdat = {min:0,max:10,mean:3,sigma:2}
  data = [pdat,pdat,pdat]
  data1 = [[0.728],[0.272],[73]]
  rows = ["p1","p2","p3"]
  rows1 = ["Lambda0","OmegaM","H0"]
  rows2 = ["Phi0","L0","alpha","beta","p","q"]
  data2 = [[-2.2],[23.64],[3],[4],[3],[-6]]

  fitting_table = widget_base(r1,/column)
  l1 = widget_label(fitting_table,value="Fitting Parameters")
  t1a = widget_table(fitting_table,value=zdat,column_labels=tag_names(zdat),row_labels=["z"],uvalue='t1a',/editable,alignment=1)
  t1b = widget_table(fitting_table,value=data,/no_column_headers,row_labels=rows,uvalue='t1b',/editable,alignment=1)

  cosmo_table = widget_base(r1,/column)
  l2 = widget_label(cosmo_table,value="Cosmological Parameters")
  t2 = widget_table(cosmo_table,value=data1,column_labels=['Initial Value'],row_labels=rows1,uvalue='tab2',/editable,alignment=1,column_widths = 100)

  lum_table = widget_base(r1,/column)
  l3 = widget_label(lum_table,value="Luminosity Function Parameters")
  t3 = widget_table(lum_table,value=data2,column_labels=['Initial Value'],row_labels=rows2,uvalue='tab3',/editable,alignment=1,column_widths = 100)

  widget_control,main,/realize
  xmanager,'fit_initial',main,/no_block

END

pro fit_info_event,ev

  widget_control,ev.id,get_uvalue=uvalue

  case uvalue of
     'close': begin
        widget_control,ev.top,/destroy
     end
  endcase

end

pro fit_info
  
  main = widget_base(title="Information about Fitting Methodology",/column)
  t1 = widget_text(main,value="this is text")
  bt = widget_button(main,uvalue='close',value='Close',xsize=50,ysize=25)

  widget_control,main,/realize
  xmanager,'fit_info',main,/no_block

end
