pro SurveySim_info
  
  spawn,'cat SurveySim_info.txt',info

  main = widget_base(title="Information about SurveySim IDL Widget",/column)
  t1 = widget_text(main,/wrap,value=info,xsize=80,ysize=20)
  bthold = widget_base(main,/row,/align_center)
  bt = widget_button(bthold,uvalue='close',value='Close',xsize=50,ysize=25)

  widget_control,main,/realize
  xmanager,'SurveySim_info',main
  
end
