function ask_for_file,dialog,initial_filename
  
  common ask_for_file_com,new_filename,dialog
  new_filename = initial_filename

  main = widget_base(title="Select File",/column,/align_center)
  body = widget_base(main,/column,/align_center)
  
  text = widget_text(main,/wrap,value=dialog)
  dialog = fsc_fileselect(body,/NoMaxSize,LabelName='New File')
  widget_control,dialog,set_value=new_filename
  
  btns = widget_base(main,/row,/align_left)
  save_btn = widget_button(btns,uvalue='select',Value="Select")
  cancel_btn = widget_button(btns,uvalue='skip',Value="Skip")

  widget_control,main,/realize
  xmanager,'ask_for_file',main,/no_block
  
  return new_filename
end
