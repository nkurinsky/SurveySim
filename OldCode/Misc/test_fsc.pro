PRO Example, theObject

tlb = Widget_Base(Title='Exercise FSC_FILESELECT...', Column=1)
;button = Widget_Button(tlb, Value='Make Compound Widget As Big As Me', $
;   Event_Pro='Example_Set_Size', Scr_XSize=500)
;button = Widget_Button(tlb, Value='Set Filename to worldelv.data in Data Directory', $
;   Event_Pro='Example_Set_Filename')
;button = Widget_Button(tlb, Value='Print Filename', $
;   Event_Pro='Example_Print_Filename')
;button = Widget_Button(tlb, Value="Shrink the Text Fields", $
;   Event_Pro='Example_Shrink')
CD, Current=thisDir
filenameID = FSC_FileSelect(tlb, Directory=thisDir, Filename='fileselect.pro', $
   /NoMaxSize, filter=filter,ObjectRef=theObject)
;button = Widget_Button(tlb, Value='Quit', Event_Pro='Example_Quit')
Widget_Control, tlb, /Realize, Set_UValue=theObject
XManager, 'example', tlb, /No_Block
END ;-----------------------------------------------------------------------------------------------------------------------------
