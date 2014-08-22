function make_info_struct,parnum
  
  info={SurveySimInformation, $
        base: 0L, $
        base1: 0L, $
        base2: 0L, $
        base_out: 0L, $
        base1_out: 0L, $
        base2_out: 0L, $
        xsize1:150., $
        ysize1:100., $
        xsize2:150., $
        ysize2:100., $
        magnification:1.00, $  
        draw1:0L, $
        draw2:0L, $
        draw3:0L, $
        win_id1:0L, $
        win_id2:0L, $
        win_id3:0L, $
        ncolors:0L, $
;widget bases
        p_main:0L, $
        obs_table:0L, $
        lum_table:0L, $
        sed_table:0L, $
        sim_table:0L, $
        dbase:0L, $
        button_base:0L, $
;input tables and files
        obsname:0L, $
        sfile:0L, $
        oname:0L, $
        mname:0L, $
        ot:0L, $
        fd1:0L, $
        fd2:0L, $
        fd3:0L, $
        t1:0L, $
        fixinfo: make_array(parnum,value=0L), $
        t2:0L, $
        t3:0L, $
        tset:0L, $
        dprint:0L, $
        print:1}
  
  return,info
end
