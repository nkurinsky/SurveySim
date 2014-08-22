;==================================================================
; Writen by Noah Kurinsky, version recent as of 8/21/14
;==================================================================

PRO SurveySim_event,ev
  COMMON simulation_com
  
  ; get event identifier
  widget_control,ev.id,get_uvalue=uvalue

  CASE uvalue OF
     'save' : save,parameters,filename='params.save' ;save parameters
     'go'  : begin              ;save settings, intialize FITS file, and pass to C++ fitting routine
        save,parameters,filename='params.save' ;save parameters
        
                                ;update observation FITS file
        widget_control,info.ot,get_value=bvals
        flux_min = [bvals[0].fmin,bvals[1].fmin,bvals[2].fmin]
        flux_err = [bvals[0].ferr,bvals[1].ferr,bvals[2].ferr]

        hdr = headfits(files.ofile)
        hnum = sxpar(hdr,"FHDU")
        hdr = headfits(files.ofile,exten=hnum)

        ;add filters
        sxaddpar, hdr, 'F1MIN',flux_min[0],'Flux cutoff, first column'
        sxaddpar, hdr, 'F2MIN',flux_min[1],'Flux cutoff, second column'
        sxaddpar, hdr, 'F3MIN',flux_min[2],'Flux cutoff, third column'
      
        modfits,files.ofile,0,hdr,exten_no=hnum

        ;make model fits file

        widget_control,info.t1,get_value=lparam
        pars = lparam(0)
        fixed = lparam(1)
        min = lparam(2)
        max = lparam(3)
        widget_control,info.t2,get_value=sparam
        sdat = sparam
        widget_control,info.t3,get_value=cparam
        cdat = cparam

        sxaddpar,hdr2,'DATE',systime(),'Date of creation'

        sxaddpar,hdr2,'PHI0',pars.phi0,'Luminosity Function Normalization'
        sxaddpar,hdr2,'PHI0_FIX',fixed.phi0,'Fix Phi0 (Y=1/N=0)'
        sxaddpar,hdr2,'PHI0_MIN',min.phi0,'Minimum Phi0 value'
        sxaddpar,hdr2,'PHI0_MAX',max.phi0,'Maximum Phi0 value'

        sxaddpar,hdr2,'L0',pars.lo,'Luminosity Function Knee'
        sxaddpar,hdr2,'L0_FIX',fixed.lo,'Fix L0 (Y=1/N=0)'
        sxaddpar,hdr2,'L0_MIN',min.lo,'Minimum L0 value'
        sxaddpar,hdr2,'L0_MAX',max.lo,'Maximum L0 value'

        sxaddpar,hdr2,'ALPHA',pars.alpha,'Luminosity Function upper slope'
        sxaddpar,hdr2,'ALPHA_FIX',fixed.alpha,'Fix Alpha (Y=1/N=0)'
        sxaddpar,hdr2,'ALPHA_MIN',min.alpha,'Minimum Alpha value'
        sxaddpar,hdr2,'ALPHA_MAX',max.alpha,'Maximum Alpha value'

        sxaddpar,hdr2,'BETA',pars.beta,'Luminosity Function lower slope'
        sxaddpar,hdr2,'BETA_FIX',fixed.beta,'Fix Beta (Y=1/N=0)'
        sxaddpar,hdr2,'BETA_MIN',min.beta,'Minimum Beta value'
        sxaddpar,hdr2,'BETA_MAX',max.beta,'Maximum Beta value'

        sxaddpar,hdr2,'P',pars.p,'Luminosity Function PHI evolution term'
        sxaddpar,hdr2,'P_FIX',fixed.p,'Fix P (Y=1/N=0)'
        sxaddpar,hdr2,'P_MIN',min.p,'Minimum P value'
        sxaddpar,hdr2,'P_MAX',max.p,'Maximum P value'

        sxaddpar,hdr2,'Q',pars.q,'Luminosity Function L evolution term'
        sxaddpar,hdr2,'Q_FIX',fixed.q,'Fix Q (Y=1/N=0)'
        sxaddpar,hdr2,'Q_MIN',min.q,'Minimum Q value'
        sxaddpar,hdr2,'Q_MAX',max.q,'Maximum Q value'

        sxaddpar,hdr2,'ZCUT',pars.zcut,'Luminosity Function z evolution limit'
        sxaddpar,hdr2,'ZCUT_FIX',fixed.zcut,'Fix ZCUT (Y=1/N=0)'
        sxaddpar,hdr2,'ZCUT_MIN',min.zcut,'Minimum ZCUT value'
        sxaddpar,hdr2,'ZCUT_MAX',max.zcut,'Maximum ZCUT value'

        sxaddpar,hdr2,'CEXP',cdat.a0,'Intrinsic luminosity evolution term'
        sxaddpar,hdr2,'CEXP_FIX',cdat.fixed,'Fix CEXP (Y=1/N=0)'
        sxaddpar,hdr2,'CEXP_MIN',cdat.amin,'Minimum CEXP value'
        sxaddpar,hdr2,'CEXP_MAX',cdat.amax,'Maximum CEXP value'

        sxaddpar,hdr2,'RUNS',sdat.runs,'Number of Runs'
        sxaddpar,hdr2,'ZMIN',sdat.zmin,'Minimum Redshift Value'
        sxaddpar,hdr2,'ZMAX',sdat.zmax,'Maximum Redshift Value'
        sxaddpar,hdr2,'DZ',sdat.dz,'Redshit Bin Width'
        sxaddpar,hdr2,'AREA',sdat.area,'Observed Solid Angle'

        sxaddpar,hdr2,'NCHAIN',msettings.nchain,'Chain Number'
        sxaddpar,hdr2,'TMAX',msettings.tmax,'Starting Anneal Temperature'
        sxaddpar,hdr2,'ANN_PCT',msettings.acceptpct,'Ideal Acceptance Percentage'
        sxaddpar,hdr2,'CONV_CONF',msettings.conv_conf,'Convergence CI Setting'
        sxaddpar,hdr2,'CONV_RMAX',msettings.conv_rmax,'Convergence Rmax Criterion'
        sxaddpar,hdr2,'CONV_STEP',msettings.conv_step,'Iterations between convergence checks'
        sxaddpar,hdr2,'BURN_STEP',msettings.burn_step,'Iterations between anneal calls in burn-in'
        sxaddpar,hdr2,'BURNVRUN',msettings.burn_ratio,'Ratio of normal to burn-in steps'
        sxaddpar,hdr2,'ANN_RNG',msettings.pct_range,'Range within which to maintain acceptance, from ideal'
        sxaddpar,hdr2,'PRINT',info.print,'Whether to Print Debug MSGs'

        templates = [0]
        mwrfits,templates,files.mfile,hdr2,/create

;===============================================================
;Run the actual simulation
;---------------------------------------------------------------
        args = files.ofile+' '+files.mfile+' '+files.sedfile+' '+files.oname
        spawn,'fitter '+args

        read_output,parameters.files.oname
        simulation_results,parameters.files.oname
     end
     'settings': begin
        settings
        plot_seds
     end
     'diag': simulation_diagnostics
     'replot': begin
        read_output,parameters.files.oname
        simulation_results,parameters.files.oname
     end
     'quit': widget_control,ev.top,/destroy
     'ot': widget_control,info.ot,get_value=bands
     't1': widget_control,info.t1,get_value=ldata
     't2': widget_control,info.t2,get_value=sdat
     't3': widget_control,info.t3,get_value=cdat
     'info': SurveySim_info
     ELSE:
  ENDCASE

END

