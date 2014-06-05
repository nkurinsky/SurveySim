;Matches BLAST bands to each other for use in later observation
;comparison. Reads in but as of yet does nothing with the model
;data. Also creates an idl save file which splits the BLAST data into
;the bands as opposed to being mixed as in the FITS image
;Current as of 6/18

pro spire_match

  cat_250 = mrdfits('/data/spire/FLS_SCAT250_v3_1_cat.fits',1,/silent)
  cat_350 = mrdfits('/data/spire/FLS_SCAT350_v3_1_cat.fits',1,/silent)
  cat_500 = mrdfits('/data/spire/FLS_SCAT500_v3_1_cat.fits',1,/silent)
  
  print,n_elements(cat_250.ra),n_elements(cat_350.ra),n_elements(cat_500.ra)

  temp = cat_500.flux

  tags = ['B250','B350','B500']
  fluxes = create_struct(tags,cat_250.flux,cat_350.flux,cat_500.flux,NAME='SPIRE Flux Densities')
  efluxes = create_struct(tags,cat_250.fluxpluserr,cat_350.fluxpluserr,cat_500.fluxpluserr,NAME='SPIRE Flux Density Errors')
  ra = create_struct(tags,cat_250.ra,cat_350.ra,cat_500.ra,NAME='Ra Coordinates (J2000)')
  dec = create_struct(tags,cat_250.dec,cat_350.dec,cat_500.dec,NAME='Dec Coordinates (J2000)')

  save,filename='save_files/SPIRE_by_band.save',description='Blast Fits File Split Into Bands',fluxes,efluxes,ra,dec

  bs = [18.1,25.2,36.6]
  bands = ['250','350','500']

  matched_fluxes = fltarr(3,n_elements(ra.B500))
  matched_fluxes_unique = fltarr(3,n_elements(ra.B500))
  matched = {B250:0,B350:0}
  matched_unique = {B250:0,B350:0}
  multiples = matched

  step = fix(n_elements(ra.B500)/20d0)

  for i=0,n_elements(ra.B500)-1 do begin
     if(i mod step eq 0) then begin
        if((i/step) mod 2 eq 0) then begin
           print,strcompress(string(fix((i/step)*5)),/remove_all),format='(a,$)'
        endif else begin
           print,'...',format='(a,$)'
        endelse
     endif
     tra = ra.B500[i]
     tdec = dec.B500[i]
     matched_fluxes(2,i) = fluxes.B500[i]
     matched_fluxes_unique(2,i) = fluxes.B500[i]
     for j=0,1 do begin
        gcirc,2,tra,tdec,ra.(j),dec.(j),dist
        gpts = where(dist lt ((bs[j]+bs[2])*0.75))
        
        inds = sort(dist[gpts])
        gpts = gpts[inds]

        if(gpts[0] ne -1) then begin 
           matched_fluxes(j,i)=fluxes.(j)[gpts[0]]
           matched.(j)++
           if(n_elements(gpts) eq 1) then begin
              matched_fluxes_unique(j,i)=fluxes.(j)[gpts[0]]
              matched_unique.(j)++
           endif else begin
              multiples.(j)++
           endelse
        endif 
     endfor
  endfor

  print,''
  print,matched
  print,multiples
  gpts = where((matched_fluxes(0,*) gt 0) and (matched_fluxes(1,*) gt 0))
  gpts2 = where((matched_fluxes_unique(0,*) gt 0) and (matched_fluxes_unique(1,*) gt 0))

  print,n_elements(gpts),n_elements(gpts2)

  range = n_elements(gpts)
  range2 = n_elements(gpts2)
  f250 = fltarr(range)
  f350 = fltarr(range)
  f500 = fltarr(range)

  for i=0,range-1 do begin
     f250[i] = matched_fluxes(0,gpts[i])
     f350[i] = matched_fluxes(1,gpts[i])
     f500[i] = matched_fluxes(2,gpts[i])
  endfor

  save,filename='save_files/matched_fluxes.save',f250,f350,f500

  f1 = fltarr(range2)
  f2 = fltarr(range2)
  f3 = fltarr(range2)

  for i=0,range2-1 do begin
     f1[i] = matched_fluxes_unique(0,gpts2[i])
     f2[i] = matched_fluxes_unique(1,gpts2[i])
     f3[i] = matched_fluxes_unique(2,gpts2[i])
  endfor

  save,filename='save_files/matched_fluxes_unique.save',f1,f2,f3

end
