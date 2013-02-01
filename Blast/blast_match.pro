;Matches BLAST bands to each other for use in later observation
;comparison. Reads in but as of yet does nothing with the model
;data. Also creates an idl save file which splits the BLAST data into
;the bands as opposed to being mixed as in the FITS image
;Current as of 6/22

pro blast_match

  ;data = mrdfits('model.fits',1,/silent)
  data2 = mrdfits('/data/blast_goods/devlin09_catalog.fit',1,/silent)
  
  ;wave = data.wavelength
  ;models = data.models

  name = data2.blast
  band = data2.band
  f_obs = data2.flux
  ef_obs = data2.noise
  ra = data2.raJ2000
  dec = data2.deJ2000

  b1 = where(band eq 250)
  b2 = where(band eq 350)
  b3 = where(band eq 500)

  print,n_elements(b1),n_elements(b2),n_elements(b3)
  
  name = data2.blast
  band = data2.band
  f_obs = data2.flux
  ef_obs = data2.noise
  b_ra = data2.raJ2000
  b_dec = data2.deJ2000

  tags = ['B250','B350','B500']
  names = create_struct(tags,name[b1],name[b2],name[b3],NAME='BLAST IDs')
  fluxes = create_struct(tags,f_obs[b1],f_obs[b2],f_obs[b3],NAME='Flux Densities')
  efluxes = create_struct(tags,ef_obs[b1],ef_obs[b2],ef_obs[b3],NAME='Flux Density Errors')
  ra = create_struct(tags,b_ra[b1],b_ra[b2],b_ra[b3],NAME='Ra Coordinates (J2000)')
  dec = create_struct(tags,b_dec[b1],b_dec[b2],b_dec[b3],NAME='Dec Coordinates (J2000)')

  save,filename='save_files/BLAST_by_band.save',description='Blast Fits File Split Into Bands',names,fluxes,efluxes,ra,dec

  bs = [30,41,59]
  bands = ['250','350','500']

  sf = f_obs[b1]
  esf = ef_obs[b1]
  fmin = min(sf,imin,max=fmax,subscript_max=imax)
  
  print,fmin,sf[imin],esf[imin]
  print,fmax,sf[imax],esf[imax]

  matched_fluxes = fltarr(3,n_elements(names.B250))
  matched = {B350:0,B500:0}
  multiples = matched

  for i=0,n_elements(names.B250)-1 do begin
     tname = names.B250[i]
     tra = ra.B250[i]
     tdec = dec.B250[i]
     matched_fluxes(0,i) = fluxes.B250[i]
     for j=1,2 do begin
        tband = bands[j]
        gcirc,2,tra,tdec,ra.(j),dec.(j),dist
        gpts = where(dist lt ((bs[j]+bs[0])*0.75))
        if(gpts[0] ne -1) then begin 
           matched_fluxes(j,i)=fluxes.(j)[gpts[0]]
           matched.(j-1)++
        endif else begin
           ;print,"No Match: ",tname,"   ",tband
        endelse
        if(n_elements(gpts) gt 1) then multiples.(j-1)++
     endfor
  endfor

  print,matched
  print,multiples
  gpts = where((matched_fluxes(1,*) gt 0) and (matched_fluxes(2,*) gt 0))

  print,n_elements(gpts),gpts[n_elements(gpts)-1]

  range = n_elements(gpts)
  f250 = fltarr(range)
  f350 = fltarr(range)
  f500 = fltarr(range)

  for i=0,range-1 do begin
     f250[i] = matched_fluxes(0,gpts[i])
     f350[i] = matched_fluxes(1,gpts[i])
     f500[i] = matched_fluxes(2,gpts[i])
  endfor

  save,filename='save_files/matched_fluxes.save',f250,f350,f500

  for i=0,9 do begin
     print,f250[i],f350[i],f500[i]
  endfor

  write_fits_obs,f250,f350,f500,[250d-6,350d-6,500d-6],observatory='BLAST'

end
