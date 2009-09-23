;+
; NAME:
;      GENERATE_BLACKBODY_INPUT
;
; PURPOSE:
;      Generate a stellar black body file for input to DIRTYv2.
;
; CATEGORY:
;      DIRTYv2.
;
; CALLING SEQUENCE:
;      GENERATE_BLACKBODY_INPUT
;
; INPUTS:
;      temperature : temperature in K
;      radius : radius in solar radii
;      out_filename : output filename
;
; KEYWORD PARAMETERS:
;      tot_lum : scale output to this total luminosity
;                (used for extragalactic models)
;
; OUTPUTS:
;      ASCII file suitable for input to DIRTYv2
;
; RESTRICTIONS:
;      This program written for research use.  No warrenties are
;      given.  Use at your own risk.
;
; PROCEDURE:
;     
; EXAMPLE:
;
; MODIFICATION HISTORY:
;     Written by  : Karl D. Gordon (15 May 2008)
;     31 Aug 2008 : added tot_lum keyword
;-
pro generate_blackbody_input,temperature,radius,out_filename, $
  silent=silent,tot_lum=tot_lum

solar_lum = 3.845d33

min_wave = 0.09
max_wave = 5000.0
wave_res = 1000.

;min_wave = 0.1
;max_wave = 1000.0
;wave_res = 100.

; get number of wavelength points
n_waves = fix(alog10(max_wave/min_wave)/ $
              alog10((1.0 + 2.0*wave_res)/(2.0*wave_res - 1.0)) + 1.0)
if (not keyword_set(silent)) then print,'n_waves = ',n_waves

; get the wavelength delta in log10 units
delta_log_wave = (alog10(max_wave) - alog10(min_wave))/(n_waves - 1)

; make the wavelength grid
waves = alog10(min_wave) + findgen(n_waves)*delta_log_wave
waves = 10^waves

;readcol,'wavegrid_codeval.dat',waves
;n_waves = n_elements(waves)

if (not keyword_set(silent)) then begin
    print,'min wave = ',waves[0]
    print,'min wave = ',waves[n_waves-1]
endif

; get the blackbody
bb = blackbody_dirtyv2(waves,temperature)

; convert radius from solar radii
radius_m = radius*6.95508d10

; convert to ergs s^-1 Hz^-1
bb_output = bb*1d7*4.d0*!PI*radius_m^2

; integrate to determine the total luminosity
if (keyword_set(tot_lum)) then begin
    freq = 2.998e14/waves
    bb_int = bb_output
    lum = int_tabulated(freq,bb_int,/double,/sort)
    bb_output *= tot_lum/(lum/solar_lum)
endif

freq = 2.998e14/waves
bb_int = bb_output
lum = int_tabulated(freq,bb_int,/double,/sort)
print,'total luminosity [solar] = ', lum/solar_lum

; output file
openw,unit1,out_filename,/get_lun
printf,unit1,'#Blackbody from generate_blackbody_input.pro'
printf,unit1,'#temperature = ', temperature
printf,unit1,'#radius (solar radii) = ',  radius
printf,unit1,'#'
printf,unit1,'#wavelength [microns], luminosity [ergs s^-1 Hz^-1]'
for i = 0,(n_waves-1) do begin
    printf,unit1,waves[i],bb_output[i]
endfor
free_lun,unit1

end
