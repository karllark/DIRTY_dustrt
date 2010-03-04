;+
; NAME:
;      PLOT_DIRTY_ATTEN
;
; PURPOSE:
;      Plot the attenuation curve using the output of DIRTYv2.
;
; CATEGORY:
;      DIRTYv2.
;
; CALLING SEQUENCE:
;      PLOT_DIRTY_ATTEN, filename
;
; INPUTS:
;      filename: name of ASCII FITS Table file output by DIRTYv2.
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
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
;     Written by : Karl D. Gordon (5 May 2009)
;-
pro plot_dirty_atten,filename,norm=norm,eps=eps, $
  xrange=xrange,yrange=yrange,table=table,log=log

table = mrdfits(filename,1,header)
table_tagnames = tag_names(table)

xtitle = '!4k!3 [!4l!3m]'
ytitle = 'A(!4k!3)'

; make attenuation curve
atten = -2.5*alog10(table.flux/table.flux_input)

if (keyword_set(norm)) then begin
    sindxs = sort(abs(table.wavelength - 0.55))
    print,'Att(0.55) = ', atten[sindxs[0]]
    atten = atten/atten[sindxs[0]]
endif

if (keyword_set(log)) then begin
    ykplot_type = 'o'
endif else begin
    ykplot_type = 'i'
endelse

if (not keyword_set(xrange)) then xrange = krange(table.wavelength,kplot_type='o')
if (not keyword_set(yrange)) then yrange = krange([atten,table.tau_norm],kplot_type=ykplot_type)
kplot_type = 'o' + ykplot_type

setup_ps_output,repstr(filename,'.fits','_atten'),eps=eps,bw=bw

setup_colors,base_color,back_color,blue_color,red_color,green_color, $
             yel_color,purple_color,light_blue_color,grey_color, $
             line_color=line_color,red=red,green=green,blue=blue, $
             bw=bw

; setup the plot
kplot,[1],[1],/no_data,xrange=xrange,yrange=yrange,kplot_type=kplot_type, $
      xtitle=xtitle,ytitle=ytitle,color=base_color,background=back_color

; plot the attenuation curve
koplot,table.wavelength,atten,psym=100,color=base_color,linestyle=0

; plot the input extinction curve
koplot,table.wavelength,table.tau_norm,psym=100,color=base_color,linestyle=2

close_ps_output,eps=eps

end
