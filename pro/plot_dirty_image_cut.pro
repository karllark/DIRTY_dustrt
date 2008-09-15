;+
; NAME:
;      PLOT_DIRTY_IMAGE_CUT
;
; PURPOSE:
;      Plot a cut of the DIRTYv2 model output images
;
; CATEGORY:
;      DIRTYv2.
;
; CALLING SEQUENCE:
;      PLOT_DIRTY_IMAGE_CUT, filename
;
; INPUTS:
;      filename: name of file with image to plot (one wavelength)
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
;     Written by : Karl D. Gordon (30 Aug 2008)
;-

pro plot_dirty_image_cut,filename,y_val=y_val, $
  ps=ps,eps=eps,xrange=xrange,yrange=yrange,atten=atten,tot_cut=tot_cut

fits_open,filename,fcb

fits_read,fcb,tot_image,tot_header,exten_no=1
fits_read,fcb,stell_image,stell_header,exten_no=3
fits_read,fcb,scat_image,scat_header,exten_no=5
fits_read,fcb,input_image,input_header,exten_no=7

fits_close,fcb

input_image /= total(input_image)

image_size = size(input_image)
if (not keyword_set(y_val)) then y_val = fix(image_size[2]/2.)

pix_vals = indgen(image_size[1])

ps_file = repstr(filename,'.fits','') + '_' + strtrim(string(y_val),2)
setup_ps_output,ps_file,eps=eps,bw=bw

setup_colors,base_color,back_color,blue_color,red_color,green_color, $
             yel_color,purple_color,light_blue_color,grey_color, $
             line_color=line_color,red=red,green=green,blue=blue, $
             bw=bw

if (keyword_set(atten)) then begin
    y_vals = -2.5*alog10(tot_image[*,y_val]/input_image[*,y_val])
    indxs = where(finite(y_vals),n_indxs)
    if (n_indxs GT 0) then begin
        if (not keyword_set(xrange)) then xrange = krange(pix_vals[indxs])
        if (not keyword_set(yrange)) then yrange = krange(y_vals[indxs],kplot_type='i')
    endif
    
    ; setup the plot
    kplot,[1],[1],/no_data,xrange=xrange,yrange=yrange,kplot_type='ii', $
          xtitle=xtitle,ytitle=ytitle,color=base_color,background=back_color
    
    koplot,pix_vals,y_vals,psym=100,color=blue_color

endif else begin
    y_vals = input_image[*,y_val]
    indxs = where(finite(y_vals) AND (y_vals GT 0.0),n_indxs)
    if (n_indxs GT 0) then begin
        if (not keyword_set(xrange)) then xrange = krange(pix_vals[indxs])
        if (not keyword_set(yrange)) then yrange = krange(y_vals[indxs],kplot_type='o')
    endif
    
    ; setup the plot
    kplot,[1],[1],/no_data,xrange=xrange,yrange=yrange,kplot_type='io', $
          xtitle=xtitle,ytitle=ytitle,color=base_color,background=back_color
    
    koplot,pix_vals,scat_image[*,y_val],psym=100,color=blue_color
    koplot,pix_vals,stell_image[*,y_val],psym=100,color=red_color
    koplot,pix_vals,input_image[*,y_val],psym=100,color=green_color
    koplot,pix_vals,tot_image[*,y_val],psym=100,color=base_color

    ; temp stuff
    mod_val = exp(-abs(pix_vals-275.0)/62.625)
    mod_val = 1.0*max(input_image[*,y_val])*mod_val
    koplot,pix_vals,mod_val,color=base_color,psym=100,linestyle=2
endelse

tot_cut = input_image[*,y_val]

close_ps_output,eps=eps

end
