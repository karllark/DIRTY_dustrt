; ----------------------------------------------------------------------

pro klegend,pos,ktype,klabel,charsize=charsize,twocol=twocol,$
            kplot_type=kplot_type,symsize=symsize,color=color, $
            line_color=line_color,thick=thick,box=box, $
            charthick=charthick,delta_y=delta_y

if (not keyword_set(color)) then begin
    color = !P.COLOR
endif else if (color EQ -1) then begin
    color = 0
endif

if (not keyword_set(line_color)) then begin
    line_color = replicate(color,n_elements(klabel))
endif

if (not keyword_set(thick)) then thick = replicate(2,n_elements(klabel))
if (not keyword_set(charthick)) then  $
  charthick = replicate(2,n_elements(klabel))

if (not keyword_set(charsize)) then begin
   charsize = 1.0
endif

if (not keyword_set(symsize)) then begin
   symsize = 1.5
endif

if (not keyword_set(kplot_type)) then begin
   kplot_type = 'ii'
endif

kxrange = !x.crange
kyrange = !y.crange

npts = n_elements(klabel)

min_x = kxrange(0)
Dx = kxrange(1) - kxrange(0)
min_y = kyrange(0)
Dy = kyrange(1) - kyrange(0)

per_x = 0.05

x_begin = pos(0)
delta_x = per_x + 0.02
y_begin = pos(1)
if (not keyword_set(delta_y)) then begin
    if (charsize LT 1.2) then begin
        delta_y = 0.04
    endif else begin
        delta_y = 0.06
    endelse
endif
max_width = 0.0
yp = y_begin + (1.2)*delta_y
yp2 = y_begin + (1.0)*delta_y
save_yp = yp
save_x_begin = x_begin
if ((npts mod 2) EQ 1) then delt_i = 0 else delt_i = 1
for i = 0,(npts-1) do begin
   x = x_begin
   x = min_x + Dx*x
   yp = yp - delta_y
   y = min_y + Dy*yp

   beg_x = x - (per_x*Dx)
   end_x = x + (per_x*Dx)
   if ((kplot_type EQ 'oo') OR (kplot_type EQ 'oi')) then begin
      beg_x = 10^(x - (per_x*Dx))
      end_x = 10^(x + (per_x*Dx))
      x = 10^x
   endif
   if ((kplot_type EQ 'io') OR (kplot_type EQ 'oo')) then begin
      y = 10^y
   endif
   
   if (ktype(1,i) EQ 100) then begin
       koplot,[beg_x,end_x],[y,y],linestyle=abs(ktype(0,i)), $
         color=line_color(i),thick=thick[i]
   endif else if (ktype(0,i) LE 0) then begin
       koplot,[beg_x,end_x],[y,y],linestyle=abs(ktype(0,i)), $
         color=line_color(i),thick=thick[i]
       koplot,[x],[y],psym=ktype(1,i),symsize=symsize,color=line_color(i)
   endif else begin
       koplot,[x],[y],psym=ktype(1,i),symsize=symsize,color=line_color(i)
   endelse

   x = x_begin + delta_x
   x = min_x + Dx*x
   yp2 = yp2 - delta_y
   y2 = min_y + Dy*yp2

   if ((kplot_type EQ 'oo') OR (kplot_type EQ 'oi')) then begin
      x = 10^x
   endif
   if ((kplot_type EQ 'io') OR (kplot_type EQ 'oo')) then begin
      y2 = 10^y2
   endif
   
   xyouts,x,y2,klabel(i),width=tst_width,charsize=charsize,color=color, $
     charthick=charthick[i]
   max_width = max([max_width,tst_width])

   if ((keyword_set(twocol)) AND ((i+delt_i) EQ fix(float(npts)/2.0))) then begin
      x_begin = x_begin + max_width + delta_x + 2.5*per_x
      yp = y_begin + (1.2)*delta_y
      yp2 = y_begin + (1.0)*delta_y
   endif

endfor

if (keyword_set(box)) then begin

   if ((keyword_set(twocol)) AND ((i mod 2) EQ 1)) then begin
       yp = yp - delta_y
   endif

    box_x1 = min_x + Dx*save_x_begin - 1.4*per_x*Dx
    box_x1_2 = min_x + Dx*x_begin - 1.4*per_x*Dx
    box_x2 = box_x1_2 + 2*delta_x*Dx + 1.2*per_x*Dx + 1.07*max_width*Dx
    box_y1 = min_y + Dy*save_yp
    box_y2 = min_y + Dy*(yp - delta_y)
    if ((kplot_type EQ 'oo') OR (kplot_type EQ 'oi')) then begin
        box_x1 = 10^box_x1
        box_x2 = 10^box_x2
    endif
    if ((kplot_type EQ 'io') OR (kplot_type EQ 'oo')) then begin
        box_y1 = 10^box_y1
        box_y2 = 10^box_y2
    endif

    koplot,[box_x1,box_x1,box_x2,box_x2,box_x1], $
      [box_y1,box_y2,box_y2,box_y1,box_y1],psym=100,color=color, $
      thick=2
endif

end

