; ----------------------------------------------------------------------

pro kplot,x,y,title=title,xtitle=xtitle,ytitle=ytitle,psym=psym,$
          xrange=xrange,yrange=yrange,position=position, $
          linestyle_in=linestyle_in,nodata=nodata,$
          charsize=charsize,kplot_type=kplot_type,error_bars=error_bars,$
          symsize=symsize,xstyle=xstyle,ystyle=ystyle,noerase=noerase, $
          color=color,xtickname=xtickname,ytickname=ytickname,  $
          xticks=xticks,yticks=yticks,background=background, $
          xerror=xerror,yerror=yerror,thick=thick, $
          no_position=no_position,_extra=extra,charthick=charthick,errstyle=errstyle

if (not keyword_set(thick)) then begin
    thick = 2
    xthick = thick
    ythick = thick
    charthick = thick
endif else begin
;    thick = !P.thick
;    xthick = !P.xthick
;    ythick = !P.ythick
;    charthick = !P.charthick
    xthick = thick
    ythick = thick
    if (not keyword_set(charthick)) then charthick = thick
endelse

if (not keyword_set(xticks)) then xticks = 0
if (not keyword_set(yticks)) then yticks = 0
if (not keyword_set(xtickname)) then xtickname = ''
if (not keyword_set(ytickname)) then ytickname = ''

if (not keyword_set(color)) then begin
    color = !P.COLOR
endif else if (color EQ -1) then begin
    color = 0
endif

if (not keyword_set(background)) then background = !P.BACKGROUND
if (not keyword_set(kplot_type)) then kplot_type = 'ii'
if (not keyword_set(title)) then title = ''
if (not keyword_set(xtitle)) then xtitle = ''
if (not keyword_set(ytitle)) then ytitle = ''
if (not keyword_set(psym)) then psym = 0
if (not keyword_set(xstyle)) then xstyle = 1
if (not keyword_set(ystyle)) then ystyle = 1

if (not keyword_set(xrange)) then begin
   if ((kplot_type EQ 'oo') OR (kplot_type EQ 'oi')) then begin
      xrange = krange(x,kplot_type='o')
   endif else begin
      xrange = krange(x)
   endelse
endif

if (not keyword_set(yrange)) then begin
   if ((kplot_type EQ 'io') OR (kplot_type EQ 'oo')) then begin
      yrange = krange(y,kplot_type='o')
   endif else begin
      yrange = krange(y)
   endelse
endif

if (not keyword_set(position)) then begin
   position = [0.15,0.15,0.95,0.95]
endif
if (keyword_set(no_position)) then begin
    position = 0
endif

if (not keyword_set(charsize)) then charsize = 1.4
if (not keyword_set(symsize)) then symsize = 1.5
if (not keyword_set(linestyle_in)) then linestyle_in = 0
linestyle = linestyle_in
if (not keyword_set(error_bars)) then error_bars = 'no'

set_symbol,psym,tpsym

linestyle = linestyle - 5*fix(float(linestyle)/5.0)

if (kplot_type EQ 'io') then begin
    plot_io,x,y,psym=tpsym,title=title,xtitle=xtitle,ytitle=ytitle,$
      position=position,xrange=xrange,yrange=yrange,nodata=nodata,$
      xstyle=xstyle,ystyle=ystyle,linestyle=linestyle,charsize=charsize,$
      symsize=symsize,noerase=noerase,color=color,background=background, $
      xtickname=xtickname,ytickname=ytickname,charthick=charthick, $
      xticks=xticks,yticks=yticks,thick=thick,xthick=xthick,ythick=ythick, $
      _extra=extra
endif else if (kplot_type EQ 'oi') then begin
    plot_oi,x,y,psym=tpsym,title=title,xtitle=xtitle,ytitle=ytitle,$
      position=position,xrange=xrange,yrange=yrange,nodata=nodata,$
      xstyle=xstyle,ystyle=ystyle,linestyle=linestyle,charsize=charsize,$
      symsize=symsize,noerase=noerase,color=color,background=background, $
      xtickname=xtickname,ytickname=ytickname,charthick=charthick, $
      xticks=xticks,yticks=yticks,thick=thick,xthick=xthick,ythick=ythick, $
      _extra=extra
endif else if (kplot_type EQ 'oo') then begin
    plot_oo,x,y,psym=tpsym,title=title,xtitle=xtitle,ytitle=ytitle,$
      position=position,xrange=xrange,yrange=yrange,nodata=nodata,$
      xstyle=xstyle,ystyle=ystyle,linestyle=linestyle,charsize=charsize,$
      symsize=symsize,noerase=noerase,color=color,background=background, $
      xtickname=xtickname,ytickname=ytickname,charthick=charthick, $
      xticks=xticks,yticks=yticks,thick=thick,xthick=xthick,ythick=ythick, $
      _extra=extra
endif else begin
    plot,x,y,psym=tpsym,title=title,xtitle=xtitle,ytitle=ytitle,$
      position=position,xrange=xrange,yrange=yrange,nodata=nodata,$
      xstyle=xstyle,ystyle=ystyle,linestyle=linestyle,charsize=charsize,$
      symsize=symsize,noerase=noerase,color=color,background=background, $
      xtickname=xtickname,ytickname=ytickname,charthick=charthick, $
      xticks=xticks,yticks=yticks,thick=thick,xthick=xthick,ythick=ythick, $
      _extra=extra
endelse
if (keyword_set(yerror)) then begin
    if (keyword_set(xerror)) then begin
        koploterr,x,y,xerror,yerror,psym=3 ;,errcolor=color,thick=thick,errstyle=errstyle
    endif else begin
        koploterr,x,y,yerror,psym=3,errcolor=color,thick=thick,errstyle=errstyle
    endelse
endif

end

