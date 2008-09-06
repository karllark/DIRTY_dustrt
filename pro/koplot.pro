; -----------------------------------------------------------------------

pro koplot,x,y,psym=psym,linestyle=linestyle,symsize=symsize,color=color, $
           xerror=xerror,yerror=yerror,thick=thick,errstyle=errstyle

if (not keyword_set(thick)) then thick = 2
if (not keyword_set(color)) then begin
    color = !P.COLOR
endif else if (color EQ -1) then begin
    color = 0
endif

if (not keyword_set(psym)) then begin
   psym = 0
endif

if (not keyword_set(symsize)) then begin
   symsize = 1.5
endif

if (not keyword_set(linestyle)) then begin
   linestyle = 0
endif

linestyle = linestyle - 5*fix(float(linestyle)/5.0)

set_symbol,psym,tpsym

oplot,x,y,psym=tpsym,linestyle=linestyle,symsize=symsize,color=color, $
  thick=thick
if (keyword_set(yerror)) then begin
    if (keyword_set(xerror)) then begin
        koploterr,x,y,xerror,yerror,psym=3,errcolor=color, $
          linestyle=linestyle,thick=thick,errstyle=errstyle
    endif else begin
        koploterr,x,y,yerror,psym=3,errcolor=color,thick=thick,errstyle=errstyle
    endelse
endif
    

end

