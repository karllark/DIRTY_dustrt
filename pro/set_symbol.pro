; ----------------------------------------------------------------------

pro set_symbol,psym,tpsym

; psym = 0 ==> bin (solid)
;        1 ==> circle (open)
;        2 ==> circle (filled)
;        3 ==> box (open)
;        4 ==> box (filled)
;        5 ==> diamond (open)
;        6 ==> diamond (filled)
;        7 ==> triangle (open)
;        8 ==> triangle (filled)
;        9 ==> star (open)
;       10 ==> star (filled)
;       11 ==> upside down triangle (open)
;       12 ==> upside down triangle (filled)
;       13 ==> top diagonal half of square (open)
;       14 ==> top diagonal half of square (filled)
;       15 ==> bottom diagonal half of square (open)
;       16 ==> bottom diagonal half of square (filled)
;       17 ==> right facing triangle (open)
;       18 ==> right facing triangle (filled)
;       19 ==> left facing triangle (open)
;       20 ==> left facing triangle (filled)
;       21 ==> small point
;       99 ==> cross
;      100 ==> line (solid)

if (abs(psym) EQ 0) then begin
   tpsym = 10
endif
if (abs(psym) EQ 1) then begin
   a = fltarr(17)
   a(0:15) = findgen(16)*(!PI*2.0/16.0)
   a(16) = 0.0
   usersym,cos(a),sin(a)
   tpsym = 8
endif
if (abs(psym) EQ 2) then begin
   a = findgen(16)*(!PI*2.0/16.0)
   usersym,cos(a),sin(a),/fill
   tpsym = 8
endif
if (abs(psym) EQ 3) then begin
   usersym,[1.0,1.0,-1.0,-1.0,1.0],[1.0,-1.0,-1.0,1.0,1.0]
   tpsym = 8
endif
if (abs(psym) EQ 4) then begin
   usersym,[1.0,1.0,-1.0,-1.0,1.0],[1.0,-1.0,-1.0,1.0,1.0],/fill
   tpsym = 8
endif
if (abs(psym) EQ 5) then begin
   usersym,1.5*[0.0,1.0,0.0,-1.0,0.0],1.5*[1.0,0.0,-1.0,0.0,1.0]
   tpsym = 8
endif
if (abs(psym) EQ 6) then begin 
   usersym,1.5*[0.0,1.0,0.0,-1.0,0.0],1.5*[1.0,0.0,-1.0,0.0,1.0],/fill
   tpsym = 8
endif
if (abs(psym) EQ 7) then begin
   usersym,[1.0,-1.0,0.0,1.0],[-1.0,-1.0,1.0,-1.0]
   tpsym = 8
endif
if (abs(psym) EQ 8) then begin
   usersym,[1.0,-1.0,0.0,1.0],[-1.0,-1.0,1.0,-1.0],/fill
   tpsym = 8
endif
if (abs(psym) EQ 9) then begin
   usersym,[0.67,0.0,-0.67,1.0,-1.0,0.67],[-1.0,1.0,-1.0,0.33,0.33,-1.0]
   tpsym = 8
endif
if (abs(psym) EQ 10) then begin
   usersym,[0.67,0.0,-0.67,1.0,-1.0,0.67],[-1.0,1.0,-1.0,0.33,0.33,-1.0],/fill
   tpsym = 8
endif
if (abs(psym) EQ 11) then begin
   usersym,[1.0,-1.0,0.0,1.0],[1.0,1.0,-1.0,1.0]
   tpsym = 8
endif
if (abs(psym) EQ 12) then begin
   usersym,[1.0,-1.0,0.0,1.0],[1.0,1.0,-1.0,1.0],/fill
   tpsym = 8
endif
if (abs(psym) EQ 13) then begin
   usersym,[1.0,-1.0,-1.0,1.0],[1.0,-1.0,1.0,1.0]
   tpsym = 8
endif
if (abs(psym) EQ 14) then begin
   usersym,[1.0,-1.0,-1.0,1.0],[1.0,-1.0,1.0,1.0],/fill
   tpsym = 8
endif
if (abs(psym) EQ 15) then begin
   usersym,[1.0,1.0,-1.0,1.0],[1.0,-1.0,-1.0,1.0]
   tpsym = 8
endif
if (abs(psym) EQ 16) then begin
   usersym,[1.0,1.0,-1.0,1.0],[1.0,-1.0,-1.0,1.0],/fill
   tpsym = 8
endif
if (abs(psym) EQ 17) then begin
   usersym,[-1.0,-1.0,1.0,-1.0],[1.0,-1.0,0.0,1.0]
   tpsym = 8
endif
if (abs(psym) EQ 18) then begin
   usersym,[-1.0,-1.0,1.0,-1.0],[1.0,-1.0,0.0,1.0],/fill
   tpsym = 8
endif
if (abs(psym) EQ 19) then begin
   usersym,[1.0,1.0,-1.0,1.0],[1.0,-1.0,0.0,1.0]
   tpsym = 8
endif
if (abs(psym) EQ 20) then begin
   usersym,[1.0,1.0,-1.0,1.0],[1.0,-1.0,0.0,1.0],/fill
   tpsym = 8
endif
if (abs(psym) EQ 21) then begin
   usersym,[0.1,0.0,-0.1,0.0,0.1],[0.0,0.1,0.0,-0.1,0.0],/fill
   tpsym = 8
endif
if (abs(psym) GT 21) then begin
   tpsym = 1
endif
if (abs(psym) EQ 99) then begin
;   usersym,[0.5,0.5,0.5,0.0,0.1],[1.0,0.0,0.5,0.5,0.5]
   tpsym = 1
endif
if (abs(psym) EQ 100) then begin
   tpsym = 0
endif

if (psym LT 0) then begin
   tpsym = -1*tpsym
endif

end

