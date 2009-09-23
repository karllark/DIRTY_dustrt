; output in W m^-2 Hz^-1

function blackbody_dirtyv2,wave,temp

; compute B(nu,T)
h = 6.6260693d-34  ; J s
c = 2.99792458d8   ; m/s
k = 1.380650d-23  ; J/K
;h = 6.63d-34  ; J s
;c = 2.998d8   ; m/s
;k = 1.38d-23  ; J/K
BnuT = 4.0*!PI*2.*(h*c/(wave*1e-6)^3)*(exp(h*c/(k*wave*1e-6*temp)) - 1.)^(-1.)

return,BnuT

end

