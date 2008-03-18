// CGS UNITS!
// Wave in cm 
// ISRF in J(lambda) erg/cm^2/s/cm/st
// HC in erg/K/gm
// Size in cm - need to compute the Heat Capacity from specific Heat capacity
// ISRF, CAbs all on wave SCALE PLEASE!


#include "TimeEvolution.h"

void TimeEvolution(vector <float>& wave, vector <float>& J, vector <float>& CAbs, 
		   vector <float>& Temperature, vector <float>& SpecificHeatCapacity,
		   vector <float>& SpecificEnthalpy, float density, float size,
		   vector <float>& Time, vector <float>& TemperatureHistory ) 
  
{
  
  // Convert specific heat capacity (erg/gm/K) to heat capacity (erg/K).
  // Convert specific enthalpy (erg/gm) to enthalpy (erg)
  vector <float> HeatCapacity,Enthalpy; 
  float GrainMass = (Constant::VFAC)*size*size*size*density;  // Grain mass in gm
  
  transform(SpecificHeatCapacity.begin(),SpecificHeatCapacity.end(),
	    back_inserter(HeatCapacity),bind1st(multiplies<float>(),GrainMass)); 
  transform(SpecificEnthalpy.begin(),SpecificEnthalpy.end(),
	    back_inserter(Enthalpy),bind1st(multiplies<float>(),GrainMass)); 

  int nWave=wave.size(); 
  int nTemp=Temperature.size(); 

  //for (int i=0;i<nTemp;i++) cout << Temperature[i] << " " << HeatCapacity[i] << endl; 
  // Compute the mean absorption time
  // (tabs)^-1 = (4pi/(hc))*int [ wave*J*Cabs ]
  float mean_photon_energy;
  vector <float> integrand(nWave);
  
  // Get tau_abs -> heating timescale
  float t_abs=HeatUtils::getTauAbs(J,CAbs,wave,mean_photon_energy); 
  
  cout << "MEAN HEATING TIME:  " << t_abs << " SEC" << endl; 
  cout << "Mean absorbed photon energy is " << mean_photon_energy << " ergs == " 
       << mean_photon_energy*(Constant::ERG_EV) << " eV." << endl; 

  // Compute the mean cooling time.   Use mean absorbed photon energy as internal energy.  
  // Cooling time = radiative (cooling) timescale.
  float Tu = HeatUtils::getTmean(Temperature,Enthalpy,mean_photon_energy); 
  cout << "TEMPERATURE CORRESPONDING TO BEING HIT BY PHOTON WITH mean_photon_energy: " << Tu << endl; 
  float t_rad = HeatUtils::getTauRad(CAbs,wave,mean_photon_energy,Tu); 
  cout << "MEAN RADIATIVE TIME SCALE:  " << t_rad << " SEC" << endl; 
  // Start time evolution. 

  // Generate TimeScale and TimeRes based on cooling and heating timescales.
  float TimeRes = t_rad/500; // 10 samples during decay of average heating.
  float SimulationTime = 100.0*t_abs; // 100 times mean absorption timescale. 
  long nTimeSamples = long(SimulationTime/TimeRes); 
  if (nTimeSamples > 5000000l) { 
    nTimeSamples = 5000000l; 
    SimulationTime = TimeRes*(float)nTimeSamples; 
  }
  //cout << "TIMERES, SIMULATIONTIME, NTIMESAMPLES: " << TimeRes << " " << SimulationTime << " " << nTimeSamples << endl; 

  // Compute number of photons absorbed in a given wavelength interval
  vector <float> nPhot(nWave); 
  vector <float> dWave(nWave);
  float Constants = TimeRes*4.0*(Constant::PI)/(Constant::PLANCKLIGHT);
  dWave[0] = (wave[1]+wave[0])/2.0 - wave[0];
  nPhot[0] = Constants*J[0]*CAbs[0]*wave[0]*dWave[0];
  for (int j=1;j<nWave-1;j++) {
    dWave[j] = (wave[j+1]-wave[j-1])/2.0;  
    nPhot[j] = Constants*J[j]*CAbs[j]*wave[j]*dWave[j];
  }
  
  dWave[nWave-1] = (wave[nWave-1]+wave[nWave-2])/2.0 - wave[nWave-2]; 
  nPhot[nWave-1] = Constants*J[nWave-1]*CAbs[nWave-1]*wave[nWave-1]*dWave[nWave-1];
  
  // Instatiate mt_ran object...
  mtran::mtrangen RanQueue( (unsigned)time ( NULL ) );
  //mtran::mtrangen RanQueue(10001l);

  float DTdtConst=TimeRes*(-4.0)*(Constant::PI); 
  vector <float> thea; 
  vector <float> thisTemp(nTimeSamples);
  // Start the grain at the temperature achieved after an "average" photon hit.
  thisTemp[0] = Tu; // Start it off hot!

  // Initialize energy accums to zero
  float finit=0.0; 
  // Total energy absorbed in a given time step
  float DU; 
  // Track whether a photon was absorbed in ANY wavelength interval 
  // during a given time step.
  bool whapped = false; 

  typedef vector <float>::iterator iter; 
  iter pnPhot=nPhot.begin(); 
  iter pwave=wave.begin();
  iter pthisTemp=thisTemp.begin(); 
  
  float thisU,thisC; 
  int uil; 

  // Hide the Runge-Kutta?
  float k1,k2,k3,k4; 
  ofstream myfile; 
  myfile.open("output.dat"); 
  myfile << thisTemp[0] << endl; 

  for (long t=1;t<nTimeSamples;t++) { // Time step loop
    //for (long t=1;t<10;t++) { // Time step loop
    
    // Generate vector of random numbers. 
    thea=RanQueue.uDev_series(nWave); 
    ////for (int i=0;i<nWave;i++) cout << i << " " << thea[i] << endl; 
    //exit(0);  
    // If thea <= nphot thea=hc/wave else thea = 0.0/wave == energy absorbed 
    // in each wavelength interval.
    // Rather than transit thea a twice with transforms (can't pass third vect), 
    // use a for loop and transit only once.   
    for (iter pthea=thea.begin();pthea!=thea.end();pthea++) {
      // *ithea = ((*ithea <= *inPhot)?PLANCKLIGHT:0.0)/(*iwave);
      // Use standard if-else so we can avoid un-necessary accums
      if (*pthea <= *pnPhot) {
	*pthea = Constant::PLANCKLIGHT/(*pwave); 
	whapped = true; 
      } else *pthea = 0.0; 
      pnPhot++; 
      pwave++;
    }

    // Here is where we sit in the tabulated U/HC
    uil = NumUtils::index(*pthisTemp,Temperature);
 
    // Total energy absorbed in this time step: 
    if (whapped) {  // We got hit.
      
      DU = accumulate(thea.begin(),thea.end(),finit); 

      //cout << uil << " " << *pthisTemp << " " << Temperature[0] << " " << Temperature[20] << endl; 
      // If temperature is above bounds, extend with T ("Dulong-Petit")
      // Below (that's a wimpy photon!) extend with T^4
      // Otherwise, interpolate table
      if (*pthisTemp >= *(Temperature.end()-1)) {
	thisU = NumUtils::lextra(*(Enthalpy.end()-1),
				 *pthisTemp,
				 *(Temperature.end()-1)); 
      } else if (*pthisTemp <= *(Temperature.begin())) {
	thisU = NumUtils::qextra(*(Enthalpy.begin()),
				 *(pthisTemp),
				 *(Temperature.begin()));
      } else {
	thisU = NumUtils::line(Temperature[uil],Temperature[uil+1], 
			       Enthalpy[uil],Enthalpy[uil+1],
			       *(pthisTemp)); 
      }

      
      cout << "OUR CURRENT TEMPERATURE IS " << *pthisTemp << " CORRESPONDING TO AN INTERNAL ENERGY/HEATC of " << thisU << " " << thisC << endl; 
      // thisU gives the Enthalpy (internal Energy) of the grains at thisTemp.
      // thisC gives the Heat Capacity of the grains at thisTemp. 

      // thisU + DU gives the new enthalpy of the grain after whapping. 
      // Find the index of (thisU+DU) in Temperature to find the new 
      // Temperature of the grain. 
      thisU += DU; 
      uil = NumUtils::index(thisU,Enthalpy); 
      if (thisU >= *(Enthalpy.end()-1)) 
	*pthisTemp = *(Temperature.end()-1)*(thisU/(*(Enthalpy.end()-1))); 
      else if (thisU <= *(Enthalpy.begin())) 
	*pthisTemp = NumUtils::frootextra(*(Temperature.begin()), 
					  thisU, 
					  *(Enthalpy.begin())); 
      else 
	*pthisTemp = NumUtils::line(Enthalpy[uil],Enthalpy[uil+1],
				    Temperature[uil],Temperature[uil+1],
				    thisU); 

      cout << "Got whapped with " << DU << " which increased our internal energy to  " << thisU << " CORRESPONDING TO A NEW TEMPERATURE OF " << *pthisTemp << endl; 
      //cout << "DTdtConst: " << DTdtConst << endl;
    }
    // Now, *pthisTemp holds the new temperature after absorption of photon(s)

    // Still need to compute thisC for non-whapped time-steps as we cool. 
    uil = NumUtils::index(*pthisTemp,Temperature); 
    //cout << "UIL: " << uil <<endl; 
    //cout << *pthisTemp << " " << *(Temperature.end()-1) << "  " << *(Temperature.begin()) << endl; 
    if (*pthisTemp >= *(Temperature.end()-1)) {
      thisC = *(HeatCapacity.end()-1); 
      //cout << "ABOVE RANGE " << thisC <<  endl; 
    } else if (*pthisTemp <= *(Temperature.begin())) {
      
      thisC = NumUtils::cextra(*(HeatCapacity.begin()), 
			       *(pthisTemp), 
			       *(Temperature.begin()));
      //cout << "BELOW RANGE " << thisC << endl;  
    } else {
      //cout << "IN RANGE" << endl; 
      //cout << "Ts, HCs, current T: " << Temperature[uil] << " " << Temperature[uil+1] << " " << 
      //	HeatCapacity[uil] << " " << HeatCapacity[uil+1]<< " " << 
      //	*(pthisTemp) << endl; 
      thisC = NumUtils::line(Temperature[uil],Temperature[uil+1], 
			     HeatCapacity[uil],HeatCapacity[uil+1],
			     *(pthisTemp));
      //cout << "IN RANGE " << thisC << endl; 
    } 

    // Four point Runge-Kutta numerical integration
    k1 = *pthisTemp;
    //cout << thisC << " " <<  "K1 is " << k1 << endl; 
    integrand = NumUtils::bbodyCGS(wave,k1);
    transform(integrand.begin(),integrand.end(),CAbs.begin(),integrand.begin(),multiplies<float>());
    //cout << NumUtils::integrate(wave,integrand) << " " << DTdtConst << " " << thisC << endl; 
    //cout << "DTCONST/THISC " << DTdtConst/thisC << endl; 
    k1 = (DTdtConst/thisC)*NumUtils::integrate(wave,integrand); 
    
    k2 = (*pthisTemp) + k1/2.0; 
    integrand = NumUtils::bbodyCGS(wave,k2);
    //for (int ii=0;ii<nWave;ii++) cout << wave[ii] << " " << integrand[ii] << endl;  
    transform(integrand.begin(),integrand.end(),CAbs.begin(),integrand.begin(),multiplies<float>());
    //for (int ii=0;ii<nWave;ii++) cout << wave[ii] << " " << integrand[ii] << endl;    
    //cout << NumUtils::integrate(wave,integrand) << " " << DTdtConst << " " << thisC << endl; 
    //cout << "DTCONST/THISC " << DTdtConst/thisC << endl; 
    k2 = (DTdtConst/thisC)*NumUtils::integrate(wave,integrand); 
    
    k3 = (*pthisTemp) + k2/2.0;
    integrand = NumUtils::bbodyCGS(wave,k3); 
    transform(integrand.begin(),integrand.end(),CAbs.begin(),integrand.begin(),multiplies<float>());
    k3 = (DTdtConst/thisC)*NumUtils::integrate(wave,integrand); 
    
    k4 = (*pthisTemp) + k3; 
    integrand = NumUtils::bbodyCGS(wave,k4); 
    transform(integrand.begin(),integrand.end(),CAbs.begin(),integrand.begin(),multiplies<float>());
    k4 = (DTdtConst/thisC)*NumUtils::integrate(wave,integrand); 
 
    //cout << "*pthisTemp " << *pthisTemp << " " << DTdtConst << " " << thisC << " " << k1 << " " <<  k2 << " " <<  k3 << " " << k4 << endl; 
    thisTemp[t] = thisTemp[t-1] + k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0; 
    //cout << t << " " << thisTemp[t] << " " << k1 << " " << k2 << " " << k3 << " " << k4 << endl; 
    //*pthisTemp++ = *pthisTemp + k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0; 

    myfile << thisTemp[t] << endl; 
    //exit(8);
    // Make sure that we don't go below equilibrium with cosmic background. 
    if (*pthisTemp < 3.0) *pthisTemp=3.0; 

    whapped = false; 
    pthisTemp++; 
    pnPhot = nPhot.begin();
    pwave = wave.begin();

  } 
  myfile.close(); 

  //for (int t=0;t<nTimeSamples;t++) cout << thisTemp[t] << endl; 

//     thea=RanQueue.uDev_series(nWave); 
//     //cout << "SIZE OF thea: " << thea.size() << endl; 
//     //for (int j=0;j<nWave;j++) cout << j << " " << thea[j] << endl; 
    
    

//   }

  //thea=RanQueue.uDev_series(10000); 
  //for (int i=0;i<2000;i++) cout << i << " " << thea[i] << endl; 
  //for (int i=0;i<10000;i++) cout << i << " " << thea[i] << endl; 
  
   // May need a way to tune these...
  //long nTimeSamp = 50000;                    // Number of time samples
  //float TimeMax = 50000.0;                     // How long (sec) to track the time evolution.
  //float TimeRes = TimeMax/float(nTimeSamp-1);  // Time Resolution 
 
  //for (long i=0;i<nTimeSamp;i++) Time.push_back(float(i)*TimeRes); 
  
  //for (long i=0;i<nTimeSamp;i++) cout << i << " " << Time[i] << endl; 
  //vector <float> freq; 
  //reverse_copy(wave.begin(),wave.end(),back_inserter(freq));
  //transform(freq.begin(),freq.end(),freq.begin(),NumUtils::inverse<float>);
  //transform(freq.begin(),freq.end(),freq.begin(),bind2nd(multiplies<float>(),(LIGHT)));
  //int nWave=wave.size(); 
  //cout << "SET UP FREQUENCY " << endl; 
  //vector <float> DeltaWave; ; 
  //DeltaWave.push_back(0.5*(wave[1]-wave[0])); 
  //for (int i=1;i<nWave-1;i++) DeltaWave.push_back( 0.5*(wave[i+1] - wave[i-1]) );
  //DeltaWave.push_back( 0.5*(wave[nWave-1]-wave[nWave-2]) );
  //cout << "SET UP DELTA WAVE" << endl; 
  //cout << "SIZE OF DELTAW " << DeltaWave.size() << endl; 
  //cout << "NWAVE " << nWave << endl; 
  //vector <float> nPhot; 
  //float Constants=TimeRes*4.0*(PI)/(PLANCKLIGHT); 
  //cout << "DEFINED NPHOT AND CONSTANTS" << endl; 
  //cout << "SIZE OF J AND CABS " << J.size() << " " << CAbs.size() << endl ;
  //for (int i=0;i<nWave;i++) nPhot.push_back(Constants*J[i]*CAbs[i]*DeltaWave[i]*wave[i]);
  //for (int i=0;i<nWave;i++) cout << i << " " << nPhot[i] << endl; 
  
  
  //cout << "Size is "  << size*1.0e4 << endl; 
  //xcout << "TIME RES IS : " << TimeRes << endl; 
}
