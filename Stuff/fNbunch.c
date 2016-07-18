
#include "TRandom2.h"

void fNbunch(){
    
    double  Rate        = 4.;    // MHz
    double  BeamRate    = 50.;   // MHz
    Double_t  mean      = Rate / BeamRate;
    Double_t  Average   = 0.;
    
    for ( int n = 0 ; n < 10 ; n++ ){
        double  probability = ((TMath::Exp(-mean))*( TMath::Power(mean,n)))/(TMath::Factorial(n));
        printf("n = %d \t probability = %f \n", n , probability);
    }
    
    printf("events:\n");
    Double_t fNbunch;
    TRandom2 * rand = new TRandom2();
    for (int event = 0 ; event < 10 ; event++ ){
        fNbunch = rand -> PoissonD(mean);
        printf("event %d \t fNbunch = %f \n", event , fNbunch);
        Average += fNbunch/10.;
    }
    printf("Average =  %f \n", Average);

}