
void fibers_position(){
    const int Nfiber = 40;
    TPlots * plot = new TPlots();
    
    TCanvas * cVY = new TCanvas("cVY");
    plot -> H1Frame( 10 , -10,10 , "" , "x" , "y" , -10,10 );
    for (int fiber = 1 ; fiber < Nfiber+1 ; fiber++){
        TF1 * fV = new TF1(Form("fV%d",fiber),Form("2.*sqrt(2.)*(%d-21)+x",fiber),-40,40);
                fV -> SetLineColor(38);
                fV -> Draw("same");
//        TF1 * fY1 = new TF1(Form("fY1%d",fiber),Form("2.*(%d-21)+0.*x",fiber),-40,40);
//        fY1 -> SetLineColor(46);
//        fY1 -> Draw("same");
        TF1 * fY2 = new TF1(Form("fY2%d",fiber),Form("2.*(%d-21)-1.+0.*x",fiber),-11,11);
        fY2 -> SetLineColor(46);
        fY2 -> Draw("same");
//        TF1 * fU = new TF1(Form("fU%d",fiber),Form("2.*sqrt(2.)*(%d-21)-x",fiber),-11,11);
//        fU -> SetLineColor(38);
//        fU -> Draw("same");
        plot -> Line( fiber-21,-10,fiber-21,10,1,2);
        plot -> Text( 21-fiber , (2.*sqrt(2.)-1)*(fiber-21) , Form("V=%d",fiber) , 38, 0.03);
        plot -> Text( 8 , 2.*(fiber-21)-1 , Form("Y=%d",fiber) , 46, 0.03);
    }
}