#include "TCanvas.h"
#include "TROOT.h"
#include "TText.h"
#include "TVirtualPad.h"

//https://root.cern.ch/phpBB3/viewtopic.php?f=3&t=21433&start=15

UInt_t GetTextWidth(TText *t) {
   UInt_t w,h;
   Int_t f = t->GetTextFont();

   if (f%10<=2) {
      t->GetTextExtent(w,h,t->GetTitle());
   } else {
      w = 0;
      TText t2 = *t;
      t2.SetTextFont(f-1);
      TVirtualPad *pad = gROOT->GetSelectedPad();
      if (!pad) return w;
      Float_t dy = pad->AbsPixeltoY(0) - pad->AbsPixeltoY((Int_t)(t->GetTextSize()));
      Float_t tsize = dy/(pad->GetY2() - pad->GetY1());
      t2.SetTextSize(tsize);
      t2.GetTextExtent(w,h,t2.GetTitle());
   }
   return w;
}

void TextExtentTest(){
   TCanvas *c1 = new TCanvas("c1","c1",500,500);

   TText *tp2 = new TText(0.3,0.5,"HELLO");
   tp2->SetTextFont(42);
   tp2->SetTextSize(0.1);
   tp2->Draw();
   printf("tp2 width ---> %d\n",GetTextWidth(tp2));

   TText *tp3 = new TText(0.1,0.5,"HELLO");
   tp3->SetTextFont(43);
   tp3->SetTextSize(10);
   tp3->Draw();
   printf("tp3 width ---> %d\n",GetTextWidth(tp3));
}
