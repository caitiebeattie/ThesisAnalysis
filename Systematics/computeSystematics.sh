#!/bin/bash

root -x -q EvaluateSys_low.C 
root -x -q EvaluateSys_high.C 

root -x -q 'EvaluateSys_low.C(4, -1, "loin")'     #in-plane 
root -x -q 'EvaluateSys_low.C(0, 1, "loout")'     #out-of-plane 
root -x -q 'EvaluateSys_high.C(4, -1, "hiin")'    #in-plane 
root -x -q 'EvaluateSys_high.C(0, 1, "hiout")'    #out-of-plane 

root -x -q 'plotCompareSystematics_Charged.C("lo")'
root -x -q 'plotCompareSystematics_Charged.C("hi")'
root -x -q 'plotCompareSystematics_Charged.C("loin")'
root -x -q 'plotCompareSystematics_Charged.C("loout")'
root -x -q 'plotCompareSystematics_Charged.C("hiin")'
root -x -q 'plotCompareSystematics_Charged.C("hiout")'
