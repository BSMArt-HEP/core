WriteBSMArtData[fname_]:=Block[{},

WriteString[fname,"\"BSMArtmodelname\":\""<>ModelName<>"\""];
WriteString[fname,",\"SPhenoOutputDir\":\""<>$sarahCurrentSPhenoDir<>"\""];
WriteString[fname,",\"SARAH Version\":\""<>SA`Version<>"\""];

];



WriteHBData[fname_]:=Block[{},

WriteString[fname,",\"Hneutral\":"<>ToString[SA`NumberNeutralHiggs]<>",\"Hcharged\":"<>ToString[SA`NumberChargedHiggs]];
tlist=Select[ParticleDefinitions[EWSB],#[[1]]==HiggsBoson&];
If[Length[tlist]!=0,hbpdgs=DeleteCases[PDG/.tlist[[1,2]],0],
hbpdgs={}];
tlist=Select[ParticleDefinitions[EWSB],#[[1]]==PseudoScalar&];
If[Length[tlist]!=0,pspdgs=DeleteCases[PDG/.tlist[[1,2]],0],
pspdgs={}];
neutralpdgs=Join[hbpdgs,pspdgs];
tlist=Select[ParticleDefinitions[EWSB],#[[1]]==ChargedHiggs&];
If[Length[tlist]!=0,chargedpdgs=Abs/@DeleteCases[PDG/.tlist[[1,2]],0],
chargedpdgs={}];
WriteString[fname,",\"Neutral Higgs\":["<>StringJoin[Riffle[ToString /@ neutralpdgs, ","]]<>"]\n"];
WriteString[fname,",\"Charged Higgs\":["<>StringJoin[Riffle[ToString /@ chargedpdgs, ","]]<>"]\n"];

];

MakeQNUMBERSBSMArt[qfile_, oddprtcls_ : {}] := Module[{OF, plist, i, j, part, tpart, goldstones, pdg, SMpdgs, 
   allparts, Q3, spin, colour, isreal, getanti, suffix, suffix1, pname1,pname2,
   suffix2, isQCDconjugated},
  
  (* PART[S], PART[F] and PART[
  V] store the necessary mass eigenstates (plus PART[A] and PART[
  G] which we don't need). 
  They are recomputed every time the vertices are computed, 
  so we are necessarily working with the latest ones. *)
  
  (*
  Block QNUMBERS 9900213 # rho+1 3 # 3 times electric charge
  1 3 # 3 times electric charge
  2 3 # number of spin states (2S+1)
  3 1 # colour rep (1:singlet,3:triplet,8:octet)
  4 1 # Particle/Antiparticle distinction (0=own anti)
    11  0 # Z2 parity: 0 even, 1 odd (for SModelS) 
  *)
  SMpdgs = {0, 1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25};
  OF = OpenWrite[qfile];
  plist = Join[Transpose[PART[S]][[1]], Transpose[PART[F]][[1]], Transpose[PART[V]][[1]]];
  goldstones = Transpose[GoldstoneGhost][[2]];
  For[i = 1, i <= Length[plist], i++,
   part = plist[[i]];
   (*Print[part];*)
   For[j = 1, j <= getGen[part], j++,
    If[!SubsetQ[goldstones, {part}] && FreeQ[goldstones, part[{j}]],
      pdg = getPDG[part, j];
      If[FreeQ[SMpdgs, Abs[pdg]],
       getanti = Sign[pdg];
       If[getanti < 0,
        tpart = AntiField[part],
        tpart = part;
        ];
       (*Print["BSM: ",pdg];*)
       isQCDconjugated = 1;
       isreal = "1";
       If[AntiField[part] === part,
        isreal = "0";
        ];
       colour = getColorDim[part];
       If[colour > 1 && ConjugatedRepQ[SA`DynL[tpart, color], SU[3]],
        colour *= -1;
        ];

        If[getanti < 0,
        pname1=getOutputNameAnti[part,j];
        pname2=getOutputName[part,j];
        ,
        pname1=getOutputName[part,j];
        pname2=getOutputNameAnti[part,j];
        ];


        (*

       If[getGen[part] > 1,
        suffix = ToString[j];
        ,
        suffix = "";
        ];
       
       If[getanti < 0,
        suffix1 = suffix <> "bar";
        suffix2 = suffix <> "";
        ,
        suffix1 = suffix <> "";
        suffix2 = suffix <> "bar"];
       
      *)

       (*WriteString[OF, "BLOCK QNUMBERS " <> ToString[Abs[pdg]] <> " # " <> ToString[part] <>suffix1];
       If[isreal == "1", WriteString[OF, " " <> ToString[part] <> suffix2];];
       *)
       WriteString[OF, "BLOCK QNUMBERS " <> ToString[Abs[pdg]] <> " # " <> ToString[pname1]];
        If[isreal == "1", WriteString[OF, " " <> ToString[pname2]];];
       WriteString[OF, "\n"];
       Q3 = 3*getElectricCharge[tpart]; 
       WriteString[OF, "   1  " <> ToString[Q3] <> " # 3 times electric charge\n"];
       Switch[getType[part],
        S,
        spin = 1,
        F,
        spin = 2,
        V,
        spin = 3,
        _,
        Print["Something went wrong with the type of particle ", part];
        Return[];
        ];
       WriteString[OF, 
        "   2  " <> ToString[spin] <> 
         " # number of spin states (2S+1)\n"];
       
       WriteString[OF, 
        "   3  " <> ToString[colour] <> 
         " # colour rep (1: singlet, 3: triplet, 8: octet)\n"];
       
       WriteString[OF, 
        "   4  " <> isreal <> 
         " # Particle/Antiparticle distinction (0=own anti)\n"];  
       If[(! FreeQ[oddprtcls, part]) || getRParity[part, SA`CurrentStates] == -1,
        WriteString[OF, 
          "  11  1 # Z2 parity: 0 even, 1 odd (for SModelS) \n"];
        ,
        WriteString[OF, 
          "  11  0 # Z2 parity: 0 even, 1 odd (for SModelS) \n"];
        ];
       ]
      ];
    ];
   
   
   ];
  Close[qfile];
  ];
  
