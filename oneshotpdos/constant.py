originaldata=[{1:"H",18:"He"},
{1:"Li",2:"Be",13:"B",14:"C",15:"N",16:"O",17:"F",18:"Ne"},
{1:"Na",2:"Mg",13:"Al",14:"Si",15:"P",16:"S",17:"Cl",18:"Ar"},
{1:"K",2:"Ca",3:"Sc",4:"Tl",5:"V",6:"Cr",7:"Mn",8:"Fe",
9:"Co",10:"Ni",11:"Cu",12:"Zn",13:"Ga",14:"Ge",15:"As",16:"Se",17:"Br",18:"Kr"},
{1:"Rb",2:"Sr",3:"Y",4:"Zr",5:"Nb",6:"Mo",7:"Tc",8:"Ru",
9:"Rh",10:"Pd",11:"Ag",12:"Cd",13:"In",14:"Sn",15:"Sb",16:"Te",17:"I",18:"Xe"},
{1:"Cs",2:"Ba",4:"Hf",5:"Ta",6:"W",7:"Re",8:"Os",
9:"Ir",10:"Pt",11:"Au",12:"Hg",13:"Tl",14:"Pb",15:"Bi",16:"Po",17:"At",18:"Rn"},
{1:"Fr",2:"Ra",4:"Rf",5:"Db",6:"Sg",7:"Bh",8:"Hs",
9:"Mt",10:"Ds",11:"Rg",12:"Cn",13:"Nh",14:"Fl",15:"Mc",16:"Lv",17:"Ts",18:"Og"},
{},
{4:"La",5:"Ce",6:"Pr",7:"Nd",8:"Pm",
9:"Sm",10:"Eu",11:"Gd",12:"Tb",13:"Dy",14:"Ho",15:"Er",16:"Tm",17:"Yb",18:"Lu"},
{4:"Ac",5:"Ta",6:"Pa",7:"U",8:"Np",
9:"Pu",10:"Am",11:"Cm",12:"Bk",13:"Cf",14:"Es",15:"Fm",16:"Md",17:"No",18:"Lr"}
]
d_swap = [{v: k for k, v in d.items()} for d in originaldata]
element_group=dict()
for d in d_swap:
    element_group.update(d)
colordata=['r','b']

ORBITAL=['s','p','d','f','g']