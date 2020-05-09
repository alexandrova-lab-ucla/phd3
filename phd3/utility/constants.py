#!/usr/bin/env python3
"""
Author  ==>> Matthew R. Hennefarth
Date    ==>> April 16, 2020
"""

__all__=[
    'AMINO_ACID_RESIDUES',
    'METALS',
    'THERMOSTATS',
    'PROTONATED',
    'DEPROTONATED',
    'HEAVY_ATOMS',
    'A_TO_BOHR',
    'AVAILABLE_BASIS',
    'AVAILABLE_GRIDS',
    'AVAILABLE_FUNCS',
    'MINN_FUNCS',
    'MFILE',
    'SUBMIT_FILE_NAME',
    'AVAILABLE_CALCULATIONS',
    'QUOTES',
    'ATOM_MASS',
    'Kb',
    'H_TO_KCAL',
    'MO_FILES',
    'PROTONATED_STANDARD',
    'DEPROTONATED_STANDARD'
]

MO_FILES = ['mos', 'alpha', 'beta']

AMINO_ACID_RESIDUES = ('ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL')

THERMOSTATS = ('ANDERSON')

# Given a Residue, it can get the heavy atom and the proton that can be either protonated or deprotonated
PROTONATED = {
    "ASP": [("OD1", "2HND"), ("OD2", "1HND")],
    "GLU": [("OE1", "2HNE"), ("OE2", "1HNE")],
    "HIS": [("NE2", "2HNE")]
}

PROTONATED_STANDARD = {
    "ASP" : [("OD1", "HD1"), ("OD2", "HD2")],
    "GLU" : [("OE1", "HE2"), ("OE2", "HE1")],
    "HIS" : [("NE2", "HE2")]
}

#Need to add arganine and lysine
DEPROTONATED = {
    "SER": [("OG", "HO")],
    "CYS": [("SG", "HG1")],
    "THR": [("OG1", "HO")],
    "ASN": [("ND2", "1HND"), ("ND2", "2HND")],
    "GLN": [("NE2", "1HNE"), ("NE2", "2HNE")],
    "TYR": [("OH", "HO")],
    "TRP": [("NE1", "HE1")],
    "HIS": [("ND1", "HD1")],
    "ARG": [("NH1", "2HH1"), ("NH1", "1HH1"), ("NH2", "2HH2"), ("NH2", "1HH2"), ("NE", "HE")],
    "LYS": [("NZ", "HZ1"), ("NZ", "HZ2"), ("NZ", "HZ3")]
    }

DEPROTONATED_STANDARD = {
    "SER" : [("OG", "HG")],
    "CYS" : [("SG", "HG")],
    "THR" : [("OG1", "HG1")],
    "ASN" : [("ND2", "HD21"), ("ND2", "HD22")],
    "GLN" : [("NE2", "HE21"), ("NE2", "HE22")],
    "TYR" : [("OH", "HH")],
    "TRP" : [("NE1", "HE1")],
    "HIS" : [("ND1", "ND1")],
    "ARG" : [("NH1", "HH21"), ("NH1", "HH11"), ("NH2", "HH22"), ("NH2", "HH12"), ("NE", "HE")],
    "LYS" : [("NZ", "HZ1"), ("NZ", "HZ2"), ("NZ", "HZ3")]
}

HEAVY_ATOMS = ['n', 'o', 's', 'se']

A_TO_BOHR = 1.8897259886

#All metals through bismuth
METALS = ('li','be','na','mg','al','k','ca','sc','ti','v','cr','mn','fe','co','ni','cu','zn','ga','rb','sr','y','zr','nb','mo','tc','ru','rh','pd','ag','cd','in','sn','cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy','ho','er','tm','yb','lu','hf','ta','w','re','os','ir','pt','au','hg','tl','pb','bi')

H_TO_KCAL = 627.509

ATOM_MASS = {
    'cl' : 35.453,
    'na' : 22.989769,
    'k' : 39.0983,
    'rb' : 85.4678,
    'cs' : 132.90545
}

Kb = 1.3806e-23

AVAILABLE_FUNCS = [
    'tpss', 'tpssh', 's-vwn', 'b97-d', 'pbe0',
    'slater-dirac-exchange', 'b2-plyp', 'vwn',
    's-vwn_Gaussian', 'pwlda', 'becke-exchange',
    'b-lyp', 'b-vwn', 'lyp', 'b-p', 'pbe', 'bh-lyp',
    'b3-lyp', 'b3-lyp_Gaussian', 'lhf', 'm06-2x', 'm06'
]

AVAILABLE_BASIS = [
    'def2-SVP', 'def2-TZVPP', 'def-SVP', 'def-SV(P)', 'def-TZVP',
    'def2-SV(P)', 'def2-TZVP', 'def-TZVPP', 'def2-QZVP', 'SVP',
    'SV', 'SV(P)', 'TZVP'
]

AVAILABLE_GRIDS = [
    'm3', 'm4', 'm5', '1', '2', '3', '4', '5'
]

MINN_FUNCS = ['m06-2x', 'm06']

MFILE = "MFILE"

SUBMIT_FILE_NAME = "submit.sh"

AVAILABLE_CALCULATIONS = ["woelfling", "geo", "numforce", "trans", "escf", "forceopt", "singlepoint", "sp", "nf"]

QUOTES = [
    R'''
                             .-"""-.    __                         
                            /       \.-"  "-.                      
                         __:  :\     ;       `.                    
                  _._.-""  :  ; `.   :   _     \                   
                .'   "-.  "   :   \ /;    \    .^.              .-,
    .-".       :        `.     \_.' \'   .'; .'   `.            `dP
 ,-"    \      ;\         \  '.     /". /  :/       `.      __ dP_,
 :    '. \_    ; `.  __.   ;_  `-._/   Y    \         `.   ( dP".';
 ;      \  `.  :   "-._    ; ""-./      ;    "-._       `--dP .'  ;
:    .--.;   \  ;      l   '.    `.     ;        ""--.   dP  /   / 
;   /    :    \/       ;\  . "-.   \___/            __\dP .-"_.-"  
;  /     L_    \`.    :  "-.J   "-._/  ""-._       ( dP\ /   /     
: :      ; \    `.`.  ;     /"+.     ""-.   ""--.._dP-, `._."      
 \;     :   \     `.`-'   _/ /  "-.___   "         \`-'            
  `.    ;    \      `._.-"  (     ("--..__..---g,   \              
    `. :      ;             /\  .-"\       ,-dP ;    ;             
      \;   .-';    _   _.--"  \/    `._,-.-dP-' |    ;             
       :     :---"" """        `.     _:'.`.\   :    ;\            
        \  , :              bug  "-. (,j\ ` /   ;\(// \\           
         `:   \                     "dP__.-"    '-\\   \;          
           \   :                .--dP,             \;              
            `--'                `dP`-'                             
                              .-j                                  
                              `-:_                                 
                                 \)                                
                                  `--'
''',
    "You gotta grep your greps and sed your seds - mike nechay",
    "If you installed Anaconda, i will never speak to you again - Kirill Shumilov",
    "Yo, this [stuff] shmacks - Zerina Mehmedovic",
    "As David alwas says, 'If you ain't cheating you aint trying' - Jack Terrell Fuller III",
    "Happy Science Day - Jack Terell Fuller III",
    "Its coke zero time - Hang Yin",
    "Your face has two barriers",
    "It is Wednesday my dudes",
    "Does it look like I know what a jpeg is? - Hank Hill",
    "Use the force general Kobi, you're my only hope (The Emperor's Strikes Back)",
    "Wack Jack in da Shack",
    "Guess you just gotta axsk'em",
    "Boi",
    R'''
                     ,----..            ,--.               ,---.'|                                         
    ,---,       /   /   \         ,--.'|   ,---,       |   | :       ,---,                             
  .'  .' `\    /   .     :    ,--,:  : |  '  .' \      :   : |     .'  .' `\                           
,---.'     \  .   /   ;.  \,`--.'`|  ' : /  ;    '.    |   ' :   ,---.'     \                          
|   |  .`\  |.   ;   /  ` ;|   :  :  | |:  :       \   ;   ; '   |   |  .`\  |                         
:   : |  '  |;   |  ; \ ; |:   |   \ | ::  |   /\   \  '   | |__ :   : |  '  |                         
|   ' '  ;  :|   :  | ; | '|   : '  '; ||  :  ' ;.   : |   | :.'||   ' '  ;  :                         
'   | ;  .  |.   |  ' ' ' :'   ' ;.    ;|  |  ;/  \   \'   :    ;'   | ;  .  |                         
|   | :  |  ''   ;  \; /  ||   | | \   |'  :  | \  \ ,'|   |  ./ |   | :  |  '                         
'   : | /  ;  \   \  ',  / '   : |  ; .'|  |  '  '--'  ;   : ;   '   : | /  ;                          
|   | '` ,/    ;   :    /  |   | '`--'  |  :  :        |   ,/    |   | '` ,/                           
;   :  .'       \   \ .'   '   : |      |  | ,'        '---'     ;   :  .'                             
|   ,.' ,----,   `---`     ;   |.'      `--''                    |   ,.'                               
'---' ,/   .`|             '---'                 ____  ,-.----.  '---'                                 
    ,`   .'  :,-.----.                         ,'  , `.\    /  \                                       
  ;    ;     /\    /  \           ,--,      ,-+-,.' _ ||   :    \                                      
.'___,/    ,' ;   :    \        ,'_ /|   ,-+-. ;   , |||   |  .\ :                                     
|    :     |  |   | .\ :   .--. |  | :  ,--.'|'   |  ;|.   :  |: |                                     
;    |.';  ;  .   : |: | ,'_ /| :  . | |   |  ,', |  ':|   |   \ :                                     
`----'  |  |  |   |  \ : |  ' | |  . . |   | /  | |  |||   : .   /                                     
    '   :  ;  |   : .  / |  | ' |  | | '   | :  | :  |,;   | |`-'                                      
    |   |  '  ;   | |  \ :  | | :  ' ; ;   . |  ; |--' |   | ;                                         
    '   :  |  |   | ;\  \|  ; ' |  | ' |   : |  | ,    :   ' |                                         
    ;   |.'   :   ' | \.':  | : ;  ; | |   : '  |/     :   : :                                         
    '---'     :   : :-'  '  :  `--'   \;   | |`-'      |   | :                                         
              |   |.'    :  ,      .-./|   ;/          `---'.|                                         
              `---'       `--`----'    '---'             `---`

    ''',
    R'''
        ___           ___           ___           ___           ___              
     /\  \         /\  \         /\  \         /\  \         /\__\             
    /::\  \       /::\  \       /::\  \       /::\  \       /:/  /             
   /:/\:\  \     /:/\:\  \     /:/\:\  \     /:/\:\  \     /:/__/              
  /:/  \:\  \   /::\~\:\  \   /::\~\:\  \   /:/  \:\  \   /::\__\____          
 /:/__/ \:\__\ /:/\:\ \:\__\ /:/\:\ \:\__\ /:/__/ \:\__\ /:/\:::::\__\         
 \:\  \  \/__/ \/_|::\/:/  / \/__\:\/:/  / \:\  \  \/__/ \/_|:|~~|~            
  \:\  \          |:|::/  /       \::/  /   \:\  \          |:|  |             
   \:\  \         |:|\/__/        /:/  /     \:\  \         |:|  |             
    \:\__\        |:|  |         /:/  /       \:\__\        |:|  |             
     \/__/         \|__|         \/__/         \/__/         \|__|             
      ___           ___           ___           ___                            
     /\__\         /\  \         /\  \         /\__\                           
    /:/ _/_       /::\  \       /::\  \       /:/  /                           
   /:/ /\__\     /:/\:\  \     /:/\:\  \     /:/__/                            
  /:/ /:/ _/_   /::\~\:\  \   /:/  \:\  \   /::\__\____                        
 /:/_/:/ /\__\ /:/\:\ \:\__\ /:/__/ \:\__\ /:/\:::::\__\                       
 \:\/:/ /:/  / \/__\:\/:/  / \:\  \  \/__/ \/_|:|~~|~                          
  \::/_/:/  /       \::/  /   \:\  \          |:|  |                           
   \:\/:/  /        /:/  /     \:\  \         |:|  |                           
    \::/  /        /:/  /       \:\__\        |:|  |                           
     \/__/         \/__/         \/__/         \|__|                           
       ___         ___           ___           ___                             
      /\  \       /\  \         /\  \         /\__\                            
      \:\  \     /::\  \       /::\  \       /:/  /                            
  ___ /::\__\   /:/\:\  \     /:/\:\  \     /:/__/                             
 /\  /:/\/__/  /::\~\:\  \   /:/  \:\  \   /::\__\____                         
 \:\/:/  /    /:/\:\ \:\__\ /:/__/ \:\__\ /:/\:::::\__\                        
  \::/  /     \/__\:\/:/  / \:\  \  \/__/ \/_|:|~~|~                           
   \/__/           \::/  /   \:\  \          |:|  |                            
                   /:/  /     \:\  \         |:|  |                            
                  /:/  /       \:\__\        |:|  |      
    ''',
    R'''
                                                                                                                                              
                                                                                                                        .,.   .(**(,                     
                                                                                                                ,(((((/.   .,. .   ,                     
                                                                                                  .    (,                      ***./                     
                                                                                        ///.   (.    .                        .**.  #                    
                                                                                      /,      .                          (.          (                   
                                                                                     /....*. ****.          .(.      @ .*%*. *,.     (                   
                                                                               .(..........( .****.          (      /#       *%      *.                  
                                                                               (......,(((,*.  ***,     ,/@(...,%%,#           .%.    ,                  
                                                                              (.......((((..,        , (.          &     &///&  ,.     (                 
                                                                            (........*(((...,(..**.   @     %#///&  .,  .(/@@/@  ,     (                 
                                                                             ,*................(.**  &     .(/@@%/%  #   %/((/@ %%/    (                 
                                                                               (...............*,    @     .(/@@(/&  (.   #@,    @   ,                 
                                                                               (....,((........./    /,     ,&///@   @@, .&*      &.    (                
                                                                                (.....(/........(     /.            &.*#.    ..(@.      ,                
                                                                                  (,............,,     .@         #(   ###(,.  .#@&(    (                
   / @ ( %@                            .@  /&           @          * @ (            /.............*,       ./(//.          . ./@&* ..,. (                
    ##&  %@.&%    & &/  *&%.    .%  .&@& /& %%*       &  #*.%%.   ##&             (........((....(.#. .        ../(&&@%&&      @ .***..                
         %@   #@  @,  .@   .@  @,   @. ,@  /@   (@      @  &@                     /.......((((..*,   .(.  ..%&**       @ @  @.%%.,***,.                
         %&    @  @   %@/////, %@@&(@. ,@  /&   *@      @  &(   %&                   (*....,((((../        #.     /    ,@/@.*( .%@  *** /                
         %@.  @#  @    @.  *@ ,@   ,@. ,@  /&   *@      @  &(   %&                     (....(((.../               @@%, #.@ ,, .@..&.    (                
                                                                                       /.....,....,/  .**.    .*/.../..%#../  & ,*%..**.,(.              
                                                                                      .,............(            .  ..@ . ...@.(* & .**  (  (&           
                                                                                       ,*...........*               .(        .*  &.    *.     #         
                                                                                         *.........(    .***        .%       ./   &&    *@.@ /,          
                                                                                        *...((,..../.  .***.           @ ,,/&@&@#@.. /@,*.@..&           
                                                                                        *...,,..&     *@(*.          .(. #@.  &  .@, . .(&.  &,        
                                                                                        /....../.          %*. .,&   @( %@%@   .      ,%@&,      #       
                                                                                           *&...&         ,  &     (,# .%###@.(&@#/////((.    */&%       
                                                                                           .,....#      %,&%,(      %..@####@(&@@%/@@@@%/.               
                                                                                            (//(%@@    ,(    ..#@#/& .&@%##@@@/////////*               
                                                         (((. ,(,...(                       #//@@&(/@, &  @#@#.   . .&@%(&@///////(%@(.                  
                                  ...       .(..../((/,........******.../                   %////(@@@(/&//@@@@@&(/(&@/////(&@%/////(                     
                     ,/////...//........*//.*******....,,*****,,.........,/                   /@#//////%///////////#@&/  /&(/////@                       
              ((............****,..........................,................(                      /@(/%////%@&((           ,.@                          
             (...,(((*./*..........*********,..........*/. .,(**(*.   ... .  .,                     ,///////////@           # &                          
             (............***,./(../(((*..,/(,.   .*/*.         .           ,(                       (@%#((##&@.            & &                          
              /..................*. .                               .**.  .*                              & &               (/(                          
              *................../ .**.    .                 .      ***   (                               @&@              */ &                          
             (............*((((.../****     %   ./.      @   @      ..   ./                               &/#            @@@#,@%  .@@@@@@(               
             /..........*(((((...( ****,(.%&,  .,@(  .....#@(.,/.     .(                                @ &            @@@@@@@@@@@(%@@@@@              
             ,..........(((/..../  .***.%           @%         .@     ,/                                  ( &            *@@@@@@@@@@@@@@@@@              
             ,..................(     /.      %%#&/  //    ,&%   @  ./                                 &@@#(&    *@@@@&  *@@@@.@@@@@@@@@@,               
            *,...,(*.............(.  ,/     ,%/#@(/@  *  *%/#%/#  * ..                                 @@@@@@@@@@@@@@@@@@           ..                   
            ,....((............./.** %      @//@@&/@   * @/#@@/@ .,  (                                 &@@@@@@@@@..#@@@@@                                
            ...................*.,**.,(     /(////&.  (  @////&  @  /                                  ,@@@&@@@@@@@@@@@&                                 
            /...................  ,   /.      @@@/   **    @&, /% /.                                       (    *&@@&.                                   
            /..................(      ..%           @../@#((&&.  (                                                                                       
            (..................(     ...  *@#,.,#@,  .,@,.       .*                                                                                      
           *.........../(((...(   ,..,(              /&   .,*.  /./                                                                                      
           (.........,(((((..(                         ,% .    .#(                                                                                       
            /.......,((((,...*.        ..*(%#/*..        .#%%#/,,                                                                                        
            (........((,......(     .%#..  .,     .#,.    /    ,,                                                                                        
           (.................*.            &      **,...,/(*    *                                                                                        
           *................(.   .,.      . .     .       ./   ./                                                                                        
           .,...............(.  **.      .,/,...//.,,,,,*,...  (                                                                                         
            (..((/.........../                           .**.. ,.,@.                                                                                     
            ...,............./  .,*                      ,(((((( *  #                                                                                    
           (.......@    *@%./,,***,             .(***(@       *@ @  @                                (@               &(                                 
           ,,.....#.         .@/*.  .*((/. .*(((&. @ #       .*#@.,@                                 (@ //     ./(    (*                                 
              &(..,#            **   .@  @  %(#(&    ,%@&(////%.&(, %                                (@.  (@  @   ,@                                   
              /....*.       %(.,*%&@@%%*,  @(((((@///#(/@@@@@%&    %@..@.                            (@    @ #&    %%                                  
              //%@#**,     @/*.. ..            ,,(#&@@@%#/////&   .%    *.&/...                      (@.  #@  @   ,@                                   
              ,///////&,   ,  %////##/#%@@&&((..     .         ..**.  ./ *& ...@                     ,( /*     .//    /*                                 
               (//@@@@@@%&@%@@@@@@@@@/////////&////#@@##*/%&@@#(*.      %,  @&@/                                                                         
               @%//////#%//(%/////////////&@&(///////&          ((  #, .%, .(#                                                                         
                   /@(/////%////////(@@/       %@@&&@%              %.  &.  %%  ,#                                                                       
                        %@#&(&@%////(           ( ,(                  ,@. .& . ,@@                                                                       
                        %////////////           *..@                     %@,%&@**(                                                                       
                         ./(&@@#@.               /*#                                                                                                     
                             *. @               *,.&                                                                                                     
                              &&@              #&  %                                                                                                     
                              * &            &@@@@@@@  &@@@@@@@                                                                                          
                              , @            @@@@@@@@@@@%#%@@@@@#                                                                                        
                             #  @             @@@@@@@@@@@@@@@@@@,                                                                                        
                          ,@@@(&%    (@@@@@   @@@@**@@@@@@@@@@#                                                                                          
                          @@@@@@@@@@@@@@@@@@@           ,**,                                                                                             
                          *@@@@@@@@@#  *@@@@@*                                                                                                           
                           @@@@&@@@@@@@@@@@@,                                                                                                            
                            .#@&& (@@@@@@@,                                                                                                              
                                
    ''',
    "There is a great disturbance in the force (The Emperor's Strikes Back) - GROMACS",
    "Take me to SCIENCE TOWN",
    "If you believe science, you can achieve science - Alfert Einstein",
    "It do be like that",
    "YEEEEEET",
    "Water is a fruit",
    "You cannot have watermelon on your top ten fruit list - Daniel Bim",
    "QMDMDMDMDMDMDMDMDMDMDMDMDMD - various",
    "Please call Hoffmanisburningsopleaseresubmit.py",
    "Lava their house cause I'm  a griefer (Im a griefer Im a griefer baby) - Minecraftcito by Zerina aka Reptilelegit",
    "WHEN THE MOON HITS YOUR EYE LIKE A BIG PIZZZZZZA-PIE, THATS AMOREEEEEEEEEEEEEE - Papa Franku"
]
