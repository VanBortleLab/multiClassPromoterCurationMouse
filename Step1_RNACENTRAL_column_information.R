require(here)
setwd(here::here())
###############IMPORTANT_COLUMN_FOR_EACH_DATABASE_(writing manually)
###########Please put the column number you want to extract from each databases .
########### Name of variables below should be same as the database name in grch38.bed file(case sensitive)
###########If there is no important column, put c()
#First element = URSID
#Second element = RNA type
#Third element  = RNA information
imp_pirbase=c(1,5,3)
#urs0000497d67  pirbase pir-mmu-27954321        10090   mirna   
#urs0000497d6d  pirbase pir-mmu-49702659        10090   pirna
imp_genecards=c(1,5,6)   #has info
#urs0000342bf5  genecards       genecards:pir-61088:urs0000342bf5_9606  9606    pirna   pir-61088
#urs0000342c87  genecards       genecards:linc01405:urs0000342c87_9606  9606    lncrna  linc01405
imp_lncipedia=c(1,5,3)
#urs0000d5b14f  lncipedia       lnc-large1-1:2  9606    lncrna  lnc-large1-1
#urs0000d5b150  lncipedia       bves-as1:9      9606    lncrna  bves-as1
imp_noncode=c(1,5,6)
#urs000075fbae  noncode nonmmut013074.2 10090   lncrna  nonmmug008288.2
#urs000075fbc9  noncode nonmmut012771.2 10090   antisense_rna   nonmmug008111.2
imp_ena=c(1,5,6)
#urs0000022706#ena#jf127899.1:1..1353:rrna#77133#rrna#
#urs0000022707#ena#dq586838.1:1..31:ncrna#9606#pirna#pir-53950
#urs0000022708#ena#jf190256.1:1..1375:rrna#77133#rrna#
#urs0000022709#ena#ay695528.1:370..533:rrna#306495#rrna#
#urs000002270a#ena#hq412925.1:464..553:trna#1037421#trna#trnl
#urs000002270a#ena#hq412928.1:471..560:trna#218143#trna#trnl
imp_malacards=c(1,5,6)
#urs00000011df  malacards       malacards:mir660:urs00000011df_9606     9606    mirna   mir660
#urs00000012ec  malacards       malacards:mt-ti:urs00000012ec_9606      9606    trna    mt-ti
imp_lncbook=c(1,5,3)
#urs0001bd2402  lncbook hsalnt0226895   9606    trna    hsalng0109131
#urs0001bd2403  lncbook hsalnt0172624   9606    lncrna  hsalng0082949
imp_mirbase=c(1,5,3)
#urs000075badc  mirbase mimat0018569    7460    mirna   
#urs000075badd  mirbase mi0016863       9606    pre_mirna
imp_refseq=c(1,5,6)
#urs00001bcdb7  refseq  nr_050243       6239    pirna   21ur-15127
#urs00001bce76  refseq  nr_051674       6239    srna    t06d8.16
imp_ensembl=c(1,5,6)
#urs0000aa14e0  ensembl ensppyt00000024603      9601    snorna  ensppyg00000021124.1
#urs0000aa14e4  ensembl enscaft00020034729      286419  rrna    enscafg00020023483.1
imp_ensembl_gencode=c(1,5,6)
#urs0000208a18  gencode enst00000413954 9606    lncrna  ensg00000224643.5
#urs0000208a34  gencode enst00000384906 9606    pre_mirna       mir99a
imp_hgnc=c(1,5,3)
#urs0000063a17  hgnc    hgnc:55765      9606    lncrna  
#urs0000064325  hgnc    hgnc:32613      9606    snorna  
imp_lncbase=c(1,5,3)
#urs000014d9c1  lncbase hsa-mir-200b-3p 9606    mirna   
#urs000014d9c1  lncbase mmu-mir-200b-3p 10090   pirna
imp_rfam=c(1,5,3)
#urs0000be993e  rfam    rf01695 1792290 antisense_rna   
#urs0000be993f  rfam    rf02096 4537    pre_mirna 
#urs0000be9940  rfam    rf00005 9940    trna    
#urs0000be9941  rfam    rf02541 1759557 ncrna
imp_tarbase=c(1,5,3)
#urs00004a6e26  tarbase hsa-mir-3659    9606    mirna   
imp_mirgenedb=c(1,5,3)
#urs0000efb973  mirgenedb       lan-mir-190-p8  7574    precursor_rna   
#urs0000efb974  mirgenedb       efe-mir-133     6396    pre_mirna
imp_intact=c(1, 5, 3)
#urs00004f2bde  intact  intact:urs00004f2bde_9606       9606    snorna  
#urs00004f859b  intact  intact:urs00004f859b_9606       9606    mirna
imp_psicquic=c(1, 5, 3)
#urs000075b8b8#psicquic#psicquic:urs000075b8b8_9606#9606#mirna#
#urs000075bb59#psicquic#psicquic:urs000075bb59_9606#9606#mirna#
imp_snodb=c(1, 5, 3)
#urs00004acfcf  snodb   189     9606    snorna  
imp_snopy=c(1,5,6) 
#urs00006bd550  snopy   danio_rerio300130       7955    snorna  snord15
imp_pdbe=c(1,5,3)
#urs000080de1d  pdb     2lbs_a  32630   misc_rna        
#urs000080de1e  pdb     1l9a_b  9606    srp_rna
imp_gtrnadb=c(1,5,6)
#urs00003ca2ea  gtrnadb gtrnadb:trna-leu-taa-1-2:hg530135.1:159414-159504       1231072 trna    trna-leu-taa-1-2
#urs00003ca382  gtrnadb gtrnadb:trna-glu-ttc-1-1:ch408076.1:114386-114457       306902  trna    trna-glu-ttc-1-1
imp_modomics=c(1,5,3)
#urs00005f0912#modomics#trna#562#trna#
#urs00005f2c2d#modomics#ssu#4932#rrna#
imp_lncrnadb=c(1,5,6)
#urs000019b796  lncrnadb        146     9606    lncrna  lust
#urs00001a3907  lncrnadb        102     10116   srna    khps1a
imp_crw=c(1,5,3)
#urs0001bc209d  crw     crw:d.5.e.s.aggregatum  4773    rrna    
#urs0001bc3a50  crw     crw:d.16.e.s.lophii     51541   rrna
imp_5srrnadb=c(1,5,3)
#urs0000f00d53#5srrnadb#b06299#1449070#rrna#
#urs0000f00d54#5srrnadb#b03555#1515746#rrna#
imp_silva=c(1, 5, 3)
#urs00009ddb5b  silva   silva:lrif02000002.1:complement(935006..936566) 669     rrna
imp_srpdb=c(1,5,3)
#urs0000069949  srpdb   baci.cere._cp000764     315749  srp_rna srp rna
#urs000006c5ed  srpdb   meth.mari._cp000609.1   402880  srp_rna srp rna
#urs000007022e  srpdb   stre.sang._cp000387     388919  srna    srp rna
imp_evlncrnas=c(1,5,3)
# urs000000f15c evlncrnas       evlncrnas:el1163        9606    lncrna  
# urs0000010eb4 evlncrnas       evlncrnas:el2452        9606    lncrna  
# urs0000021c91 evlncrnas       evlncrnas:el0929        9606    lncrna  
imp_ribocentre=c(1,5,3)
# urs0000034135 ribocentre      seq211045       582514  rnase_p_rna     
# urs0000034a2e ribocentre      seq014259       188913  ncrna   
# urs00000368f8 ribocentre      seq010876       33032   ncrna   
# urs0000036f85 ribocentre      seq218160       45790   rnase_p_rna     
imp_ribovision=c(1,5,3)
# urs0000005270 ribovision      ribovision:4v6x_a8      9606    rrna    
# urs000000ee76 ribovision      ribovision:4v6u_b3      186497  rrna    
# urs00000f9d45 ribovision      ribovision:4v6x_a7      9606    rrna
imp_expression_atlas=c(1, 5, 3)
imp_mgi=c(1, 5, 3) #CHECK
