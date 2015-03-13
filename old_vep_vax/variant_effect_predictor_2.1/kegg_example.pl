#!/usr/bin/perl -w
use strict;


use SOAP::Lite;

sub SOAP::Serializer::as_ArrayOfstring{
  my ($self, $value, $name, $type, $attr) = @_;
  return [$name, {'xsi:type' => 'array', %$attr}, $value];
}

sub SOAP::Serializer::as_ArrayOfint{
  my ($self, $value, $name, $type, $attr) = @_;
  return [$name, {'xsi:type' => 'array', %$attr}, $value];
}

my $wsdl = 'http://soap.genome.jp/KEGG.wsdl';
open OUT,">","/Users/myourshaw/lab/pypeline/variant_effect_predictor_2.1/kegg_pathways.txt";
my $paths = SOAP::Lite
             -> service($wsdl)
             -> list_pathways("hsa");
foreach my $path (@{$paths}) {
    my $genes = SOAP::Lite
             -> service($wsdl)
             -> get_genes_by_pathway($path->{entry_id});
    foreach my $gene(@{$genes}){
        print OUT "$path->{entry_id}\t$path->{definition}\t$gene\n";
    }
}
close OUT
#my $pathsT = SOAP::Lite
#         -> service($wsdl)
#         -> get_pathways_by_genes(['hsa:10327' , 'hsa:124']);
#foreach my $pathT (@{$pathsT}){
#    print;
#}
#
#my $genesT = SOAP::Lite
#         -> service($wsdl)
#         -> get_genes_by_pathway('path:hsa00010');
#foreach my $geneT(@{$genesT}){
#   print  "$geneT\n";
#}
#hsa:10327
#hsa:124
#hsa:125
#hsa:126
#hsa:127
#hsa:128
#hsa:130
#hsa:130589
#hsa:131
#hsa:160287
#hsa:1737
#hsa:1738
#hsa:2023
#hsa:2026
#hsa:2027
#hsa:217
#hsa:218
#hsa:219
#hsa:220
#hsa:2203
#hsa:221
#hsa:222
#hsa:223
#hsa:224
#hsa:226
#hsa:229
#hsa:230
#hsa:2538
#hsa:2597
#hsa:26330
#hsa:2645
#hsa:2821
#hsa:3098
#hsa:3099
#hsa:3101
#hsa:3939
#hsa:3945
#hsa:3948
#hsa:441531
#hsa:501
#hsa:5105
#hsa:5106
#hsa:5160
#hsa:5161
#hsa:5162
#hsa:5211
#hsa:5213
#hsa:5214
#hsa:5223
#hsa:5224
#hsa:5230
#hsa:5232
#hsa:5236
#hsa:5313
#hsa:5315
#hsa:55276
#hsa:55902
#hsa:57818
#hsa:669
#hsa:7167
#hsa:80201
#hsa:83440
#hsa:84532
#hsa:8789
#hsa:92483


#my $pathse = SOAP::Lite
#         -> service($wsdl)
#         -> get_pathways_by_genes(['eco:b0077' , 'eco:b0078']);
#foreach my $pathe (@{$pathse}){
#    print;
#}
#my $limit = 100;
#for (my $offset=0;$offset<300;$offset+=100){
#    my $genes = SOAP::Lite
#             -> service($wsdl)
#             -> get_genes_by_organism("hsa", $offset. $limit);
#    if (scalar(@{$genes})==0){
#       last;
#    }
#    my $genesa = SOAP::Data->type(array => $genes);
#    my $paths = SOAP::Lite
#             -> service($wsdl)
#             -> get_pathways_by_genes($genesa);
#    
#    foreach my $path(@{$paths}){
#       print;
#    }
#    foreach my $gene (@{$genes}) {
#      print "$gene\n";
#    }
#}


#print; #hsa:285033
#my $dbs = SOAP::Lite
#             -> service($wsdl)
#             -> list_databases();
#foreach my $db (@{$dbs}) {
#  print "$db->{entry_id}\t$db->{definition}\n";
#}
#nt	Non-redundant nucleic acid sequence database
#aa	Non-redundant protein sequence database
#gb	GenBank nucleic acid sequence database
#gbu	Cumulative daily updates of GenBank since the latest release
#rs	NCBI Reference Sequences
#rs	NCBI Reference Sequences
#rs	NCBI Reference Sequences
#emb	EMBL nucleic acid sequence database
#embu	Cumulative daily updates of EMBL since the latest release
#est	dbEST database of Expressed Sequence Tags
#gss	dbGSS database of Genome Survey Sequences
#sts	dbSTS database of Sequence Taged Sites
#htgs	HTGs database of High Throughput Genomic Sequences
#sp	SWISS-PROT protein sequence database
#tr	TrEMBL protein sequence database
#prf	PRF protein sequence database
#gp	Translated protein sequences from genbank
#gpu	Translated protein sequences from genbank-upd
#pdb	RCSB Protein Data Bank
#pdbstr	Re-organized Protein Data Bank
#epd	Eukaryotic Promoter Database
#ps	Dictionary of Protein Sites and Patterns
#pdoc	PROSITE ducumentation file
#bl	Blocks Database
#pr	Protein Motif Fingerprint Database
#pd	ProDom protein domain database
#pf	Protein families database of alignments and HMMs
#pmd	Protein Mutant Database
#aax	Amino Acid Index Database, Combination of AAindex1, AAindex2 and AAindex3
#ex	KEGG Expression Database
#lit	Literature Database compiled by PRF
#omim	Online Mendelian Inheritance in Man
#path	KEGG Pathway Database
#md	KEGG Module Database
#ds	KEGG Disease Database
#br	KEGG Brite Database
#ko	KEGG Orthology Database
#gn	KEGG Genome Database
#mgnm	KEGG Meta Genome Database
#genes	KEGG Genes Database
#dg	KEGG Draft/Partial Genomes Genes Database
#eg	KEGG EST Contigs Genes Database
#mg	KEGG Metagenome Genes Database
#vg	KEGG Virus Genes Database
#og	KEGG Organelle Genes Database
#cpd	KEGG Compound Database
#dr	KEGG Drug Database
#ev	KEGG Environ Database
#gl	KEGG Glycan Database
#rn	KEGG Reaction Database
#rp	KEGG Reactant Pair Database
#rc	KEGG Reactant Class Database
#ec	KEGG Enzyme Database
#ld	Database of Link Information

#my $orgs = SOAP::Lite
#             -> service($wsdl)
#             -> list_organisms();
#foreach my $org (@{$orgs}) {
#  print "$org->{entry_id}\t$org->{definition}\n";
#}
#hsa	Homo sapiens (human)
#ptr	Pan troglodytes (chimpanzee)
#pon	Pongo abelii (Sumatran orangutan)
#mcc	Macaca mulatta (rhesus monkey)
#mmu	Mus musculus (mouse)
#rno	Rattus norvegicus (rat)
#cfa	Canis familiaris (dog)
#aml	Ailuropoda melanoleuca (giant panda)
#bta	Bos taurus (cow)
#ssc	Sus scrofa (pig)
#ecb	Equus caballus (horse)
#mdo	Monodelphis domestica (opossum)
#oaa	Ornithorhynchus anatinus (platypus)
#gga	Gallus gallus (chicken)
#mgp	Meleagris gallopavo (turkey)
#tgu	Taeniopygia guttata (zebra finch)
#xla	Xenopus laevis (African clawed frog)
#xtr	Xenopus tropicalis (western clawed frog)
#dre	Danio rerio (zebrafish)
#bfo	Branchiostoma floridae (Florida lancelet)
#cin	Ciona intestinalis (sea squirt)
#spu	Strongylocentrotus purpuratus (purple sea urchin)
#dme	Drosophila melanogaster (fruit fly)
#dpo	Drosophila pseudoobscura pseudoobscura
#dan	Drosophila ananassae
#der	Drosophila erecta
#dpe	Drosophila persimilis
#dse	Drosophila sechellia
#dsi	Drosophila simulans
#dwi	Drosophila willistoni
#dya	Drosophila yakuba
#dgr	Drosophila grimshawi
#dmo	Drosophila mojavensis
#dvi	Drosophila virilis
#aga	Anopheles gambiae (mosquito)
#aag	Aedes aegypti (yellow fever mosquito)
#cqu	Culex quinquefasciatus (southern house mosquito)
#ame	Apis mellifera (honey bee)
#nvi	Nasonia vitripennis (jewel wasp)
#tca	Tribolium castaneum (red flour beetle)
#api	Acyrthosiphon pisum (pea aphid)
#phu	Pediculus humanus corporis (human body louse)
#isc	Ixodes scapularis (black-legged tick)
#cel	Caenorhabditis elegans (nematode)
#cbr	Caenorhabditis briggsae
#bmy	Brugia malayi (filaria)
#smm	Schistosoma mansoni
#nve	Nematostella vectensis (sea anemone)
#hmg	Hydra magnipapillata
#tad	Trichoplax adhaerens
#ath	Arabidopsis thaliana (thale cress)
#aly	Arabidopsis lyrata (lyrate rockcress)
#pop	Populus trichocarpa (black cottonwood)
#rcu	Ricinus communis (castor bean)
#vvi	Vitis vinifera (wine grape)
#osa	Oryza sativa japonica (Japanese rice)
#sbi	Sorghum bicolor (sorghum)
#zma	Zea mays (maize)
#smo	Selaginella moellendorffii
#ppp	Physcomitrella patens subsp. patens
#cre	Chlamydomonas reinhardtii
#vcn	Volvox carteri f. nagariensis
#olu	Ostreococcus lucimarinus
#ota	Ostreococcus tauri
#cme	Cyanidioschyzon merolae
#sce	Saccharomyces cerevisiae (budding yeast)
#ago	Ashbya gossypii (Eremothecium gossypii)
#kla	Kluyveromyces lactis
#lth	Lachancea thermotolerans
#ppa	Pichia pastoris
#vpo	Vanderwaltozyma polyspora
#zro	Zygosaccharomyces rouxii
#cgr	Candida glabrata
#dha	Debaryomyces hansenii
#pic	Scheffersomyces stipitis
#pgu	Meyerozyma guilliermondii
#lel	Lodderomyces elongisporus
#cal	Candida albicans
#ctp	Candida tropicalis
#cdu	Candida dubliniensis
#yli	Yarrowia lipolytica
#clu	Clavispora lusitaniae
#ncr	Neurospora crassa
#pan	Podospora anserina
#mgr	Magnaporthe oryzae
#fgr	Fusarium graminearum
#ssl	Sclerotinia sclerotiorum
#bfu	Botryotinia fuckeliana
#ani	Aspergillus nidulans
#afm	Aspergillus fumigatus
#nfi	Neosartorya fischeri
#aor	Aspergillus oryzae
#ang	Aspergillus niger
#afv	Aspergillus flavus
#act	Aspergillus clavatus
#pcs	Penicillium chrysogenum
#cim	Coccidioides immitis
#cpw	Coccidioides posadasii
#ure	Uncinocarpus reesii
#pno	Phaeosphaeria nodorum
#tml	Tuber melanosporum
#spo	Schizosaccharomyces pombe (fission yeast)
#cne	Cryptococcus neoformans JEC21
#cnb	Cryptococcus neoformans B-3501A
#ppl	Postia placenta
#lbc	Laccaria bicolor
#mpr	Moniliophthora perniciosa
#cci	Coprinopsis cinerea
#scm	Schizophyllum commune
#uma	Ustilago maydis
#mgl	Malassezia globosa
#ecu	Encephalitozoon cuniculi
#mbr	Monosiga brevicollis
#ngr	Naegleria gruberi
#ddi	Dictyostelium discoideum (cellular slime mold)
#ehi	Entamoeba histolytica
#edi	Entamoeba dispar
#pfa	Plasmodium falciparum 3D7
#pfd	Plasmodium falciparum Dd2
#pfh	Plasmodium falciparum HB3
#pyo	Plasmodium yoelii
#pcb	Plasmodium chabaudi
#pbe	Plasmodium berghei
#pkn	Plasmodium knowlesi
#pvx	Plasmodium vivax
#tan	Theileria annulata
#tpv	Theileria parva
#bbo	Babesia bovis
#cpv	Cryptosporidium parvum
#cho	Cryptosporidium hominis
#tgo	Toxoplasma gondii
#tet	Tetrahymena thermophila
#ptm	Paramecium tetraurelia
#tbr	Trypanosoma brucei
#tcr	Trypanosoma cruzi
#lma	Leishmania major
#lif	Leishmania infantum
#lbz	Leishmania braziliensis
#gla	Giardia lamblia
#tva	Trichomonas vaginalis
#pti	Phaeodactylum tricornutum
#tps	Thalassiosira pseudonana
#pif	Phytophthora infestans
#eco	Escherichia coli K-12 MG1655
#ecj	Escherichia coli K-12 W3110
#ecd	Escherichia coli K-12 DH10B
#ebw	Escherichia coli K-12 MC4100(MuLac) BW2952
#ece	Escherichia coli O157 H7 EDL933 (EHEC)
#ecs	Escherichia coli O157 H7 Sakai (EHEC)
#ecf	Escherichia coli O157 H7 EC4115 (EHEC)
#etw	Escherichia coli O157 H7 TW14359 (EHEC)
#eoj	Escherichia coli O26 H11 11368 (EHEC)
#eoi	Escherichia coli O111 H- 11128 (EHEC)
#eoh	Escherichia coli O103 H2 12009 (EHEC)
#ecg	Escherichia coli O127 H6 E2348/69 (EPEC)
#eok	Escherichia coli O55 H7 CB9615 (EPEC)
#ecc	Escherichia coli O6 K2 H1 CFT073 (UPEC)
#ecp	Escherichia coli O6 K15 H31 536 (UPEC)
#eci	Escherichia coli O18 K1 H7 UTI89 (UPEC)
#ecv	Escherichia coli O1 K1 H7 (APEC)
#ecx	Escherichia coli O9 HS (commensal)
#ecw	Escherichia coli O139 H28 E24377A (ETEC)
#ecm	Escherichia coli SMS-3-5 (environmental)
#ecy	Escherichia coli O152 H28 SE11 (commensal)
#ecr	Escherichia coli O8 IAI1 (commensal)
#ecq	Escherichia coli O81 ED1a (commensal)
#eck	Escherichia coli 55989 (EAEC)
#ect	Escherichia coli O7 K1 IAI39 (ExPEC)
#eum	Escherichia coli O17 K52 H18 UMN026 (ExPEC)
#ecz	Escherichia coli O45 K1 H7 S88 (ExPEC)
#ecl	Escherichia coli C ATCC 8739
#ebr	Escherichia coli B REL606
#ebd	Escherichia coli BL21-Gold(DE3)pLysS AG
#efe	Escherichia fergusonii
#sty	Salmonella enterica subsp. enterica serovar Typhi CT18
#stt	Salmonella enterica subsp. enterica serovar Typhi Ty2
#stm	Salmonella enterica subsp. enterica serovar Typhimurium LT2
#spt	Salmonella enterica subsp. enterica serovar Paratyphi A ATCC9150
#sek	Salmonella enterica subsp. enterica serovar Paratyphi A AKU12601
#spq	Salmonella enterica subsp. enterica serovar Paratyphi B
#sei	Salmonella enterica subsp. enterica serovar Paratyphi C
#sec	Salmonella enterica subsp. enterica serovar Choleraesuis
#seh	Salmonella enterica subsp. enterica serovar Heidelberg
#see	Salmonella enterica subsp. enterica serovar Newport
#sew	Salmonella enterica subsp. enterica serovar Schwarzengrund
#sea	Salmonella enterica subsp. enterica serovar Agona
#sed	Salmonella enterica subsp. enterica serovar Dublin
#seg	Salmonella enterica subsp. enterica serovar Gallinarum
#set	Salmonella enterica subsp. enterica serovar Enteritidis
#ses	Salmonella enterica subsp. arizonae
#sbg	Salmonella bongori
#ype	Yersinia pestis CO92 (biovar Orientalis)
#ypk	Yersinia pestis KIM 10 (biovar Mediaevalis)
#ypa	Yersinia pestis Antiqua (biovar Antiqua)
#ypn	Yersinia pestis Nepal516 (biovar Antiqua)
#ypm	Yersinia pestis 91001 (biovar Microtus)
#ypp	Yersinia pestis Pestoides F
#ypg	Yersinia pestis Angola
#ypz	Yersinia pestis Z176003
#yps	Yersinia pseudotuberculosis IP32953 (serotype I)
#ypi	Yersinia pseudotuberculosis IP31758 (serotype O 1b)
#ypy	Yersinia pseudotuberculosis YPIII
#ypb	Yersinia pseudotuberculosis PB1/+
#yen	Yersinia enterocolitica subsp. enterocolitica 8081
#yep	Yersinia enterocolitica subsp. palearctica 105.5R(r)
#sfl	Shigella flexneri 301 (serotype 2a)
#sfx	Shigella flexneri 2457T (serotype 2a)
#sfv	Shigella flexneri 8401 (serotype 5b)
#ssn	Shigella sonnei
#sbo	Shigella boydii Sb227
#sbc	Shigella boydii CDC 3083-94
#sdy	Shigella dysenteriae
#eca	Pectobacterium atrosepticum
#pct	Pectobacterium carotovorum
#pwa	Pectobacterium wasabiae
#eta	Erwinia tasmaniensis
#epy	Erwinia pyrifoliae
#eam	Erwinia amylovora CFBP1430
#eay	Erwinia amylovora ATCC 49946
#ebi	Erwinia billingiae
#plu	Photorhabdus luminescens
#pay	Photorhabdus asymbiotica
#buc	Buchnera aphidicola APS
#bas	Buchnera aphidicola Sg
#bab	Buchnera aphidicola Bp
#bcc	Buchnera aphidicola Cc
#bap	Buchnera aphidicola 5A
#bau	Buchnera aphidicola Tuc7
#baj	Buchnera aphidicola (Cinara tujafilina)
#wbr	Wigglesworthia glossinidia
#sgl	Sodalis glossinidius
#ent	Enterobacter sp. 638
#enc	Enterobacter cloacae subsp. cloacae ATCC 13047
#esc	Enterobacter cloacae SCF1
#eae	Enterobacter aerogenes
#esa	Cronobacter sakazakii
#ctu	Cronobacter turicensis
#kpn	Klebsiella pneumoniae
#kpe	Klebsiella pneumoniae 342
#kpu	Klebsiella pneumoniae NTUH-K2044
#kva	Klebsiella variicola
#cko	Citrobacter koseri ATCC BAA-895
#cro	Citrobacter rodentium
#spe	Serratia proteamaculans
#srs	Serratia sp. AS12
#srr	Serratia sp. AS9
#pmr	Proteus mirabilis
#eic	Edwardsiella ictaluri
#etr	Edwardsiella tarda
#bfl	Candidatus Blochmannia floridanus
#bpn	Candidatus Blochmannia pennsylvanicus
#bva	Candidatus Blochmannia vafer
#hde	Candidatus Hamiltonella defensa
#dda	Dickeya dadantii Ech703
#ddc	Dickeya dadantii Ech586
#ddd	Dickeya dadantii 3937
#dze	Dickeya zeae
#xbo	Xenorhabdus bovienii
#xne	Xenorhabdus nematophila
#pam	Pantoea ananatis
#pva	Pantoea vagans
#pao	Pantoea sp. At-9b
#rip	Candidatus Riesia pediculicola
#rah	Rahnella sp. Y9602
#men	Candidatus Moranella endobia
#hin	Haemophilus influenzae Rd KW20 (serotype d)
#hit	Haemophilus influenzae 86-028NP (nontypeable)
#hip	Haemophilus influenzae PittEE
#hiq	Haemophilus influenzae PittGG
#hif	Haemophilus influenzae F3031
#hil	Haemophilus influenzae F3047
#hdu	Haemophilus ducreyi
#hap	Haemophilus parasuis
#hso	Haemophilus somnus 129PT
#hsm	Haemophilus somnus 2336
#pmu	Pasteurella multocida
#msu	Mannheimia succiniciproducens
#apl	Actinobacillus pleuropneumoniae L20 (serotype 5b)
#apj	Actinobacillus pleuropneumoniae JL03 (serotype 3)
#apa	Actinobacillus pleuropneumoniae AP76 (serotype 7)
#asu	Actinobacillus succinogenes
#aap	Aggregatibacter aphrophilus
#aat	Aggregatibacter actinomycetemcomitans
#gan	Gallibacterium anatis
#xfa	Xylella fastidiosa 9a5c
#xft	Xylella fastidiosa Temecula1
#xfm	Xylella fastidiosa M12
#xfn	Xylella fastidiosa M23
#xcc	Xanthomonas campestris pv. campestris ATCC 33913
#xcb	Xanthomonas campestris pv. campestris 8004
#xca	Xanthomonas campestris pv. campestris B100
#xcv	Xanthomonas campestris pv. vesicatoria
#xac	Xanthomonas axonopodis
#xoo	Xanthomonas oryzae KACC10331
#xom	Xanthomonas oryzae MAFF311018
#xop	Xanthomonas oryzae PXO99A
#xal	Xanthomonas albilineans
#sml	Stenotrophomonas maltophilia K279a
#smt	Stenotrophomonas maltophilia R551-3
#psu	Pseudoxanthomonas suwonensis
#vch	Vibrio cholerae O1
#vco	Vibrio cholerae O395
#vcm	Vibrio cholerae M66-2
#vcj	Vibrio cholerae MJ-1236
#vvu	Vibrio vulnificus CMCP6
#vvy	Vibrio vulnificus YJ016
#vvm	Vibrio vulnificus MO6-24/O
#vpa	Vibrio parahaemolyticus
#vha	Vibrio harveyi
#vsp	Vibrio splendidus
#vex	Vibrio sp. Ex25
#vfi	Vibrio fischeri
#vfm	Vibrio fischeri MJ11
#vsa	Aliivibrio salmonicida LFI1238
#ppr	Photobacterium profundum
#van	Vibrio anguillarum
#pae	Pseudomonas aeruginosa PAO1
#pau	Pseudomonas aeruginosa UCBPP-PA14
#pap	Pseudomonas aeruginosa PA7
#pag	Pseudomonas aeruginosa LESB58
#ppu	Pseudomonas putida KT2440
#ppf	Pseudomonas putida F1
#ppg	Pseudomonas putida GB-1
#ppw	Pseudomonas putida W619
#ppt	Pseudomonas putida S16
#pst	Pseudomonas syringae pv. tomato DC3000
#psb	Pseudomonas syringae pv. syringae B728a
#psp	Pseudomonas syringae pv. phaseolicola 1448A
#pfl	Pseudomonas fluorescens Pf-5
#pfo	Pseudomonas fluorescens Pf0-1
#pfs	Pseudomonas fluorescens SBW25
#pen	Pseudomonas entomophila
#pmy	Pseudomonas mendocina ymp
#pmk	Pseudomonas mendocina NK-01
#psa	Pseudomonas stutzeri
#psz	Pseudomonas stutzeri ATCC 17588
#pba	Pseudomonas brassicacearum
#pfv	Pseudomonas fulva
#cja	Cellvibrio japonicus
#cga	Cellvibrio gilvus
#avn	Azotobacter vinelandii
#par	Psychrobacter arcticum
#pcr	Psychrobacter cryohalolentis
#prw	Psychrobacter sp. PRwf-1
#aci	Acinetobacter sp. ADP1
#acd	Acinetobacter sp. DR1
#acb	Acinetobacter baumannii ATCC 17978
#abm	Acinetobacter baumannii SDF
#aby	Acinetobacter baumannii AYE
#abc	Acinetobacter baumannii ACICU
#abn	Acinetobacter baumannii AB0057
#abb	Acinetobacter baumannii AB307-0294
#mct	Moraxella catarrhalis
#son	Shewanella oneidensis
#sdn	Shewanella denitrificans
#sfr	Shewanella frigidimarina
#saz	Shewanella amazonensis
#sbl	Shewanella baltica OS155
#sbm	Shewanella baltica OS185
#sbn	Shewanella baltica OS195
#sbp	Shewanella baltica OS223
#slo	Shewanella loihica
#spc	Shewanella putrefaciens
#sse	Shewanella sediminis
#spl	Shewanella pealeana
#she	Shewanella sp. MR-4
#shm	Shewanella sp. MR-7
#shn	Shewanella sp. ANA-3
#shw	Shewanella sp. W3-18-1
#shl	Shewanella halifaxensis
#swd	Shewanella woodyi ATCC 51908
#swp	Shewanella piezotolerans WP3
#svo	Shewanella violacea
#ilo	Idiomarina loihiensis
#cps	Colwellia psychrerythraea
#pha	Pseudoalteromonas haloplanktis
#pat	Pseudoalteromonas atlantica
#psm	Pseudoalteromonas sp. SM9913
#sde	Saccharophagus degradans
#maq	Marinobacter aquaeolei
#amc	Alteromonas macleodii
#alt	Alteromonas sp. SN2
#gag	Glaciecola sp. 4H-3-7+YE-5
#pin	Psychromonas ingrahamii
#ttu	Teredinibacter turnerae
#fbl	Ferrimonas balearica
#cbu	Coxiella burnetii RSA 493
#cbs	Coxiella burnetii RSA 331
#cbd	Coxiella burnetii Dugway 5J108-111
#cbg	Coxiella burnetii CbuG_Q212
#cbc	Coxiella burnetii CbuK_Q154
#lpn	Legionella pneumophila Philadelphia 1
#lpf	Legionella pneumophila Lens
#lpp	Legionella pneumophila Paris
#lpc	Legionella pneumophila Corby
#lpa	Legionella pneumophila 2300/99 Alcoy
#llo	Legionella longbeachae
#mca	Methylococcus capsulatus
#mmt	Methylomonas methanica
#ftu	Francisella tularensis subsp. tularensis SCHU S4
#ftf	Francisella tularensis subsp. tularensis FSC 198
#ftw	Francisella tularensis subsp. tularensis WY96-3418
#ftl	Francisella tularensis subsp. holarctica LVS
#fth	Francisella tularensis subsp. holarctica OSU18
#fta	Francisella tularensis subsp. holarctica FTNF002-00
#ftm	Francisella tularensis subsp. mediasiatica FSC147
#ftn	Francisella novicida U112
#fph	Francisella philomiragia
#frt	Francisella sp. TX077308
#tcx	Thiomicrospira crunogena
#tcy	Thioalkalimicrobium cyclicum
#noc	Nitrosococcus oceani
#nhl	Nitrosococcus halophilus
#nwa	Nitrosococcus watsonii
#alv	Allochromatium vinosum
#aeh	Alkalilimnicola ehrlichei
#hha	Halorhodospira halophila
#tgr	Thioalkalivibrio sp. HL-EbGR7
#tkm	Thioalkalivibrio sp. K90mix
#hna	Halothiobacillus neapolitanus
#hch	Hahella chejuensis
#csa	Chromohalobacter salexigens
#hel	Halomonas elongata
#abo	Alcanivorax borkumensis
#kko	Kangiella koreensis
#mmw	Marinomonas sp. MWYL1
#mme	Marinomonas mediterranea
#mpc	Marinomonas posidonica
#aha	Aeromonas hydrophila
#asa	Aeromonas salmonicida
#avr	Aeromonas veronii
#tau	Tolumonas auensis
#dno	Dichelobacter nodosus
#afe	Acidithiobacillus ferrooxidans ATCC 53993
#afr	Acidithiobacillus ferrooxidans ATCC 23270
#acu	Acidithiobacillus caldus
#afi	Acidithiobacillus ferrivorans
#bci	Baumannia cicadellinicola
#crp	Candidatus Carsonella ruddii
#rma	Candidatus Ruthia magnifica
#vok	Candidatus Vesicomyosocius okutanii
#gpb	Gamma proteobacterium HdN1
#nma	Neisseria meningitidis Z2491 (serogroup A)
#nme	Neisseria meningitidis MC58 (serogroup B)
#nmc	Neisseria meningitidis FAM18 (serogroup C)
#nmn	Neisseria meningitidis 053442 (serogroup C)
#nmi	Neisseria meningitidis alpha14
#ngo	Neisseria gonorrhoeae FA 1090
#ngk	Neisseria gonorrhoeae NCCP11945
#nla	Neisseria lactamica
#cvi	Chromobacterium violaceum
#lhk	Laribacter hongkongensis
#rso	Ralstonia solanacearum GMI1000
#rsc	Ralstonia solanacearum CFBP2957
#rsl	Ralstonia solanacearum PSI07
#rpi	Ralstonia pickettii 12J
#rpf	Ralstonia pickettii 12D
#reu	Ralstonia eutropha JMP134
#reh	Ralstonia eutropha H16
#rme	Cupriavidus metallidurans
#cti	Cupriavidus taiwanensis
#cnc	Cupriavidus necator
#bma	Burkholderia mallei ATCC 23344
#bmv	Burkholderia mallei SAVP1
#bml	Burkholderia mallei NCTC 10229
#bmn	Burkholderia mallei NCTC 10247
#bps	Burkholderia pseudomallei K96243
#bpm	Burkholderia pseudomallei 1710b
#bpl	Burkholderia pseudomallei 1106a
#bpd	Burkholderia pseudomallei 668
#bpr	Burkholderia pseudomallei MSHR346
#bte	Burkholderia thailandensis
#bvi	Burkholderia vietnamiensis
#bur	Burkholderia sp. 383
#bcn	Burkholderia cenocepacia AU1054
#bch	Burkholderia cenocepacia HI2424
#bcm	Burkholderia cenocepacia MC0-3
#bcj	Burkholderia cenocepacia J2315
#bam	Burkholderia cepacia
#bac	Burkholderia ambifaria MC40-6
#bmu	Burkholderia multivorans ATCC 17616 (JGI)
#bmj	Burkholderia multivorans ATCC 17616 (Tohoku)
#bxe	Burkholderia xenovorans
#bph	Burkholderia phymatum
#bpy	Burkholderia phytofirmans
#bgl	Burkholderia glumae
#bug	Burkholderia sp. CCGE1001
#bge	Burkholderia sp. CCGE1002
#bgf	Burkholderia sp. CCGE1003
#brh	Burkholderia rhizoxinica
#bgd	Burkholderia gladioli
#buj	Burkholderia sp. JV3
#pnu	Polynucleobacter sp. QLW-P1DMWA-1
#pne	Polynucleobacter necessarius
#bpe	Bordetella pertussis
#bpa	Bordetella parapertussis
#bbr	Bordetella bronchiseptica
#bpt	Bordetella petrii
#bav	Bordetella avium
#axy	Achromobacter xylosoxidans
#teq	Taylorella equigenitalis
#put	Pusillimonas sp. T7-7
#rfr	Rhodoferax ferrireducens
#pol	Polaromonas sp. JS666
#pna	Polaromonas naphthalenivorans
#aav	Acidovorax avenae
#ajs	Acidovorax sp. JS42
#dia	Acidovorax ebreus
#aaa	Acidovorax avenae subsp. avenae ATCC 19860
#vei	Verminephrobacter eiseniae
#dac	Delftia acidovorans
#del	Delftia sp. Cs1-4
#vap	Variovorax paradoxus S110
#vpe	Variovorax paradoxus EPS
#ctt	Comamonas testosteroni
#adn	Alicycliphilus denitrificans BC
#adk	Alicycliphilus denitrificans K601
#rta	Ramlibacter tataouinensis
#mpt	Methylibium petroleiphilum
#har	Herminiimonas arsenicoxydans
#mms	Minibacterium massiliensis
#hse	Herbaspirillum seropedicae
#zin	Candidatus Zinderia insecticola CARI
#cfu	Collimonas fungivorans
#lch	Leptothrix cholodnii
#tin	Thiomonas intermedia
#neu	Nitrosomonas europaea
#net	Nitrosomonas eutropha
#nit	Nitrosomonas sp. AL212
#nii	Nitrosomonas sp. Is79A3
#nmu	Nitrosospira multiformis
#eba	Aromatoleum aromaticum EbN1
#azo	Azoarcus sp. BH72
#dar	Dechloromonas aromatica
#tmz	Thauera sp. MZ1T
#tbd	Thiobacillus denitrificans
#mfa	Methylobacillus flagellatus
#mmb	Methylotenera mobilis
#meh	Methylotenera sp. 301
#mei	Methylovorus sp. SIP3-4
#mep	Methylovorus sp. MP688
#app	Accumulibacter phosphatis
#tpn	Candidatus Tremblaya princeps
#slt	Sideroxydans lithotrophicus
#gca	Gallionella capsiferriformans
#hpy	Helicobacter pylori 26695
#hpj	Helicobacter pylori J99
#hpa	Helicobacter pylori HPAG1
#hps	Helicobacter pylori Shi470
#hpg	Helicobacter pylori G27
#hpp	Helicobacter pylori P12
#hpb	Helicobacter pylori B38
#hpl	Helicobacter pylori B8
#hpc	Helicobacter pylori PeCan4
#hpm	Helicobacter pylori SJM180
#hhe	Helicobacter hepaticus
#hac	Helicobacter acinonychis
#hms	Helicobacter mustelae
#hfe	Helicobacter felis
#hbi	Helicobacter bizzozeronii
#wsu	Wolinella succinogenes
#tdn	Sulfurimonas denitrificans
#sua	Sulfurimonas autotrophica
#sku	Sulfuricurvum kujiense
#cje	Campylobacter jejuni NCTC11168
#cjr	Campylobacter jejuni RM1221
#cjj	Campylobacter jejuni 81-176
#cju	Campylobacter jejuni 81116
#cjn	Campylobacter jejuni ICDCCJ07001
#cjd	Campylobacter jejuni subsp. doylei 269.97
#cff	Campylobacter fetus
#ccv	Campylobacter curvus
#cha	Campylobacter hominis ATCC BAA-381
#cco	Campylobacter concisus 13826
#cla	Campylobacter lari
#abu	Arcobacter butzleri
#ant	Arcobacter nitrofigilis
#sdl	Sulfurospirillum deleyianum
#nis	Nitratiruptor sp. SB155-2
#sun	Sulfurovum sp. NBC37-1
#nsa	Nitratifractor salsuginis
#nam	Nautilia profundicola
#gsu	Geobacter sulfurreducens
#gme	Geobacter metallireducens
#gur	Geobacter uraniumreducens
#glo	Geobacter lovleyi
#gbm	Geobacter bemidjiensis
#geo	Geobacter sp. FRC-32
#gem	Geobacter sp. M21
#geb	Geobacter sp. M18
#pca	Pelobacter carbinolicus
#ppd	Pelobacter propionicus
#dvu	Desulfovibrio vulgaris Hildenborough
#dvl	Desulfovibrio vulgaris DP4
#dvm	Desulfovibrio vulgaris Miyazaki F
#dde	Desulfovibrio desulfuricans G20
#dds	Desulfovibrio desulfuricans ATCC 27774
#dma	Desulfovibrio magneticus
#dsa	Desulfovibrio salexigens
#das	Desulfovibrio aespoeensis
#lip	Lawsonia intracellularis
#dba	Desulfomicrobium baculatum
#drt	Desulfohalobium retbaense
#bba	Bdellovibrio bacteriovorus
#dps	Desulfotalea psychrophila
#dak	Desulfurivibrio alkaliphilus
#dpr	Desulfobulbus propionicus
#dol	Candidatus Desulfococcus oleovorans
#dal	Desulfatibacillum alkenivorans
#dat	Desulfobacterium autotrophicum
#ade	Anaeromyxobacter dehalogenans 2CP-C
#acp	Anaeromyxobacter dehalogenans 2CP-1
#afw	Anaeromyxobacter sp. Fw109-5
#ank	Anaeromyxobacter sp. K
#mxa	Myxococcus xanthus
#mfu	Myxococcus fulvus
#sur	Stigmatella aurantiaca
#scl	Sorangium cellulosum
#hoh	Haliangium ochraceum
#sat	Syntrophus aciditrophicus
#dao	Desulfobacca acetoxidans DSM 11109
#sfu	Syntrophobacter fumaroxidans
#dbr	Desulfarculus baarsii
#hmr	Hippea maritima
#rpr	Rickettsia prowazekii
#rty	Rickettsia typhi
#rcm	Rickettsia canadensis
#rco	Rickettsia conorii
#rfe	Rickettsia felis
#rak	Rickettsia akari
#rri	Rickettsia rickettsii Sheila Smith
#rrj	Rickettsia rickettsii Iowa
#rms	Rickettsia massiliae
#rpk	Rickettsia peacockii
#raf	Rickettsia africae
#rhe	Rickettsia heilongjiangensis
#rbe	Rickettsia bellii RML369-C
#rbo	Rickettsia bellii OSU 85-389
#ots	Orientia tsutsugamushi Boryong
#ott	Orientia tsutsugamushi Ikeda
#wol	Wolbachia wMel
#wbm	Wolbachia wBm
#wri	Wolbachia sp. wRi
#wpi	Wolbachia pipientis
#ama	Anaplasma marginale St. Maries
#amf	Anaplasma marginale Florida
#acn	Anaplasma centrale
#aph	Anaplasma phagocytophilum
#eru	Ehrlichia ruminantium Welgevonden (South Africa)
#erw	Ehrlichia ruminantium Welgevonden (France)
#erg	Ehrlichia ruminantium Gardel
#ecn	Ehrlichia canis
#ech	Ehrlichia chaffeensis
#nse	Neorickettsia sennetsu
#nri	Neorickettsia risticii
#pub	Candidatus Pelagibacter ubique
#pel	Candidatus Pelagibacter sp. IMCC9063
#mmn	Candidatus Midichloria mitochondrii
#mlo	Mesorhizobium loti
#mci	Mesorhizobium ciceri
#mop	Mesorhizobium opportunistum
#mes	Mesorhizobium sp. BNC1
#pla	Parvibaculum lavamentivorans
#sme	Sinorhizobium meliloti 1021
#smk	Sinorhizobium meliloti AK83
#smd	Sinorhizobium medicae
#rhi	Sinorhizobium fredii NGR234
#atu	Agrobacterium tumefaciens C58
#ara	Agrobacterium radiobacter K84
#avi	Agrobacterium vitis S4
#agr	Agrobacterium sp. H13-3
#ret	Rhizobium etli CFN 42
#rec	Rhizobium etli CIAT 652
#rle	Rhizobium leguminosarum
#rlt	Rhizobium leguminosarum bv. trifolii WSM2304
#rlg	Rhizobium leguminosarum bv. trifolii WSM1325
#las	Candidatus Liberibacter asiaticus
#lso	Candidatus Liberibacter solanacearum
#bme	Brucella melitensis bv. 1 16M
#bmi	Brucella melitensis ATCC 23457
#bmf	Brucella melitensis biovar Abortus
#bmb	Brucella abortus 9-941
#bmc	Brucella abortus S19
#bms	Brucella suis 1330
#bmt	Brucella suis ATCC 23445
#bov	Brucella ovis
#bcs	Brucella canis
#bmr	Brucella microti
#bpp	Brucella pinnipedialis
#oan	Ochrobactrum anthropi
#bja	Bradyrhizobium japonicum
#bra	Bradyrhizobium sp. ORS278
#bbt	Bradyrhizobium sp. BTAi1
#rpa	Rhodopseudomonas palustris CGA009
#rpb	Rhodopseudomonas palustris HaA2
#rpc	Rhodopseudomonas palustris BisB18
#rpd	Rhodopseudomonas palustris BisB5
#rpe	Rhodopseudomonas palustris BisA53
#rpt	Rhodopseudomonas palustris TIE-1
#rpx	Rhodopseudomonas palustris DX-1
#nwi	Nitrobacter winogradskyi
#nha	Nitrobacter hamburgensis
#oca	Oligotropha carboxidovorans
#bhe	Bartonella henselae
#bqu	Bartonella quintana
#bbk	Bartonella bacilliformis
#btr	Bartonella tribocorum
#bgr	Bartonella grahamii
#bcd	Bartonella clarridgeiae
#xau	Xanthobacter autotrophicus
#azc	Azorhizobium caulinodans
#sno	Starkeya novella
#mex	Methylobacterium extorquens
#mea	Methylobacterium extorquens AM1
#mdi	Methylobacterium extorquens DM4
#mrd	Methylobacterium radiotolerans
#met	Methylobacterium sp. 4-46
#mpo	Methylobacterium populi
#mch	Methylobacterium chloromethanicum
#mno	Methylobacterium nodulans
#bid	Beijerinckia indica
#msl	Methylocella silvestris
#hdn	Hyphomicrobium denitrificans
#hmc	Hyphomicrobium sp. MC1
#rva	Rhodomicrobium vannielii
#hci	Candidatus Hodgkinia cicadicola
#ccr	Caulobacter crescentus CB15
#ccs	Caulobacter crescentus NA1000
#cak	Caulobacter sp. K31
#cse	Caulobacter segnis
#pzu	Phenylobacterium zucineum
#bsb	Brevundimonas subvibrioides
#aex	Asticcacaulis excentricus
#sil	Silicibacter pomeroyi
#sit	Ruegeria sp. TM1040
#rsp	Rhodobacter sphaeroides 2.4.1
#rsh	Rhodobacter sphaeroides ATCC 17029
#rsq	Rhodobacter sphaeroides ATCC 17025
#rsk	Rhodobacter sphaeroides KD131
#rcp	Rhodobacter capsulatus
#jan	Jannaschia sp. CCS1
#rde	Roseobacter denitrificans
#rli	Roseobacter litoralis Och 149
#pde	Paracoccus denitrificans
#dsh	Dinoroseobacter shibae
#kvu	Ketogulonicigenium vulgare
#mmr	Maricaulis maris
#hne	Hyphomonas neptunium
#hba	Hirschia baltica
#zmo	Zymomonas mobilis
#zmn	Zymomonas mobilis subsp. mobilis NCIMB 11163
#zmp	Zymomonas mobilis subsp. pomaceae ATCC 29192
#nar	Novosphingobium aromaticivorans
#npp	Novosphingobium sp. PP1Y
#sal	Sphingopyxis alaskensis
#swi	Sphingomonas wittichii
#sjp	Sphingobium japonicum
#sch	Sphingobium chlorophenolicum
#eli	Erythrobacter litoralis
#gox	Gluconobacter oxydans
#gbe	Granulibacter bethesdensis
#acr	Acidiphilium cryptum JF-5
#amv	Acidiphilium multivorum
#gdi	Gluconacetobacter diazotrophicus PAl 5 (Brazil)
#gdj	Gluconacetobacter diazotrophicus PAl 5 (JGI)
#apt	Acetobacter pasteurianus
#rru	Rhodospirillum rubrum
#rce	Rhodospirillum centenum
#mag	Magnetospirillum magneticum
#azl	Azospirillum sp. B510
#pbr	Parvularcula bermudensis
#apb	Candidatus Puniceispirillum marinum
#pgv	Polymorphum gilvum
#mgm	Magnetococcus sp. MC-1
#din	Desulfurispirillum indicum
#bsu	Bacillus subtilis
#bss	Bacillus subtilis subsp. spizizenii
#bsn	Bacillus subtilis BSn5
#bha	Bacillus halodurans
#ban	Bacillus anthracis Ames
#bar	Bacillus anthracis Ames 0581
#bat	Bacillus anthracis Sterne
#bah	Bacillus anthracis CDC 684
#bai	Bacillus anthracis A0248
#bal	Bacillus cereus biovar anthracis CI
#bce	Bacillus cereus ATCC 14579
#bca	Bacillus cereus ATCC 10987
#bcz	Bacillus cereus ZK
#bcr	Bacillus cereus AH187
#bcb	Bacillus cereus B4264
#bcu	Bacillus cereus AH820
#bcg	Bacillus cereus G9842
#bcq	Bacillus cereus Q1
#bcx	Bacillus cereus 03BB102
#bcy	Bacillus cytotoxis NVH 391-98
#btk	Bacillus thuringiensis 97-27
#btl	Bacillus thuringiensis Al Hakam
#btb	Bacillus thuringiensis BMB171
#bwe	Bacillus weihenstephanensis
#bli	Bacillus licheniformis ATCC 14580
#bld	Bacillus licheniformis DSM13
#bay	Bacillus amyloliquefaciens FZB42
#bao	Bacillus amyloliquefaciens DSM 7
#bae	Bacillus atrophaeus
#bcl	Bacillus clausii
#bpu	Bacillus pumilus
#bpf	Bacillus pseudofirmus
#bmq	Bacillus megaterium QM B1551
#bmd	Bacillus megaterium DSM 319
#bse	Bacillus selenitireducens
#bco	Bacillus cellulosilyticus
#bck	Bacillus coagulans
#oih	Oceanobacillus iheyensis
#gka	Geobacillus kaustophilus
#gtn	Geobacillus thermodenitrificans
#gth	Geobacillus thermoglucosidasius
#gwc	Geobacillus sp. WCH70
#gyc	Geobacillus sp. Y412MC61
#gya	Geobacillus sp. Y412MC52
#gct	Geobacillus sp. C56-T3
#gmc	Geobacillus sp. Y4.1MC1
#afl	Anoxybacillus flavithermus
#sau	Staphylococcus aureus N315 (MRSA/VSSA)
#sav	Staphylococcus aureus Mu50 (MRSA/VISA)
#saw	Staphylococcus aureus Mu3 (MRSA/hetero-VISA)
#sah	Staphylococcus aureus JH1 (MRSA/VSSA)
#saj	Staphylococcus aureus JH9 (MRSA/VRSA)
#sam	Staphylococcus aureus MW2 (CA-MRSA)
#sas	Staphylococcus aureus MSSA476 (MSSA)
#sar	Staphylococcus aureus MRSA252 (MRSA)
#sac	Staphylococcus aureus COL (MRSA)
#sax	Staphylococcus aureus USA300_TCH1516 (CA-MSSA)
#saa	Staphylococcus aureus USA300_FPR3757 (CA-MRSA)
#sao	Staphylococcus aureus NCTC8325
#sae	Staphylococcus aureus Newman
#sad	Staphylococcus aureus ED98
#sab	Staphylococcus aureus RF122
#sep	Staphylococcus epidermidis ATCC 12228
#ser	Staphylococcus epidermidis RP62A
#sha	Staphylococcus haemolyticus
#ssp	Staphylococcus saprophyticus
#sca	Staphylococcus carnosus
#slg	Staphylococcus lugdunensis
#ssd	Staphylococcus pseudintermedius
#lmo	Listeria monocytogenes EGD-e
#lmf	Listeria monocytogenes F2365
#lmh	Listeria monocytogenes HCC23
#lmc	Listeria monocytogenes Clip81459
#lmn	Listeria monocytogenes 08-5578
#lmy	Listeria monocytogenes 08-5923
#lin	Listeria innocua
#lwe	Listeria welshimeri SLCC5334
#lsg	Listeria seeligeri
#lsp	Lysinibacillus sphaericus
#esi	Exiguobacterium sibiricum
#eat	Exiguobacterium sp. AT1b
#mcl	Macrococcus caseolyticus
#bbe	Brevibacillus brevis
#pjd	Paenibacillus sp. JDR-2
#gym	Geobacillus sp. Y412MC10
#ppy	Paenibacillus polymyxa E681
#ppm	Paenibacillus polymyxa SC2
#pms	Paenibacillus mucilaginosus
#aac	Alicyclobacillus acidocaldarius
#bts	Bacillus tusciae
#lla	Lactococcus lactis subsp. lactis IL1403
#llk	Lactococcus lactis subsp. lactis KF147
#llc	Lactococcus lactis subsp. cremoris SK11
#llm	Lactococcus lactis subsp. cremoris MG1363
#spy	Streptococcus pyogenes SF370 (serotype M1)
#spz	Streptococcus pyogenes MGAS5005 (serotype M1)
#spm	Streptococcus pyogenes MGAS8232 (serotype M18)
#spg	Streptococcus pyogenes MGAS315 (serotype M3)
#sps	Streptococcus pyogenes SSI-1 (serotype M3)
#sph	Streptococcus pyogenes MGAS10270 (serotype M2)
#spi	Streptococcus pyogenes MGAS10750 (serotype M4)
#spj	Streptococcus pyogenes MGAS2096 (serotype M12)
#spk	Streptococcus pyogenes MGAS9429 (serotype M12)
#spf	Streptococcus pyogenes Manfredo (serotype M5)
#spa	Streptococcus pyogenes MGAS10394 (serotype M6)
#spb	Streptococcus pyogenes MGAS6180 (serotype M28)
#soz	Streptococcus pyogenes NZ131 (serotype M49)
#spn	Streptococcus pneumoniae TIGR4 (virulent serotype 4)
#spd	Streptococcus pneumoniae D39 (virulent serotype 2)
#spr	Streptococcus pneumoniae R6 (avirulent)
#spw	Streptococcus pneumoniae CGSP14 (serotype 14)
#spx	Streptococcus pneumoniae G54 (serotype 19F)
#sne	Streptococcus pneumoniae ATCC 700669 (serotype 23F ST81 lineage)
#spv	Streptococcus pneumoniae Hungary19A 6
#snm	Streptococcus pneumoniae 70585
#sjj	Streptococcus pneumoniae JJA
#spp	Streptococcus pneumoniae P1031
#snt	Streptococcus pneumoniae Taiwan19F-14
#snc	Streptococcus pneumoniae TCH8431/19A
#snb	Streptococcus pneumoniae 670-6B
#snp	Streptococcus pneumoniae AP200
#sag	Streptococcus agalactiae 2603 (serotype V)
#san	Streptococcus agalactiae NEM316 (serotype III)
#sak	Streptococcus agalactiae A909 (serotype Ia)
#smu	Streptococcus mutans UA159
#smc	Streptococcus mutans NN2025
#stc	Streptococcus thermophilus CNRZ1066
#stl	Streptococcus thermophilus LMG18311
#ste	Streptococcus thermophilus LMD-9
#ssa	Streptococcus sanguinis
#ssu	Streptococcus suis 05ZYH33
#ssv	Streptococcus suis 98HAH33
#ssb	Streptococcus suis BM407
#ssi	Streptococcus suis P1/7
#sss	Streptococcus suis SC84
#sst	Streptococcus suis ST3
#sgo	Streptococcus gordonii
#seq	Streptococcus equi subsp. zooepidemicus H70
#sez	Streptococcus equi subsp. zooepidemicus MGCS10565
#seu	Streptococcus equi subsp. equi 4047
#sub	Streptococcus uberis
#sds	Streptococcus dysgalactiae
#sga	Streptococcus gallolyticus UCN34
#sgg	Streptococcus gallolyticus subsp. gallolyticus
#smb	Streptococcus mitis B6
#sor	Streptococcus oralis
#stk	Streptococcus parauberis
#stb	Streptococcus pasteurianus
#scp	Streptococcus parasanguinis
#ssr	Streptococcus salivarius
#std	Streptococcus pseudopneumoniae
#lpl	Lactobacillus plantarum WCFS1
#lpj	Lactobacillus plantarum JDM1
#lps	Lactobacillus plantarum subsp. plantarum ST-III
#ljo	Lactobacillus johnsonii NCC 533
#ljf	Lactobacillus johnsonii FI9785
#lac	Lactobacillus acidophilus NCFM
#lai	Lactobacillus acidophilus 30SC
#lsa	Lactobacillus sakei
#lsl	Lactobacillus salivarius
#ldb	Lactobacillus delbrueckii ATCC 11842
#lbu	Lactobacillus delbrueckii ATCC BAA-365
#lde	Lactobacillus delbrueckii subsp. bulgaricus ND02
#lbr	Lactobacillus brevis
#lca	Lactobacillus casei ATCC 334
#lcb	Lactobacillus casei BL23
#lcz	Lactobacillus casei Zhang
#lga	Lactobacillus gasseri
#lre	Lactobacillus reuteri DSM 20016
#lrf	Lactobacillus reuteri JCM 1112
#lru	Lactobacillus reuteri SD2112
#lhe	Lactobacillus helveticus
#lfe	Lactobacillus fermentum
#lrh	Lactobacillus rhamnosus GG
#lrl	Lactobacillus rhamnosus Lc 705
#lcr	Lactobacillus crispatus
#lam	Lactobacillus amylovorus
#lbh	Lactobacillus buchneri
#lke	Lactobacillus kefiranofaciens
#ppe	Pediococcus pentosaceus
#efa	Enterococcus faecalis
#mps	Melissococcus plutonius
#ooe	Oenococcus oeni
#lme	Leuconostoc mesenteroides
#lci	Leuconostoc citreum
#lki	Leuconostoc kimchii
#lgs	Leuconostoc gasicomitatum
#lec	Leuconostoc sp. C2
#aur	Aerococcus urinae
#crn	Carnobacterium sp. 17-4
#wko	Weissella koreensis
#cac	Clostridium acetobutylicum ATCC 824
#cae	Clostridium acetobutylicum DSM 1731
#cpe	Clostridium perfringens 13
#cpf	Clostridium perfringens ATCC 13124
#cpr	Clostridium perfringens SM101
#ctc	Clostridium tetani E88
#cno	Clostridium novyi
#cth	Clostridium thermocellum
#cdf	Clostridium difficile 630
#cdc	Clostridium difficile CD196
#cdl	Clostridium difficile R20291
#cbo	Clostridium botulinum A ATCC 3502
#cba	Clostridium botulinum A ATCC 19397
#cbh	Clostridium botulinum A Hall
#cby	Clostridium botulinum A2
#cbl	Clostridium botulinum A3 Loch Maree
#cbk	Clostridium botulinum B Eklund 17B
#cbb	Clostridium botulinum B1 Okra
#cbi	Clostridium botulinum Ba4
#cbn	Clostridium botulinum BKT015925
#cbt	Clostridium botulinum E3
#cbf	Clostridium botulinum F Langeland
#cbe	Clostridium beijerinckii
#ckl	Clostridium kluyveri DSM 555
#ckr	Clostridium kluyveri NBRC 12016
#cpy	Clostridium phytofermentans
#cce	Clostridium cellulolyticum
#clj	Clostridium ljungdahlii
#csh	Clostridium saccharolyticum
#ccb	Clostridium cellulovorans
#cst	Clostridium sticklandii
#cls	Clostridium sp. SY8519
#amt	Alkaliphilus metalliredigens
#aoe	Alkaliphilus oremlandii
#asf	Candidatus Arthromitus sp. SFB-mouse-Japan
#sth	Symbiobacterium thermophilum
#swo	Syntrophomonas wolfei
#slp	Syntrophothermus lipocalidus
#vpr	Veillonella parvula
#ssg	Selenomonas sputigena
#afn	Acidaminococcus fermentans
#dsy	Desulfitobacterium hafniense Y51
#dhd	Desulfitobacterium hafniense DCB-2
#drm	Desulfotomaculum reducens
#dae	Desulfotomaculum acetoxidans
#dca	Desulfotomaculum carboxydivorans
#dku	Desulfotomaculum kuznetsovii
#dru	Desulfotomaculum ruminis
#pth	Pelotomaculum thermopropionicum
#dau	Candidatus Desulforudis audaxviator
#tjr	Thermincola potens JR
#sgy	Syntrophobotulus glycolicus
#hmo	Heliobacterium modesticaldum
#fma	Finegoldia magna
#apr	Anaerococcus prevotii
#eel	Eubacterium eligens
#ere	Eubacterium rectale
#elm	Eubacterium limosum
#bpb	Butyrivibrio proteoclasticus
#cle	Clostridium lentocellum
#eha	Ethanoligenens harbinense
#ral	Ruminococcus albus
#tmr	Thermaerobacter marianensis
#say	Sulfobacillus acidophilus
#clo	Clostridiales genomosp. BVAB3
#tte	Thermoanaerobacter tengcongensis
#tex	Thermoanaerobacter sp. X514
#thx	Thermoanaerobacter sp. X513
#tpd	Thermoanaerobacter pseudethanolicus
#tit	Thermoanaerobacter italicus
#tmt	Thermoanaerobacter mathranii
#tbo	Thermoanaerobacter brockii
#chy	Carboxydothermus hydrogenoformans
#tep	Tepidanaerobacter sp. Re1
#mta	Moorella thermoacetica
#adg	Ammonifex degensii
#csc	Caldicellulosiruptor saccharolyticus
#ate	Caldicellulosiruptor bescii
#cob	Caldicellulosiruptor obsidiansis
#chd	Caldicellulosiruptor hydrothermalis
#cow	Caldicellulosiruptor owensensis
#cki	Caldicellulosiruptor kristjanssonii
#ckn	Caldicellulosiruptor kronotskyensis
#toc	Thermosediminibacter oceani
#ttm	Thermoanaerobacterium thermosaccharolyticum
#txy	Thermoanaerobacterium xylanolyticum
#cpo	Coprothermobacter proteolyticus
#tnr	Thermodesulfobium narugense
#mas	Mahella australiensis
#nth	Natranaerobius thermophilus
#hor	Halothermothrix orenii
#has	Halanaerobium hydrogeniformans
#aar	Acetohalobium arabaticum
#erh	Erysipelothrix rhusiopathiae
#mge	Mycoplasma genitalium
#mpn	Mycoplasma pneumoniae
#mpu	Mycoplasma pulmonis
#mpe	Mycoplasma penetrans
#mga	Mycoplasma gallisepticum
#mmy	Mycoplasma mycoides subsp. mycoides SC PG1
#mml	Mycoplasma mycoides subsp. capri LC 95010
#mmo	Mycoplasma mobile
#mhy	Mycoplasma hyopneumoniae 232
#mhj	Mycoplasma hyopneumoniae J
#mhp	Mycoplasma hyopneumoniae 7448
#msy	Mycoplasma synoviae
#mcp	Mycoplasma capricolum
#maa	Mycoplasma agalactiae PG2
#mal	Mycoplasma agalactiae 5632
#mat	Mycoplasma arthritidis
#mco	Mycoplasma conjunctivae
#mho	Mycoplasma hominis
#mcd	Mycoplasma crocodyli
#mhr	Mycoplasma hyorhinis
#mfr	Mycoplasma fermentans JER
#mfm	Mycoplasma fermentans M64
#mbv	Mycoplasma bovis PG45
#mbh	Mycoplasma bovis Hubei-1
#mlc	Mycoplasma leachii
#mha	Mycoplasma haemofelis
#mss	Mycoplasma suis Illinois
#msk	Mycoplasma suis KI3806
#mpf	Mycoplasma putrefaciens
#uur	Ureaplasma parvum serovar 3 ATCC 700970
#upa	Ureaplasma parvum serovar 3 ATCC 27815
#uue	Ureaplasma urealyticum serovar 10 ATCC 33699
#poy	Phytoplasma OY
#ayw	Phytoplasma AYWB
#pml	Candidatus Phytoplasma mali
#pal	Candidatus Phytoplasma australiense
#acl	Acholeplasma laidlawii
#mfl	Mesoplasma florum
#mtu	Mycobacterium tuberculosis H37Rv
#mtc	Mycobacterium tuberculosis CDC1551
#mra	Mycobacterium tuberculosis H37Ra
#mtf	Mycobacterium tuberculosis F11
#mtb	Mycobacterium tuberculosis KZN 1435
#mbo	Mycobacterium bovis AF2122/97
#mbb	Mycobacterium bovis BCG Pasteur 1173P2
#mbt	Mycobacterium bovis BCG Tokyo 172
#maf	Mycobacterium africanum
#mce	Mycobacterium canettii
#mle	Mycobacterium leprae TN
#mlb	Mycobacterium leprae Br4923
#mpa	Mycobacterium avium paratuberculosis
#mav	Mycobacterium avium 104
#msm	Mycobacterium smegmatis
#mul	Mycobacterium ulcerans
#mva	Mycobacterium vanbaalenii
#mgi	Mycobacterium gilvum
#mab	Mycobacterium abscessus ATCC 19977
#mmc	Mycobacterium sp. MCS
#mkm	Mycobacterium sp. KMS
#mjl	Mycobacterium sp. JLS
#msp	Mycobacterium sp. Spyr1
#mjd	Mycobacterium sp. JDM601
#mmi	Mycobacterium marinum M
#asd	Amycolicicoccus subflavus
#cgl	Corynebacterium glutamicum ATCC 13032 (Kyowa Hakko)
#cgb	Corynebacterium glutamicum ATCC 13032 (Bielefeld)
#cgt	Corynebacterium glutamicum R
#cef	Corynebacterium efficiens
#cdi	Corynebacterium diphtheriae
#cjk	Corynebacterium jeikeium
#cur	Corynebacterium urealyticum
#car	Corynebacterium aurimucosum
#ckp	Corynebacterium kroppenstedtii
#cpu	Corynebacterium pseudotuberculosis
#crd	Corynebacterium resistens
#cul	Corynebacterium ulcerans
#cva	Corynebacterium variabile
#nfa	Nocardia farcinica
#rha	Rhodococcus sp. RHA1
#rer	Rhodococcus erythropolis
#rop	Rhodococcus opacus
#req	Rhodococcus equi
#gbr	Gordonia bronchialis
#tpr	Tsukamurella paurometabola
#srt	Segniliparus rotundus
#sco	Streptomyces coelicolor
#sma	Streptomyces avermitilis
#sgr	Streptomyces griseus
#scb	Streptomyces scabiei
#twh	Tropheryma whipplei Twist
#tws	Tropheryma whipplei TW08/27
#lxx	Leifsonia xyli xyli CTCB07
#cmi	Clavibacter michiganensis subsp. michiganensis
#cms	Clavibacter michiganensis subsp. sepedonicus
#mts	Microbacterium testaceum
#art	Arthrobacter sp. FB24
#aau	Arthrobacter aurescens
#ach	Arthrobacter chlorophenolicus
#aai	Arthrobacter arilaitensis
#apn	Arthrobacter phenanthrenivorans
#rsa	Renibacterium salmoninarum
#krh	Kocuria rhizophila
#mlu	Micrococcus luteus
#rmu	Rothia mucilaginosa
#rdn	Rothia dentocariosa
#bcv	Beutenbergia cavernae
#bfa	Brachybacterium faecium
#jde	Jonesia denitrificans
#kse	Kytococcus sedentarius
#xce	Xylanimonas cellulosilytica
#iva	Isoptericola variabilis
#ske	Sanguibacter keddieii
#cfl	Cellulomonas flavigena
#cfi	Cellulomonas fimi
#ica	Intrasporangium calvum
#pac	Propionibacterium acnes KPA171202
#pak	Propionibacterium acnes SK137
#pfr	Propionibacterium freudenreichii
#mph	Microlunatus phosphovorus
#nca	Nocardioides sp. JS614
#kfl	Kribbella flavida
#tfu	Thermobifida fusca
#nda	Nocardiopsis dassonvillei
#tcu	Thermomonospora curvata
#sro	Streptosporangium roseum
#fra	Frankia sp. CcI3
#fre	Frankia sp. EAN1pec
#fri	Frankia sp. EuI1c
#fal	Frankia alni
#fsy	Frankia symbiont
#ace	Acidothermus cellulolyticus
#nml	Nakamurella multipartita
#gob	Geodermatophilus obscurus
#kra	Kineococcus radiotolerans
#sen	Saccharopolyspora erythraea
#svi	Saccharomonospora viridis
#tbi	Thermobispora bispora
#amd	Amycolatopsis mediterranei
#pdx	Pseudonocardia dioxanivorans
#ami	Actinosynnema mirum
#stp	Salinispora tropica
#saq	Salinispora arenicola
#mau	Micromonospora aurantiaca
#mil	Micromonospora sp. L5
#vma	Verrucosispora maris
#cai	Catenulispora acidiphila
#sna	Stackebrandtia nassauensis
#ahe	Arcanobacterium haemolyticum
#mcu	Mobiluncus curtisii
#blo	Bifidobacterium longum NCC2705
#blj	Bifidobacterium longum DJO10A
#bln	Bifidobacterium longum subsp. infantis ATCC 15697
#blf	Bifidobacterium longum subsp. infantis 157F
#bll	Bifidobacterium longum subsp. longum JDM301
#blb	Bifidobacterium longum subsp. longum BBMN68
#blm	Bifidobacterium longum subsp. longum JCM 1217
#bad	Bifidobacterium adolescentis
#bla	Bifidobacterium animalis subsp. lactis AD011
#blc	Bifidobacterium animalis subsp. lactis Bl-04
#blt	Bifidobacterium animalis subsp. lactis DSM 10140
#bde	Bifidobacterium dentium
#bbi	Bifidobacterium bifidum S17
#bbp	Bifidobacterium bifidum PRL2010
#gva	Gardnerella vaginalis
#gvg	Gardnerella vaginalis ATCC 14019
#rxy	Rubrobacter xylanophilus
#cwo	Conexibacter woesei
#afo	Acidimicrobium ferrooxidans DSM 10331
#ccu	Cryptobacterium curtum
#shi	Slackia heliotrinireducens
#apv	Atopobium parvulum
#ele	Eggerthella lenta
#eyy	Eggerthella sp. YY7918
#ols	Olsenella uli
#cgo	Coriobacterium glomerans
#ctr	Chlamydia trachomatis D/UW-3/CX
#cta	Chlamydia trachomatis A/HAR-13
#ctb	Chlamydia trachomatis L2/434/Bu
#ctl	Chlamydia trachomatis L2b/UCH-1/proctitis
#cto	Chlamydia trachomatis L2c
#ctj	Chlamydia trachomatis B/Jali20/OT
#ctz	Chlamydia trachomatis B/TZ1A828/OT
#cmu	Chlamydia muridarum
#cpn	Chlamydophila pneumoniae CWL029
#cpa	Chlamydophila pneumoniae AR39
#cpj	Chlamydophila pneumoniae J138
#cpt	Chlamydophila pneumoniae TW183
#cca	Chlamydophila caviae
#cab	Chlamydophila abortus
#cfe	Chlamydophila felis
#cpm	Chlamydophila pecorum
#chp	Chlamydophila psittaci
#pcu	Candidatus Protochlamydia amoebophila
#puv	Parachlamydia acanthamoebae
#wch	Waddlia chondrophila
#sng	Simkania negevensis
#bbu	Borrelia burgdorferi B31
#bbz	Borrelia burgdorferi ZS7
#bga	Borrelia garinii
#baf	Borrelia afzelii
#bbs	Borrelia bissettii
#btu	Borrelia turicatae
#bhr	Borrelia hermsii
#bdu	Borrelia duttonii
#bre	Borrelia recurrentis
#tpa	Treponema pallidum subsp. pallidum Nichols
#tpp	Treponema pallidum subsp. pallidum SS14
#tde	Treponema denticola
#tsu	Treponema succinifaciens
#tbe	Treponema brennaborense
#taz	Treponema azotonutricium
#tpi	Treponema primitia
#tpl	Treponema paraluiscuniculi
#ssm	Spirochaeta smaragdinae
#sta	Spirochaeta thermophila
#sbu	Spirochaeta sp. Buddy
#scc	Spirochaeta coccoides
#scd	Spirochaeta caldaria DSM 7334
#lil	Leptospira interrogans serovar lai
#lic	Leptospira interrogans serovar Copenhageni
#lbj	Leptospira borgpetersenii JB197
#lbl	Leptospira borgpetersenii L550
#lbi	Leptospira biflexa serovar Patoc Patoc 1 (Paris)
#lbf	Leptospira biflexa serovar Patoc Patoc 1 (Ames)
#bhy	Brachyspira hyodysenteriae
#brm	Brachyspira murdochii
#bpo	Brachyspira pilosicoli
#aba	Candidatus Koribacter versatilis
#aca	Acidobacterium capsulatum
#acm	Acidobacterium sp. MP5ACTX9
#tsa	Terriglobus saanensis
#sus	Candidatus Solibacter usitatus
#bth	Bacteroides thetaiotaomicron
#bfr	Bacteroides fragilis YCH46
#bfs	Bacteroides fragilis NCTC9343
#bvu	Bacteroides vulgatus
#bhl	Bacteroides helcogenes
#bsa	Bacteroides salanitronis
#pgi	Porphyromonas gingivalis W83
#pgn	Porphyromonas gingivalis ATCC 33277
#pgt	Porphyromonas gingivalis TDC60
#pah	Porphyromonas asaccharolytica
#pdi	Parabacteroides distasonis
#ppn	Paludibacter propionicigenes
#osp	Odoribacter splanchnicus
#aps	Candidatus Azobacteroides pseudotrichonymphae
#pru	Prevotella ruminicola
#pmz	Prevotella melaninogenica
#pdn	Prevotella denticola
#sru	Salinibacter ruber
#srm	Salinibacter ruber
#rmr	Rhodothermus marinus
#cpi	Chitinophaga pinensis
#phe	Pedobacter heparinus
#psn	Pedobacter saltans
#shg	Sphingobacterium sp. 21
#hhy	Haliscomenobacter hydrossis
#cmr	Cyclobacterium marinum
#chu	Cytophaga hutchinsonii
#dfe	Dyadobacter fermentans
#sli	Spirosoma linguale
#lby	Leadbetterella byssophila
#rsi	Runella slithyformis
#mtt	Marivirga tractuosa
#gfo	Gramella forsetii
#fjo	Flavobacterium johnsoniae
#fps	Flavobacterium psychrophilum
#coc	Capnocytophaga ochracea
#ccm	Capnocytophaga canimorsus
#rbi	Robiginitalea biformata
#zpr	Zunongwangia profunda
#cat	Croceibacter atlanticus
#ran	Riemerella anatipestifer
#fbc	Maribacter sp. HTCC2170
#cao	Cellulophaga algicola
#cly	Cellulophaga lytica
#wvi	Weeksella virosa
#kdi	Krokinobacter sp. 4H-3-7-5
#lan	Lacinutrix sp. 5H-3-7-4
#zga	Zobellia galactanivorans
#mrs	Muricauda ruestringensis
#fba	Flavobacteriaceae bacterium
#smg	Candidatus Sulcia muelleri GWSS
#sms	Candidatus Sulcia muelleri SMDSEM
#smh	Candidatus Sulcia muelleri DMIN
#sum	Candidatus Sulcia muelleri CARI
#bbl	Blattabacterium sp. (Blattella germanica)
#bpi	Blattabacterium sp. (Periplaneta americana)
#fte	Fluviicola taffensis
#aas	Candidatus Amoebophilus asiaticus
#fsu	Fibrobacter succinogenes
#fnu	Fusobacterium nucleatum
#lba	Leptotrichia buccalis
#str	Sebaldella termitidis
#smf	Streptobacillus moniliformis
#ipo	Ilyobacter polytropus
#ote	Opitutus terrae
#caa	Coraliomargarita akajimensis
#min	Methylacidiphilum infernorum
#amu	Akkermansia muciniphila
#gau	Gemmatimonas aurantiaca
#rba	Rhodopirellula baltica
#psl	Pirellula staleyi
#plm	Planctomyces limnophilus
#pbs	Planctomyces brasiliensis
#ipa	Isosphaera pallida
#emi	Elusimicrobium minutum
#rsd	Uncultured Termite group 1 bacterium phylotype Rs-D17
#tai	Thermanaerovibrio acidaminovorans
#aco	Aminobacterium colombiense
#syn	Synechocystis sp. PCC6803
#syw	Synechococcus sp. WH8102
#syc	Synechococcus elongatus PCC6301
#syf	Synechococcus elongatus PCC7942
#syd	Synechococcus sp. CC9605
#sye	Synechococcus sp. CC9902
#syg	Synechococcus sp. CC9311
#syr	Synechococcus sp. RCC307
#syx	Synechococcus sp. WH7803
#syp	Synechococcus sp. PCC7002
#cya	Synechococcus sp. JA-3-3Ab
#cyb	Synechococcus sp. JA-2-3B'a(2-13)
#tel	Thermosynechococcus elongatus
#mar	Microcystis aeruginosa
#cyt	Cyanothece sp. ATCC 51142
#cyp	Cyanothece sp. PCC 8801
#cyc	Cyanothece sp. PCC 7424
#cyn	Cyanothece sp. PCC 7425
#cyh	Cyanothece sp. PCC 8802
#cyj	Cyanothece sp. PCC 7822
#cyu	Cyanobacterium UCYN-A
#gvi	Gloeobacter violaceus
#ana	Anabaena sp. PCC7120
#npu	Nostoc punctiforme
#ava	Anabaena variabilis
#naz	Anabaena azollae 0708
#pma	Prochlorococcus marinus SS120
#pmm	Prochlorococcus marinus MED4
#pmt	Prochlorococcus marinus MIT 9313
#pmn	Prochlorococcus marinus NATL2A
#pmi	Prochlorococcus marinus MIT9312
#pmb	Prochlorococcus marinus AS9601
#pmc	Prochlorococcus marinus MIT 9515
#pmf	Prochlorococcus marinus MIT 9303
#pmg	Prochlorococcus marinus MIT 9301
#pmh	Prochlorococcus marinus MIT 9215
#pmj	Prochlorococcus marinus MIT 9211
#pme	Prochlorococcus marinus NATL1A
#ter	Trichodesmium erythraeum
#amr	Acaryochloris marina
#cte	Chlorobaculum tepidum
#cpc	Chlorobaculum parvum NCIB 8327
#cch	Chlorobium chlorochromatii
#cph	Chlorobium phaeobacteroides DSM 266
#cpb	Chlorobium phaeobacteroides BS1
#cli	Chlorobium limicola
#pvi	Chlorobium phaeovibrioides
#plt	Pelodictyon luteolum
#pph	Pelodictyon phaeoclathratiforme
#paa	Prosthecochloris aestuarii
#cts	Chloroherpeton thalassium
#det	Dehalococcoides ethenogenes
#deh	Dehalococcoides sp. CBDB1
#deb	Dehalococcoides sp. BAV1
#dev	Dehalococcoides sp. VS
#deg	Dehalococcoides sp. GT
#dly	Dehalogenimonas lykanthroporepellens
#rrs	Roseiflexus sp. RS-1
#rca	Roseiflexus castenholzii DSM13941
#cau	Chloroflexus aurantiacus
#cag	Chloroflexus aggregans
#chl	Chloroflexus sp. Y-400-fl
#hau	Herpetosiphon aurantiacus
#tro	Thermomicrobium roseum
#sti	Sphaerobacter thermophilus
#atm	Anaerolinea thermophila
#dra	Deinococcus radiodurans
#dge	Deinococcus geothermalis
#ddr	Deinococcus deserti
#dmr	Deinococcus maricopensis
#dpt	Deinococcus proteolyticus
#tra	Truepera radiovictrix
#tth	Thermus thermophilus HB27
#ttj	Thermus thermophilus HB8
#tsc	Thermus scotoductus
#mrb	Meiothermus ruber
#msv	Meiothermus silvanus
#opr	Oceanithermus profundus
#mhd	Marinithermus hydrothermalis
#aae	Aquifex aeolicus
#hya	Hydrogenobaculum sp. Y04AAS1
#hyd	Hydrogenobaculum sp. 3684
#hys	Hydrogenobaculum sp. SHO
#hth	Hydrogenobacter thermophilus
#tal	Thermocrinis albus
#sul	Sulfurihydrogenibium sp. YO3AOP1
#saf	Sulfurihydrogenibium azorense
#pmx	Persephonella marina
#tam	Thermovibrio ammonificans
#dte	Desulfurobacterium thermolithotrophum
#tma	Thermotoga maritima
#tpt	Thermotoga petrophila
#tle	Thermotoga lettingae
#trq	Thermotoga sp. RQ2
#tna	Thermotoga neapolitana
#tnp	Thermotoga naphthophila
#tta	Thermotoga thermarum
#tme	Thermosipho melanesiensis
#taf	Thermosipho africanus
#fno	Fervidobacterium nodosum
#pmo	Petrotoga mobilis
#kol	Kosmotoga olearia
#dth	Dictyoglomus thermophilum
#dtu	Dictyoglomus turgidum
#tye	Thermodesulfovibrio yellowstonii
#nde	Candidatus Nitrospira defluvii
#ttr	Thermobaculum terrenum
#ddf	Deferribacter desulfuricans SSM1
#dap	Denitrovibrio acetiphilus
#cni	Calditerrivibrio nitroreducens
#fsi	Flexistipes sinusarabici
#tid	Thermodesulfatator indicus
#top	Thermodesulfobacterium sp. OPB45
#mja	Methanocaldococcus jannaschii
#mfe	Methanocaldococcus fervens
#mvu	Methanocaldococcus vulcanius
#mfs	Methanocaldococcus sp. FS406-22
#mif	Methanocaldococcus infernus
#mig	Methanotorris igneus
#mmp	Methanococcus maripaludis S2
#mmq	Methanococcus maripaludis C5
#mmx	Methanococcus maripaludis C6
#mmz	Methanococcus maripaludis C7
#mmd	Methanococcus maripaludis XI
#mae	Methanococcus aeolicus
#mvn	Methanococcus vannielii
#mvo	Methanococcus voltae
#mok	Methanothermococcus okinawensis
#mac	Methanosarcina acetivorans
#mba	Methanosarcina barkeri
#mma	Methanosarcina mazei
#mbu	Methanococcoides burtonii
#mmh	Methanohalophilus mahii
#mev	Methanohalobium evestigatum
#mzh	Methanosalsum zhilinae
#mtp	Methanosaeta thermophila
#mcj	Methanosaeta concilii
#mhu	Methanospirillum hungatei
#mla	Methanocorpusculum labreanum
#mem	Methanoculleus marisnigri
#mpi	Methanoplanus petrolearius
#mbn	Candidatus Methanoregula boonei
#mpl	Candidatus Methanosphaerula palustris
#mpd	Methanocella paludicola
#mth	Methanothermobacter thermautotrophicus
#mmg	Methanothermobacter marburgensis
#mst	Methanosphaera stadtmanae
#msi	Methanobrevibacter smithii ATCC 35061
#mru	Methanobrevibacter ruminantium
#mel	Methanobacterium sp. AL-21
#mew	Methanobacterium sp. SWAN-1
#mfv	Methanothermus fervidus
#mka	Methanopyrus kandleri
#afu	Archaeoglobus fulgidus
#apo	Archaeoglobus profundus
#ave	Archaeoglobus veneficus
#fpl	Ferroglobus placidus
#hal	Halobacterium sp. NRC-1
#hsl	Halobacterium salinarum R1
#hma	Haloarcula marismortui
#hhi	Haloarcula hispanica
#hwa	Haloquadratum walsbyi
#nph	Natronomonas pharaonis
#hla	Halorubrum lacusprofundi
#hut	Halorhabdus utahensis
#hmu	Halomicrobium mukohataei
#htu	Haloterrigena turkmenica
#nmg	Natrialba magadii
#hvo	Haloferax volcanii
#hje	Halalkalicoccus jeotgali
#hbo	Halogeometricum borinquense
#hxa	Halopiger xanaduensis
#tac	Thermoplasma acidophilum
#tvo	Thermoplasma volcanium
#pto	Picrophilus torridus
#pho	Pyrococcus horikoshii
#pab	Pyrococcus abyssi
#pfu	Pyrococcus furiosus
#pyn	Pyrococcus sp. NA2
#pya	Pyrococcus yayanosii
#tko	Thermococcus kodakaraensis
#ton	Thermococcus onnurineus
#tga	Thermococcus gammatolerans
#tsi	Thermococcus sibiricus
#tba	Thermococcus barophilus
#the	Thermococcus sp. 4557
#abi	Aciduliprofundum boonei
#rci	Uncultured methanogenic archaeon RC-I
#ape	Aeropyrum pernix
#smr	Staphylothermus marinus
#shc	Staphylothermus hellenicus
#iho	Ignicoccus hospitalis
#dka	Desulfurococcus kamchatkensis
#dmu	Desulfurococcus mucosus
#tag	Thermosphaera aggregans
#iag	Ignisphaera aggregans
#hbu	Hyperthermus butylicus
#sso	Sulfolobus solfataricus
#sto	Sulfolobus tokodaii
#sai	Sulfolobus acidocaldarius
#sis	Sulfolobus islandicus L.S.2.15
#sia	Sulfolobus islandicus M.14.25
#sim	Sulfolobus islandicus M.16.27
#sid	Sulfolobus islandicus M.16.4
#siy	Sulfolobus islandicus Y.G.57.14
#sin	Sulfolobus islandicus Y.N.15.51
#sii	Sulfolobus islandicus L.D.8.5
#mse	Metallosphaera sedula
#mcn	Metallosphaera cuprina
#aho	Acidianus hospitalis
#pai	Pyrobaculum aerophilum
#pis	Pyrobaculum islandicum
#pcl	Pyrobaculum calidifontis
#pas	Pyrobaculum arsenaticum
#cma	Caldivirga maquilingensis
#tne	Thermoproteus neutrophilus
#tuz	Thermoproteus uzoniensis
#vdi	Vulcanisaeta distributa
#vmo	Vulcanisaeta moutnovskia
#tpe	Thermofilum pendens
#asc	Acidilobus saccharovorans
#nmr	Nitrosopumilus maritimus
#csy	Cenarchaeum symbiosum A
#neq	Nanoarchaeum equitans
#kcr	Candidatus Korarchaeum cryptofilum


#my $serv = SOAP::Lite->service($wsdl);
#
#my $offset = 1;
#my $limit = 5;
#
#my $top5 = $serv->get_best_neighbors_by_gene('eco:b0002', $offset, $limit);
#
#foreach my $hit (@{$top5}) {
#  print "$hit->{genes_id1}\t$hit->{genes_id2}\t$hit->{sw_score}\n";
#}

#my $paths = SOAP::Lite
#             -> service($wsdl)
#             -> list_pathways("hsa");
#foreach my $path (@{$paths}) {
#  print "$path->{entry_id}\t$path->{definition}\n";
#}
#path:hsa00010	Glycolysis / Gluconeogenesis - Homo sapiens (human)
#path:hsa00020	Citrate cycle (TCA cycle) - Homo sapiens (human)
#path:hsa00030	Pentose phosphate pathway - Homo sapiens (human)
#path:hsa00040	Pentose and glucuronate interconversions - Homo sapiens (human)
#path:hsa00051	Fructose and mannose metabolism - Homo sapiens (human)
#path:hsa00052	Galactose metabolism - Homo sapiens (human)
#path:hsa00053	Ascorbate and aldarate metabolism - Homo sapiens (human)
#path:hsa00061	Fatty acid biosynthesis - Homo sapiens (human)
#path:hsa00062	Fatty acid elongation - Homo sapiens (human)
#path:hsa00071	Fatty acid metabolism - Homo sapiens (human)
#path:hsa00072	Synthesis and degradation of ketone bodies - Homo sapiens (human)
#path:hsa00100	Steroid biosynthesis - Homo sapiens (human)
#path:hsa00120	Primary bile acid biosynthesis - Homo sapiens (human)
#path:hsa00130	Ubiquinone and other terpenoid-quinone biosynthesis - Homo sapiens (human)
#path:hsa00140	Steroid hormone biosynthesis - Homo sapiens (human)
#path:hsa00190	Oxidative phosphorylation - Homo sapiens (human)
#path:hsa00230	Purine metabolism - Homo sapiens (human)
#path:hsa00232	Caffeine metabolism - Homo sapiens (human)
#path:hsa00240	Pyrimidine metabolism - Homo sapiens (human)
#path:hsa00250	Alanine, aspartate and glutamate metabolism - Homo sapiens (human)
#path:hsa00260	Glycine, serine and threonine metabolism - Homo sapiens (human)
#path:hsa00270	Cysteine and methionine metabolism - Homo sapiens (human)
#path:hsa00280	Valine, leucine and isoleucine degradation - Homo sapiens (human)
#path:hsa00290	Valine, leucine and isoleucine biosynthesis - Homo sapiens (human)
#path:hsa00300	Lysine biosynthesis - Homo sapiens (human)
#path:hsa00310	Lysine degradation - Homo sapiens (human)
#path:hsa00330	Arginine and proline metabolism - Homo sapiens (human)
#path:hsa00340	Histidine metabolism - Homo sapiens (human)
#path:hsa00350	Tyrosine metabolism - Homo sapiens (human)
#path:hsa00360	Phenylalanine metabolism - Homo sapiens (human)
#path:hsa00380	Tryptophan metabolism - Homo sapiens (human)
#path:hsa00400	Phenylalanine, tyrosine and tryptophan biosynthesis - Homo sapiens (human)
#path:hsa00410	beta-Alanine metabolism - Homo sapiens (human)
#path:hsa00430	Taurine and hypotaurine metabolism - Homo sapiens (human)
#path:hsa00450	Selenocompound metabolism - Homo sapiens (human)
#path:hsa00460	Cyanoamino acid metabolism - Homo sapiens (human)
#path:hsa00471	D-Glutamine and D-glutamate metabolism - Homo sapiens (human)
#path:hsa00472	D-Arginine and D-ornithine metabolism - Homo sapiens (human)
#path:hsa00480	Glutathione metabolism - Homo sapiens (human)
#path:hsa00500	Starch and sucrose metabolism - Homo sapiens (human)
#path:hsa00510	N-Glycan biosynthesis - Homo sapiens (human)
#path:hsa00511	Other glycan degradation - Homo sapiens (human)
#path:hsa00512	Mucin type O-Glycan biosynthesis - Homo sapiens (human)
#path:hsa00514	Other types of O-glycan biosynthesis - Homo sapiens (human)
#path:hsa00520	Amino sugar and nucleotide sugar metabolism - Homo sapiens (human)
#path:hsa00524	Butirosin and neomycin biosynthesis - Homo sapiens (human)
#path:hsa00531	Glycosaminoglycan degradation - Homo sapiens (human)
#path:hsa00532	Glycosaminoglycan biosynthesis - chondroitin sulfate - Homo sapiens (human)
#path:hsa00533	Glycosaminoglycan biosynthesis - keratan sulfate - Homo sapiens (human)
#path:hsa00534	Glycosaminoglycan biosynthesis - heparan sulfate - Homo sapiens (human)
#path:hsa00561	Glycerolipid metabolism - Homo sapiens (human)
#path:hsa00562	Inositol phosphate metabolism - Homo sapiens (human)
#path:hsa00563	Glycosylphosphatidylinositol(GPI)-anchor biosynthesis - Homo sapiens (human)
#path:hsa00564	Glycerophospholipid metabolism - Homo sapiens (human)
#path:hsa00565	Ether lipid metabolism - Homo sapiens (human)
#path:hsa00590	Arachidonic acid metabolism - Homo sapiens (human)
#path:hsa00591	Linoleic acid metabolism - Homo sapiens (human)
#path:hsa00592	alpha-Linolenic acid metabolism - Homo sapiens (human)
#path:hsa00600	Sphingolipid metabolism - Homo sapiens (human)
#path:hsa00601	Glycosphingolipid biosynthesis - lacto and neolacto series - Homo sapiens (human)
#path:hsa00603	Glycosphingolipid biosynthesis - globo series - Homo sapiens (human)
#path:hsa00604	Glycosphingolipid biosynthesis - ganglio series - Homo sapiens (human)
#path:hsa00620	Pyruvate metabolism - Homo sapiens (human)
#path:hsa00630	Glyoxylate and dicarboxylate metabolism - Homo sapiens (human)
#path:hsa00640	Propanoate metabolism - Homo sapiens (human)
#path:hsa00650	Butanoate metabolism - Homo sapiens (human)
#path:hsa00670	One carbon pool by folate - Homo sapiens (human)
#path:hsa00730	Thiamine metabolism - Homo sapiens (human)
#path:hsa00740	Riboflavin metabolism - Homo sapiens (human)
#path:hsa00750	Vitamin B6 metabolism - Homo sapiens (human)
#path:hsa00760	Nicotinate and nicotinamide metabolism - Homo sapiens (human)
#path:hsa00770	Pantothenate and CoA biosynthesis - Homo sapiens (human)
#path:hsa00780	Biotin metabolism - Homo sapiens (human)
#path:hsa00785	Lipoic acid metabolism - Homo sapiens (human)
#path:hsa00790	Folate biosynthesis - Homo sapiens (human)
#path:hsa00830	Retinol metabolism - Homo sapiens (human)
#path:hsa00860	Porphyrin and chlorophyll metabolism - Homo sapiens (human)
#path:hsa00900	Terpenoid backbone biosynthesis - Homo sapiens (human)
#path:hsa00910	Nitrogen metabolism - Homo sapiens (human)
#path:hsa00920	Sulfur metabolism - Homo sapiens (human)
#path:hsa00970	Aminoacyl-tRNA biosynthesis - Homo sapiens (human)
#path:hsa00980	Metabolism of xenobiotics by cytochrome P450 - Homo sapiens (human)
#path:hsa00982	Drug metabolism - cytochrome P450 - Homo sapiens (human)
#path:hsa00983	Drug metabolism - other enzymes - Homo sapiens (human)
#path:hsa01040	Biosynthesis of unsaturated fatty acids - Homo sapiens (human)
#path:hsa01100	Metabolic pathways - Homo sapiens (human)
#path:hsa02010	ABC transporters - Homo sapiens (human)
#path:hsa03008	Ribosome biogenesis in eukaryotes - Homo sapiens (human)
#path:hsa03010	Ribosome - Homo sapiens (human)
#path:hsa03013	RNA transport - Homo sapiens (human)
#path:hsa03015	mRNA surveillance pathway - Homo sapiens (human)
#path:hsa03018	RNA degradation - Homo sapiens (human)
#path:hsa03020	RNA polymerase - Homo sapiens (human)
#path:hsa03022	Basal transcription factors - Homo sapiens (human)
#path:hsa03030	DNA replication - Homo sapiens (human)
#path:hsa03040	Spliceosome - Homo sapiens (human)
#path:hsa03050	Proteasome - Homo sapiens (human)
#path:hsa03060	Protein export - Homo sapiens (human)
#path:hsa03320	PPAR signaling pathway - Homo sapiens (human)
#path:hsa03410	Base excision repair - Homo sapiens (human)
#path:hsa03420	Nucleotide excision repair - Homo sapiens (human)
#path:hsa03430	Mismatch repair - Homo sapiens (human)
#path:hsa03440	Homologous recombination - Homo sapiens (human)
#path:hsa03450	Non-homologous end-joining - Homo sapiens (human)
#path:hsa03460	Fanconi anemia pathway - Homo sapiens (human)
#path:hsa04010	MAPK signaling pathway - Homo sapiens (human)
#path:hsa04012	ErbB signaling pathway - Homo sapiens (human)
#path:hsa04020	Calcium signaling pathway - Homo sapiens (human)
#path:hsa04060	Cytokine-cytokine receptor interaction - Homo sapiens (human)
#path:hsa04062	Chemokine signaling pathway - Homo sapiens (human)
#path:hsa04070	Phosphatidylinositol signaling system - Homo sapiens (human)
#path:hsa04080	Neuroactive ligand-receptor interaction - Homo sapiens (human)
#path:hsa04110	Cell cycle - Homo sapiens (human)
#path:hsa04114	Oocyte meiosis - Homo sapiens (human)
#path:hsa04115	p53 signaling pathway - Homo sapiens (human)
#path:hsa04120	Ubiquitin mediated proteolysis - Homo sapiens (human)
#path:hsa04122	Sulfur relay system - Homo sapiens (human)
#path:hsa04130	SNARE interactions in vesicular transport - Homo sapiens (human)
#path:hsa04140	Regulation of autophagy - Homo sapiens (human)
#path:hsa04141	Protein processing in endoplasmic reticulum - Homo sapiens (human)
#path:hsa04142	Lysosome - Homo sapiens (human)
#path:hsa04144	Endocytosis - Homo sapiens (human)
#path:hsa04145	Phagosome - Homo sapiens (human)
#path:hsa04146	Peroxisome - Homo sapiens (human)
#path:hsa04150	mTOR signaling pathway - Homo sapiens (human)
#path:hsa04210	Apoptosis - Homo sapiens (human)
#path:hsa04260	Cardiac muscle contraction - Homo sapiens (human)
#path:hsa04270	Vascular smooth muscle contraction - Homo sapiens (human)
#path:hsa04310	Wnt signaling pathway - Homo sapiens (human)
#path:hsa04320	Dorso-ventral axis formation - Homo sapiens (human)
#path:hsa04330	Notch signaling pathway - Homo sapiens (human)
#path:hsa04340	Hedgehog signaling pathway - Homo sapiens (human)
#path:hsa04350	TGF-beta signaling pathway - Homo sapiens (human)
#path:hsa04360	Axon guidance - Homo sapiens (human)
#path:hsa04370	VEGF signaling pathway - Homo sapiens (human)
#path:hsa04380	Osteoclast differentiation - Homo sapiens (human)
#path:hsa04510	Focal adhesion - Homo sapiens (human)
#path:hsa04512	ECM-receptor interaction - Homo sapiens (human)
#path:hsa04514	Cell adhesion molecules (CAMs) - Homo sapiens (human)
#path:hsa04520	Adherens junction - Homo sapiens (human)
#path:hsa04530	Tight junction - Homo sapiens (human)
#path:hsa04540	Gap junction - Homo sapiens (human)
#path:hsa04610	Complement and coagulation cascades - Homo sapiens (human)
#path:hsa04612	Antigen processing and presentation - Homo sapiens (human)
#path:hsa04614	Renin-angiotensin system - Homo sapiens (human)
#path:hsa04620	Toll-like receptor signaling pathway - Homo sapiens (human)
#path:hsa04621	NOD-like receptor signaling pathway - Homo sapiens (human)
#path:hsa04622	RIG-I-like receptor signaling pathway - Homo sapiens (human)
#path:hsa04623	Cytosolic DNA-sensing pathway - Homo sapiens (human)
#path:hsa04630	Jak-STAT signaling pathway - Homo sapiens (human)
#path:hsa04640	Hematopoietic cell lineage - Homo sapiens (human)
#path:hsa04650	Natural killer cell mediated cytotoxicity - Homo sapiens (human)
#path:hsa04660	T cell receptor signaling pathway - Homo sapiens (human)
#path:hsa04662	B cell receptor signaling pathway - Homo sapiens (human)
#path:hsa04664	Fc epsilon RI signaling pathway - Homo sapiens (human)
#path:hsa04666	Fc gamma R-mediated phagocytosis - Homo sapiens (human)
#path:hsa04670	Leukocyte transendothelial migration - Homo sapiens (human)
#path:hsa04672	Intestinal immune network for IgA production - Homo sapiens (human)
#path:hsa04710	Circadian rhythm - mammal - Homo sapiens (human)
#path:hsa04720	Long-term potentiation - Homo sapiens (human)
#path:hsa04721	Synaptic vesicle cycle - Homo sapiens (human)
#path:hsa04722	Neurotrophin signaling pathway - Homo sapiens (human)
#path:hsa04724	Glutamatergic synapse - Homo sapiens (human)
#path:hsa04725	Cholinergic synapse - Homo sapiens (human)
#path:hsa04727	GABAergic synapse - Homo sapiens (human)
#path:hsa04730	Long-term depression - Homo sapiens (human)
#path:hsa04740	Olfactory transduction - Homo sapiens (human)
#path:hsa04742	Taste transduction - Homo sapiens (human)
#path:hsa04744	Phototransduction - Homo sapiens (human)
#path:hsa04810	Regulation of actin cytoskeleton - Homo sapiens (human)
#path:hsa04910	Insulin signaling pathway - Homo sapiens (human)
#path:hsa04912	GnRH signaling pathway - Homo sapiens (human)
#path:hsa04914	Progesterone-mediated oocyte maturation - Homo sapiens (human)
#path:hsa04916	Melanogenesis - Homo sapiens (human)
#path:hsa04920	Adipocytokine signaling pathway - Homo sapiens (human)
#path:hsa04930	Type II diabetes mellitus - Homo sapiens (human)
#path:hsa04940	Type I diabetes mellitus - Homo sapiens (human)
#path:hsa04950	Maturity onset diabetes of the young - Homo sapiens (human)
#path:hsa04960	Aldosterone-regulated sodium reabsorption - Homo sapiens (human)
#path:hsa04961	Endocrine and other factor-regulated calcium reabsorption - Homo sapiens (human)
#path:hsa04962	Vasopressin-regulated water reabsorption - Homo sapiens (human)
#path:hsa04964	Proximal tubule bicarbonate reclamation - Homo sapiens (human)
#path:hsa04966	Collecting duct acid secretion - Homo sapiens (human)
#path:hsa04970	Salivary secretion - Homo sapiens (human)
#path:hsa04971	Gastric acid secretion - Homo sapiens (human)
#path:hsa04972	Pancreatic secretion - Homo sapiens (human)
#path:hsa04973	Carbohydrate digestion and absorption - Homo sapiens (human)
#path:hsa04974	Protein digestion and absorption - Homo sapiens (human)
#path:hsa04975	Fat digestion and absorption - Homo sapiens (human)
#path:hsa04976	Bile secretion - Homo sapiens (human)
#path:hsa04977	Vitamin digestion and absorption - Homo sapiens (human)
#path:hsa04978	Mineral absorption - Homo sapiens (human)
#path:hsa05010	Alzheimer's disease - Homo sapiens (human)
#path:hsa05012	Parkinson's disease - Homo sapiens (human)
#path:hsa05014	Amyotrophic lateral sclerosis (ALS) - Homo sapiens (human)
#path:hsa05016	Huntington's disease - Homo sapiens (human)
#path:hsa05020	Prion diseases - Homo sapiens (human)
#path:hsa05100	Bacterial invasion of epithelial cells - Homo sapiens (human)
#path:hsa05110	Vibrio cholerae infection - Homo sapiens (human)
#path:hsa05120	Epithelial cell signaling in Helicobacter pylori infection - Homo sapiens (human)
#path:hsa05130	Pathogenic Escherichia coli infection - Homo sapiens (human)
#path:hsa05131	Shigellosis - Homo sapiens (human)
#path:hsa05132	Salmonella infection - Homo sapiens (human)
#path:hsa05133	Pertussis - Homo sapiens (human)
#path:hsa05140	Leishmaniasis - Homo sapiens (human)
#path:hsa05142	Chagas disease (American trypanosomiasis) - Homo sapiens (human)
#path:hsa05143	African trypanosomiasis - Homo sapiens (human)
#path:hsa05144	Malaria - Homo sapiens (human)
#path:hsa05145	Toxoplasmosis - Homo sapiens (human)
#path:hsa05146	Amoebiasis - Homo sapiens (human)
#path:hsa05150	Staphylococcus aureus infection - Homo sapiens (human)
#path:hsa05152	Tuberculosis - Homo sapiens (human)
#path:hsa05160	Hepatitis C - Homo sapiens (human)
#path:hsa05162	Measles - Homo sapiens (human)
#path:hsa05164	Influenza A - Homo sapiens (human)
#path:hsa05166	HTLV-I infection - Homo sapiens (human)
#path:hsa05168	Herpes simplex infection - Homo sapiens (human)
#path:hsa05200	Pathways in cancer - Homo sapiens (human)
#path:hsa05210	Colorectal cancer - Homo sapiens (human)
#path:hsa05211	Renal cell carcinoma - Homo sapiens (human)
#path:hsa05212	Pancreatic cancer - Homo sapiens (human)
#path:hsa05213	Endometrial cancer - Homo sapiens (human)
#path:hsa05214	Glioma - Homo sapiens (human)
#path:hsa05215	Prostate cancer - Homo sapiens (human)
#path:hsa05216	Thyroid cancer - Homo sapiens (human)
#path:hsa05217	Basal cell carcinoma - Homo sapiens (human)
#path:hsa05218	Melanoma - Homo sapiens (human)
#path:hsa05219	Bladder cancer - Homo sapiens (human)
#path:hsa05220	Chronic myeloid leukemia - Homo sapiens (human)
#path:hsa05221	Acute myeloid leukemia - Homo sapiens (human)
#path:hsa05222	Small cell lung cancer - Homo sapiens (human)
#path:hsa05223	Non-small cell lung cancer - Homo sapiens (human)
#path:hsa05310	Asthma - Homo sapiens (human)
#path:hsa05320	Autoimmune thyroid disease - Homo sapiens (human)
#path:hsa05322	Systemic lupus erythematosus - Homo sapiens (human)
#path:hsa05323	Rheumatoid arthritis - Homo sapiens (human)
#path:hsa05330	Allograft rejection - Homo sapiens (human)
#path:hsa05332	Graft-versus-host disease - Homo sapiens (human)
#path:hsa05340	Primary immunodeficiency - Homo sapiens (human)
#path:hsa05410	Hypertrophic cardiomyopathy (HCM) - Homo sapiens (human)
#path:hsa05412	Arrhythmogenic right ventricular cardiomyopathy (ARVC) - Homo sapiens (human)
#path:hsa05414	Dilated cardiomyopathy - Homo sapiens (human)
#path:hsa05416	Viral myocarditis - Homo sapiens (human)



#my $genes = SOAP::Data->type(array => ["eco:b1002", "eco:b2388"]);

#my $result = $serv -> mark_pathway_by_objects("path:eco00010", $genes);
#
#print $result;	# URL of the generated image